%% MULTI-ASPECT + MULTI-PASS TOMOGRAPHY — บิน 4 ทิศเพื่อให้เห็นกล่องเต็มใบ
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  ต่อจาก box_tomography_multipass.m : ปัญหาคือ "เห็นแค่บางด้าน" — บินเส้นเดียวฝั่งเดียว จึง
%    - Front เห็นชัด (broadside)
%    - Left/Right เห็นจาง (มองเฉียงมากจากปลาย track)
%    - Back ถูกบังสนิท (occlusion จริง)
%
%  v3 แก้ด้วย MULTI-ASPECT: บิน track เดิม (multi-pass 21 ชั้น) ซ้ำ 4 ทิศ
%  รอบกล่อง (0/90/180/270 องศา) -> ทุกผนังได้ broadside ของตัวเอง
%    - ในแต่ละ aspect: รวม coherent (tomography -> height โฟกัส)
%    - ข้าม aspect:    รวม incoherent (max ของ |image|) เพราะมุมต่างกัน 90
%                      องศา เฟส/สเปกเคิลไม่ coherent กัน (ทำแบบ SAR จริง)
%  ผล: เห็นผนัง 4 ด้าน + หลังคา = กล่องเต็มใบ
%  (Bottom ยังไม่มีทางเห็น — กล่องวางบนพื้น เรดาร์มองจากบน = ถูกฟิสิกส์)
%
%  RUNTIME: ~4 เท่าของ v2 (4 aspects) แต่ BP เขียนใหม่ให้เร็วขึ้น
%  (manual linear interp แทน interp1) ลด aspects ได้ที่ตัวแปร aspectDeg

clear; clc; close all;

%% ===== RADAR / GEOMETRY PARAMETERS =====
c   = physconst('LightSpeed');
fc  = 4e9;      lambda = c/fc;
fs  = 120e6;
prf = 1000;
speed = 100;   flightDuration = 4;
rngres = 3;                    % m

baseH  = 500;
groundR0 = 1000;               % ระยะ ground range จาก track ถึงกล่อง (ทุก aspect)

%% ===== MULTI-PASS + MULTI-ASPECT CONFIG =====
Npass  = 21;
dzPass = 1.0;
azDecim = 5;
elevWindow = true;
aspectDeg = [0 90 180 270];    % ทิศการบินรอบกล่อง (ลดเหลือ [0] = v2 เดิม)
Nasp = numel(aspectDeg);

Naz_full = flightDuration*prf + 1;
azIdx    = 1:azDecim:Naz_full;
Naz      = numel(azIdx);
t        = (azIdx-1)/prf;
ry       = -200 + speed*t;     % along-track offset

%% ===== BOX SCATTERERS (สร้างครบ 6 หน้า แล้วค่อย mask ราย aspect) =====
boxCenter = [groundR0; 0; 0];
boxL = 30; boxW = 20; boxH = 10;
Nf = 10;
surfaceModel = 'rough';

[scatPosAll, scatAmpAll, scatFace] = buildBoxScatterers( ...
    boxCenter, boxL, boxW, boxH, Nf, surfaceModel);
faceNormal = [ -1 0 0;  1 0 0;  0 -1 0;  0 1 0;  0 0 1;  0 0 -1];
faceLabels = {'Front','Back','Left','Right','Top','Bottom'};
fprintf('Surface model: %s | scatterers (all faces): %d\n', surfaceModel, size(scatPosAll,2));

%% ===== VOXEL GRID =====
voxStep = 1.0;
xScene = (boxCenter(2)-30) : voxStep : (boxCenter(2)+30);   % cross
yScene = (boxCenter(1)-20) : voxStep : (boxCenter(1)+20);   % range
zScene = -5               : voxStep : 25;                   % height
[Xg, Yg, Zg] = meshgrid(xScene, yScene, zScene);
Ny = numel(yScene); Nx = numel(xScene); Nz = numel(zScene);
sz   = size(Xg);
voxT = [Yg(:), Xg(:), Zg(:)];             % Nvox x 3  [range x, cross y, z]
fprintf('Voxel grid: %d x %d x %d = %d voxels\n', Ny,Nx,Nz,Ny*Nx*Nz);

rbins = (1080 : c/2/fs : 1180).';

%% ===== MULTI-ASPECT TOMOGRAPHIC BP =====
%  aspect 0 = เรดาร์ฝั่ง -x มองไปทางกล่อง (เหมือน v2), หมุน track รอบกล่อง
if elevWindow, wElev = hamming(Npass); else, wElev = ones(Npass,1); end
magComb   = zeros(sz);         % รวม incoherent ข้าม aspect (max)
magAspect1 = [];               % เก็บ aspect แรกไว้เทียบ (= single aspect)
tAll = tic;

for iA = 1:Nasp
    phi = aspectDeg(iA);
    R2  = [cosd(phi) -sind(phi); sind(phi) cosd(phi)];

    % track relative to box (aspect 0): [-groundR0 ; ry], หมุนด้วย R2
    relXY = R2 * [-groundR0*ones(1,Naz); ry];
    platXY = boxCenter(1:2) + relXY;            % 2 x Naz (คงที่ทุก pass)

    % occlusion สำหรับ aspect นี้ (ใช้ pass กลาง)
    Hm = baseH + (Npass-1)/2*dzPass;
    losA = [platXY; Hm*ones(1,Naz)] - boxCenter;
    visF = false(1,6);
    for f = 1:6, visF(f) = any(faceNormal(f,:)*losA > 0); end
    ampA = scatAmpAll .* visF(scatFace);        % ปิดหน้าที่มองไม่เห็น
    keep = ampA ~= 0;
    posK = scatPosAll(:,keep);  ampK = ampA(keep);
    fprintf('Aspect %3d deg | visible: %s | scatterers: %d\n', ...
        phi, strjoin(faceLabels(visF),', '), nnz(keep));

    % tomographic BP (coherent ข้าม pass ภายใน aspect)
    bpA = zeros(sz);
    for p = 1:Npass
        Hp   = baseH + (p-1)*dzPass;
        plat = [platXY; Hp*ones(1,Naz)];
        S    = simPassAnalytic(posK, ampK, plat, rbins, rngres, lambda);
        bpA  = bpA + wElev(p) * bpOnePass(S, plat, rbins, voxT, fc, c, sz);
        fprintf('  aspect %d/%d  pass %2d/%2d\n', iA, Nasp, p, Npass);
    end
    magA = abs(bpA);
    if iA == 1, magAspect1 = magA; end
    magComb = max(magComb, magA);               % incoherent fusion
end
fprintf('Multi-aspect BP complete (%.1f s total).\n', toc(tAll));

magComb    = magComb    ./ max(magComb(:));
magAspect1 = magAspect1 ./ max(magAspect1(:));

%% ===== POST-PROCESSING FOR DISPLAY =====
magSm = smooth3(magComb, 'gaussian', [3 3 5], 1.2);
magSm = magSm ./ max(magSm(:));
magSmdB = 20*log10(magSm + eps);

[Vb,Fb] = makeBoxWire(boxCenter, boxL/2, boxW/2, boxH);
[Xe,Ye,Ze] = boxEdges(boxCenter, boxL/2, boxW/2, boxH);   % กรอบเส้น (หมุนไม่บั๊ก)
boxCr = boxCenter(2)+boxW/2*[-1 1 1 -1 -1];
boxRg = boxCenter(1)+boxL/2*[-1 -1 1 1 -1];
xlims = [-18 18]; ylims = [982 1018]; zlims = [-3 15];

%% ===== FIGURE 1: ACQUISITION GEOMETRY — 4 aspects รอบกล่อง =====
figure(1); set(gcf,'Name','Fig1: multi-aspect geometry','Color','w','Position',[60 80 900 640]);
aspCol = lines(Nasp);
patch('Vertices',Vb,'Faces',Fb,'FaceColor',[0.6 0.8 1],'FaceAlpha',0.3, ...
      'EdgeColor',[0.1 0.3 0.9],'LineWidth',1.5); hold on;
for iA = 1:Nasp
    phi = aspectDeg(iA);
    R2  = [cosd(phi) -sind(phi); sind(phi) cosd(phi)];
    relXY = R2 * [-groundR0*ones(1,Naz); ry];
    platXY = boxCenter(1:2) + relXY;
    for p = 1:5:Npass                    % วาดบางส่วนพอเห็น stack
        Hp = baseH + (p-1)*dzPass;
        plot3(platXY(2,:), platXY(1,:), Hp*ones(1,Naz), '-', ...
              'Color', aspCol(iA,:), 'LineWidth', 1);
    end
    text(platXY(2,round(Naz/2)), platXY(1,round(Naz/2)), baseH+60, ...
        sprintf('aspect %d\\circ', phi), 'Color', aspCol(iA,:), 'FontWeight','bold');
end
hold off; grid on; view(-40,25); daspect([1 1 0.5]);
xlabel('Cross / y (m)'); ylabel('Range / x (m)'); zlabel('Height (m)');
title(sprintf('Fig 1 | W8 v3 — %d aspects x %d passes around the box', Nasp, Npass));

%% ===== FIGURE 2: SINGLE ASPECT vs MULTI-ASPECT (top-down MIP) =====
mipA1 = squeeze(max(magAspect1,[],3));
mipMC = squeeze(max(magComb,[],3));

figure(2); set(gcf,'Name','Fig2: single vs multi aspect','Color','k','Position',[80 60 950 460]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

axa = nexttile;
imagesc(xScene, yScene, 20*log10(mipA1+eps)); set(axa,'YDir','normal'); clim([-30 0]); colormap(axa,'gray');
hold on; plot(boxCr,boxRg,'c--','LineWidth',1.5); hold off;
xlabel('Cross (m)','Color','w'); ylabel('Range (m)','Color','w');
title('1 aspect — front wall only','Color','w');
set(axa,'Color','k','XColor','w','YColor','w');

axb = nexttile;
imagesc(xScene, yScene, 20*log10(mipMC+eps)); set(axb,'YDir','normal'); clim([-30 0]); colormap(axb,'gray');
hold on; plot(boxCr,boxRg,'c--','LineWidth',1.5); hold off;
xlabel('Cross (m)','Color','w'); ylabel('Range (m)','Color','w');
title(sprintf('%d aspects — all 4 walls', Nasp),'Color','w');
set(axb,'Color','k','XColor','w','YColor','w');
cb = colorbar; cb.Color='w'; ylabel(cb,'dB','Color','w'); cb.Layout.Tile='east';
sgtitle('Fig 2 | W8 v3 — multi-aspect fills in the missing walls (top-down MIP)','Color','w');

%% ===== FIGURE 3: isosurface 4 มุมมอง — กล่องเต็มใบ =====
figure(3); set(gcf,'Name','Fig3: full-box isosurface (v3)','Color','w','Position',[60 60 1150 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

isoLv  = [0.4, 0.15];
isoCol = [0.85 0.25 0.15; 1.00 0.65 0.20];
isoAl  = [0.9, 0.22];
viewAz = [-35, 145,  90, -35];
viewEl = [ 22,  22,   0,  85];
viewNm = {'Perspective A','Perspective B (back)','Front view','Top-down'};

for k = 1:4
    ax = nexttile; hold(ax,'on');
    for L = 1:2
        p = patch(ax, isosurface(Xg, Yg, Zg, magSm, isoLv(L)));
        isonormals(Xg, Yg, Zg, magSm, p);
        set(p,'FaceColor',isoCol(L,:),'EdgeColor','none','FaceAlpha',isoAl(L));
    end
    plot3(ax, Xe, Ye, Ze, '-', 'Color',[0.1 0.3 0.9], 'LineWidth',1.5);   % กรอบเส้น
    hold(ax,'off');
    xlabel('Cross (m)'); ylabel('Range (m)'); zlabel('Height (m)');
    title(viewNm{k});
    xlim(xlims); ylim(ylims); zlim(zlims);
    daspect(ax,[1 1 1]); grid(ax,'on'); box(ax,'on');
    view(ax, viewAz(k), viewEl(k));
    camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
end
sgtitle(sprintf('Fig 3 | W8 v3 — full box: 4 walls + roof recovered (isosurface @ %.0f/%.0f dB)', ...
    20*log10(isoLv(1)), 20*log10(isoLv(2))));
rotate3d on;

%% ===== FIGURE 4: point cloud =====
pcThr = -18;
idxPC = find(magSmdB >= pcThr);
figure(4); set(gcf,'Name','Fig4: voxel point cloud (v3)','Color','w','Position',[120 80 780 640]);
scatter3(Xg(idxPC), Yg(idxPC), Zg(idxPC), 12, magSmdB(idxPC), 'filled', ...
    'MarkerFaceAlpha', 0.4);
hold on;
plot3(Xe, Ye, Ze, '-', 'Color',[0.1 0.3 0.9], 'LineWidth',1.8);   % กรอบเส้น
hold off;
colormap(jet); clim([pcThr 0]);
cb4 = colorbar; ylabel(cb4,'dB');
xlabel('Cross (m)'); ylabel('Range (m)'); zlabel('Height (m)');
xlim(xlims); ylim(ylims); zlim(zlims);
daspect([1 1 1]); grid on; box on; view(-35, 22);
title(sprintf('Fig 4 | W8 v3 — voxels \\geq %d dB + true box', pcThr));
rotate3d on;

%% ===== EVALUATION =====
fprintf('\n========== EVALUATION — W8 v3: Multi-aspect tomography ==========\n');

% ค่า reconstruction ที่กึ่งกลางแต่ละหน้าของกล่องจริง (ต้องเด่นทุกหน้าเว้น Bottom)
faceCtr = [boxCenter(1)-boxL/2, boxCenter(2),          boxH/2;   % Front
           boxCenter(1)+boxL/2, boxCenter(2),          boxH/2;   % Back
           boxCenter(1),        boxCenter(2)-boxW/2,   boxH/2;   % Left
           boxCenter(1),        boxCenter(2)+boxW/2,   boxH/2;   % Right
           boxCenter(1),        boxCenter(2),          boxH;     % Top
           boxCenter(1),        boxCenter(2),          0     ];  % Bottom
fprintf('\n[Recovered level at each face centre]  (0 dB = image peak)\n');
for f = 1:6
    vA1 = interp3(Xg,Yg,Zg, magAspect1, faceCtr(f,2), faceCtr(f,1), faceCtr(f,3), 'linear', 0);
    vMC = interp3(Xg,Yg,Zg, magSm,      faceCtr(f,2), faceCtr(f,1), faceCtr(f,3), 'linear', 0);
    fprintf('  %-7s : 1 aspect %7.1f dB   |  %d aspects %7.1f dB\n', ...
        faceLabels{f}, 20*log10(vA1+eps), Nasp, 20*log10(vMC+eps));
end

[e3r,e3c,e3z] = peakExtents(magComb, xScene,yScene,zScene);
fprintf('\n[-3 dB extent through peak]  range %.2f | cross %.2f | height %.2f m\n', e3r,e3c,e3z);
fprintf('\n[Note]\n');
fprintf('  Coherent ภายใน aspect (height focusing), incoherent ข้าม aspect (max).\n');
fprintf('  Bottom ไม่มีวันเห็น (กล่องวางบนพื้น) -> กล่อง 5 หน้า = ครบตามฟิสิกส์.\n');

%% ===== SAVE FIGURES =====
figDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end
exportgraphics(figure(1), fullfile(figDir,'fig_w8v3_geometry.png'),   'Resolution',150);
exportgraphics(figure(2), fullfile(figDir,'fig_w8v3_aspects.png'),    'Resolution',150);
exportgraphics(figure(3), fullfile(figDir,'fig_w8v3_iso4views.png'),  'Resolution',150);
exportgraphics(figure(4), fullfile(figDir,'fig_w8v3_pointcloud.png'), 'Resolution',150);
fprintf('\nFigures saved to %s\n', figDir);


%%%% ===== LOCAL FUNCTIONS =====

function [pos, amp, faceIdx] = buildBoxScatterers(bc, L, Wd, H, Nf, surfaceModel)
% สร้างจุดครบ 6 หน้า (ไม่ cull ที่นี่ — mask ราย aspect ข้างนอก)
%   faceIdx: 1 x N หมายเลขหน้าของแต่ละจุด (1..6)
    rng(1);
    hL=L/2; hW=Wd/2; hH=H/2; u=linspace(-1,1,Nf); pf=Nf^2;
    [V,W]=meshgrid(u*hW,u*hH);
    F=cell(1,6);
    F{1}=[(bc(1)-hL)*ones(1,pf); (bc(2)+V(:))'; (bc(3)+hH+W(:))'];   % Front
    F{2}=[(bc(1)+hL)*ones(1,pf); (bc(2)+V(:))'; (bc(3)+hH+W(:))'];   % Back
    [U,W]=meshgrid(u*hL,u*hH);
    F{3}=[(bc(1)+U(:))'; (bc(2)-hW)*ones(1,pf); (bc(3)+hH+W(:))'];   % Left
    F{4}=[(bc(1)+U(:))'; (bc(2)+hW)*ones(1,pf); (bc(3)+hH+W(:))'];   % Right
    [U,V]=meshgrid(u*hL,u*hW);
    F{5}=[(bc(1)+U(:))'; (bc(2)+V(:))'; (bc(3)+H)*ones(1,pf)];       % Top
    F{6}=[(bc(1)+U(:))'; (bc(2)+V(:))'; zeros(1,pf)];                % Bottom
    pos=[F{:}];
    faceIdx = repelem(1:6, pf);
    amp = ones(1, 6*pf);
    if strcmpi(surfaceModel,'rough')
        amp = amp .* exp(1j*2*pi*rand(1,6*pf));   % เฟสคงที่ต่อจุด ใช้ร่วมทุก aspect
    end
end


function S = simPassAnalytic(pos, amp, plat, rbins, rngres, lambda)
% analytic range-compressed data: sinc ที่ slant range + carrier phase
    Naz=size(plat,2); Nr=numel(rbins); S=zeros(Nr,Naz);
    for k=1:size(pos,2)
        R = sqrt(sum((plat-pos(:,k)).^2, 1));
        S = S + amp(k) * sinc((rbins-R)/rngres) .* exp(-1j*4*pi*R/lambda);
    end
end


function img = bpOnePass(S, plat, rbins, voxT, fc, c, sz)
% 3D BP หนึ่ง pass — manual linear interp (เร็วกว่า interp1 มาก)
%   voxT: Nvox x 3 [range x, cross y, z]
    Nr = numel(rbins); dr = rbins(2)-rbins(1); r1 = rbins(1);
    Nvox = size(voxT,1);
    acc  = zeros(Nvox,1);
    Naz  = size(plat,2);
    for a = 1:Naz
        d  = voxT - plat([1 2 3], a).';       % Nvox x 3 (plat rows: x,y,z)
        sr = sqrt(sum(d.^2, 2));              % Nvox x 1
        fi = (sr - r1)/dr + 1;
        i0 = floor(fi);
        v  = (i0 >= 1) & (i0 < Nr);
        w  = fi - i0;
        Sa = S(:,a);
        iv = zeros(Nvox,1);
        iv(v) = (1-w(v)).*Sa(i0(v)) + w(v).*Sa(i0(v)+1);
        acc = acc + iv .* exp(1j*4*pi*fc/c*sr);
    end
    img = reshape(acc, sz);
end


function [V,F] = makeBoxWire(bc, hl, hw, H)
% wireframe กล่องจริง (X=cross, Y=range, Z=height ; z=0..H) — สำหรับ patch (Fig1 geometry)
    cx=bc(2); cy=bc(1);
    V=[cx-hw,cy-hl,0; cx+hw,cy-hl,0; cx+hw,cy+hl,0; cx-hw,cy+hl,0; ...
       cx-hw,cy-hl,H; cx+hw,cy-hl,H; cx+hw,cy+hl,H; cx-hw,cy+hl,H];
    F=[1 2 3 4; 5 6 7 8; 1 2 6 5; 3 4 8 7; 1 4 8 5; 2 3 7 6];
end

function [Xe,Ye,Ze] = boxEdges(bc, hl, hw, H)
% 12 ขอบเป็น polyline เดียว (คั่น NaN) -> plot3 หมุนไม่บั๊ก (X=cross,Y=range,Z=height)
    cx=bc(2); cy=bc(1);
    C=[cx-hw,cy-hl,0; cx+hw,cy-hl,0; cx+hw,cy+hl,0; cx-hw,cy+hl,0; ...
       cx-hw,cy-hl,H; cx+hw,cy-hl,H; cx+hw,cy+hl,H; cx-hw,cy+hl,H];
    E=[1 2;2 3;3 4;4 1; 5 6;6 7;7 8;8 5; 1 5;2 6;3 7;4 8];
    Xe=nan(1,3*size(E,1)); Ye=Xe; Ze=Xe;
    for e=1:size(E,1)
        j=3*(e-1)+1;
        Xe(j:j+1)=C(E(e,:),1); Ye(j:j+1)=C(E(e,:),2); Ze(j:j+1)=C(E(e,:),3);
    end
end


function [rExt,cExt,zExt] = peakExtents(mag, xS, yS, zS)
    [~,pk]=max(mag(:)); [iy,ix,iz]=ind2sub(size(mag),pk);
    rExt=ax3(yS, squeeze(mag(:,ix,iz)));
    cExt=ax3(xS, squeeze(mag(iy,:,iz)));
    zExt=ax3(zS, squeeze(mag(iy,ix,:)));
end

function w = ax3(ax, cut)
    ax=ax(:); m=abs(cut(:)); d=20*log10(m/max(m)+eps); a=ax(d>=-3);
    if isempty(a), w=0; else, w=max(a)-min(a); end
end
