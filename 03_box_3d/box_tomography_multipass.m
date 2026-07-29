%% MULTI-PASS TOMOGRAPHIC SAR — บินซ้ำหลายระดับความสูงเพื่อกู้มิติความสูง
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  ต่อจาก box_bp_3d_singlepass.m : ส่วน simulation/BP เหมือนเดิมทุกประการ
%  แก้เฉพาะการแสดงผล 3D (แนวเดียวกับ W7 v2):
%
%  (1) เดิม Fig 3 ระบายสีบน "จุดกล่องจริง" -> ไม่ได้โชว์ volume ที่
%      tomography สร้างได้จริง  แก้: เพิ่ม isosurface ของ reconstruction
%  (2) volume ดิบมี speckle ('rough' surface) + carrier fringe ->
%      isosurface ตรง ๆ จะขรุขระ/เป็นลาย  แก้: smooth3 ก่อนวาด
%      (evaluation ยังใช้ข้อมูลดิบ)
%  (3) มุมมองเดียว + ไม่ daspect -> ดูยาก สัดส่วนเพี้ยน
%      แก้: 4 มุมมอง, daspect([1 1 1]), isosurface 2 ชั้นโปร่งแสง,
%      crop แกนรอบกล่อง, rotate3d on
%  (4) เพิ่ม point cloud ของ voxel ที่แรงกว่า threshold (อ่านง่ายอีกแบบ)
%
%  หมายเหตุ: Back/Bottom ถูก occlusion ตัดออก (เรดาร์ไม่เคยเห็น) ->
%  ภาพที่ถูกต้องทางฟิสิกส์จะโชว์ได้มากสุด 4 หน้า ไม่ใช่กล่องปิด 6 ด้าน

clear; clc; close all;

%% ===== RADAR / GEOMETRY PARAMETERS =====
c   = physconst('LightSpeed');
fc  = 4e9;      lambda = c/fc;
bw  = c/(2*3);                 % 3 m range resolution
fs  = 120e6;
prf = 1000;
speed = 100;   flightDuration = 4;
rngres = c/(2*bw);

baseH  = 500;
groundR0 = 1000;

%% ===== MULTI-PASS CONFIG =====
Npass  = 21;
dzPass = 1.0;
azDecim = 5;
elevWindow = true;

Naz_full = flightDuration*prf + 1;
azIdx    = 1:azDecim:Naz_full;
Naz      = numel(azIdx);
t        = (azIdx-1)/prf;
ry       = -200 + speed*t;

%% ===== BOX SCATTERERS + PHYSICS RCS + OCCLUSION =====
boxCenter = [groundR0; 0; 0];
boxL = 30; boxW = 20; boxH = 10;
Nf = 10;
ampMode = 'uniform';
surfaceModel = 'rough';

platMid = [zeros(1,Naz); ry; (baseH+(Npass-1)/2*dzPass)*ones(1,Naz)];
faceNormal = [ -1 0 0;  1 0 0;  0 -1 0;  0 1 0;  0 0 1;  0 0 -1];
faceLabels = {'Front','Back','Left','Right','Top','Bottom'};
losAll = platMid - boxCenter;
faceVisible = true(1,6);
for f = 1:6, faceVisible(f) = any(faceNormal(f,:)*losAll > 0); end

[scatPos, scatAmp] = buildBoxScatterers(boxCenter, boxL, boxW, boxH, Nf, lambda, faceVisible, ampMode, surfaceModel);
fprintf('Surface model: %s\n', surfaceModel);
fprintf('Occlusion | visible: %s | hidden: %s\n', ...
    strjoin(faceLabels(faceVisible),', '), strjoin(faceLabels(~faceVisible),', '));
fprintf('Scatterers (visible): %d\n', size(scatPos,2));

%% ===== VOXEL GRID =====
voxStep = 1.0;
xScene = (boxCenter(2)-30) : voxStep : (boxCenter(2)+30);   % cross
yScene = (boxCenter(1)-20) : voxStep : (boxCenter(1)+20);   % range
zScene = -5               : voxStep : 25;                   % height
[Xg, Yg, Zg] = meshgrid(xScene, yScene, zScene);
Ny = numel(yScene); Nx = numel(xScene); Nz = numel(zScene);
fprintf('Voxel grid: %d x %d x %d = %d voxels (step %.1f m)\n', Ny,Nx,Nz,Ny*Nx*Nz,voxStep);

rbins = (1080 : c/2/fs : 1180).';
Nr = numel(rbins);

%% ===== TOMOGRAPHIC 3D BACKPROJECTION (multi-pass) =====
fprintf('Running multi-pass 3D backprojection (%d passes x %d pulses)...\n', Npass, Naz);
tic;
sz  = size(Xg);
vox = [Yg(:)'; Xg(:)'; Zg(:)'];
bp3 = zeros(sz);
bp1 = zeros(sz);
if elevWindow, wElev = hamming(Npass); else, wElev = ones(Npass,1); end
midPass = round((Npass+1)/2);
for p = 1:Npass
    Hp   = baseH + (p-1)*dzPass;
    plat = [zeros(1,Naz); ry; Hp*ones(1,Naz)];
    S = simPassAnalytic(scatPos, scatAmp, plat, rbins, rngres, lambda);
    imgP = zeros(sz);
    for a = 1:Naz
        rp = plat(:,a);
        sr = sqrt(sum((vox-rp).^2, 1));
        iv = interp1(rbins, S(:,a), sr, 'linear', 0);
        imgP = imgP + reshape(iv .* exp(1j*4*pi*fc/c .* sr), sz);
    end
    bp3 = bp3 + wElev(p) * imgP;
    if p == midPass, bp1 = imgP; end
    fprintf('  pass %2d/%2d  (H = %.0f m)\n', p, Npass, Hp);
end
tBP = toc;
fprintf('Tomographic BP complete (%.1f s).\n', tBP);

mag3 = abs(bp3); mag3 = mag3 ./ max(mag3(:));
mag1 = abs(bp1); mag1 = mag1 ./ max(mag1(:));

%% ===== POST-PROCESSING FOR DISPLAY  << ใหม่ (v2) =====
%  ลบ speckle/fringe ด้วย incoherent smoothing เพื่อการแสดงผลเท่านั้น
%  kernel เบากว่า W7 เพราะ height โฟกัสแล้ว (ไม่อยากละลายผนังกล่อง)
magSm = smooth3(mag3, 'gaussian', [3 3 5], 1.2);
magSm = magSm ./ max(magSm(:));
magSmdB = 20*log10(magSm + eps);

[Vb,Fb] = makeBoxWire(boxCenter, boxL/2, boxW/2, boxH);
[Xe,Ye,Ze] = boxEdges(boxCenter, boxL/2, boxW/2, boxH);   % กรอบเส้น (หมุนไม่บั๊ก)
boxCr = boxCenter(2)+boxW/2*[-1 1 1 -1 -1];
boxRg = boxCenter(1)+boxL/2*[-1 -1 1 1 -1];

% แกน crop รอบกล่อง (ตัดพื้นที่ว่างออก -> กล่องเต็มเฟรม อ่านง่าย)
xlims = [-18 18]; ylims = [982 1018]; zlims = [-3 15];

%% ===== THEORY: elevation resolution & ambiguity =====
r0 = sqrt(groundR0^2 + baseH^2);
th = atan(baseH/groundR0);
cos2 = cos(th)^2;
B_elev  = (Npass-1)*dzPass;
dz_res  = lambda*r0 / (2*B_elev*cos2);
Ha_amb  = lambda*r0 / (2*dzPass*cos2);
fprintf('\n[Tomography design]\n');
fprintf('  Passes: %d | dzPass: %.1f m | elevation baseline B = %.0f m\n', Npass, dzPass, B_elev);
fprintf('  Height resolution (theory) : %.2f m\n', dz_res);
fprintf('  Ambiguous height span      : %.1f m  (grid height ~%.0f m)\n', Ha_amb, zScene(end)-zScene(1));

%% ===== FIGURE 1: SINGLE vs MULTI — range-height MIP (เหมือนเดิม) =====
mip1 = squeeze(max(mag1,[],2));
mip3 = squeeze(max(mag3,[],2));

figure(1); set(gcf,'Name','Fig1: single vs multi (range-height)','Color','k','Position',[60 80 950 460]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

axa = nexttile;
imagesc(zScene, yScene, 20*log10(mip1+eps)); set(axa,'YDir','normal'); clim([-40 0]); colormap(axa,'gray');
hold on; yline(boxCenter(1)-boxL/2,'c--'); yline(boxCenter(1)+boxL/2,'c--');
xline(0,'c:'); xline(boxH,'c:'); hold off;
xlabel('Height (m)','Color','w'); ylabel('Range (m)','Color','w');
title('SINGLE pass — height smeared (layover)','Color','w');
set(axa,'Color','k','XColor','w','YColor','w');

axb = nexttile;
imagesc(zScene, yScene, 20*log10(mip3+eps)); set(axb,'YDir','normal'); clim([-40 0]); colormap(axb,'gray');
hold on; yline(boxCenter(1)-boxL/2,'c--'); yline(boxCenter(1)+boxL/2,'c--');
xline(0,'c:'); xline(boxH,'c:'); hold off;
xlabel('Height (m)','Color','w'); ylabel('Range (m)','Color','w');
title(sprintf('MULTI pass (%d) — height focused',Npass),'Color','w');
set(axb,'Color','k','XColor','w','YColor','w');
cb = colorbar; cb.Color='w'; ylabel(cb,'dB','Color','w'); cb.Layout.Tile='east';
sgtitle('Fig 1 | W8 v2 — elevation aperture collapses the layover (cyan dots = true box 0..10 m)','Color','w');

%% ===== FIGURE 2: multi-pass MIP — 3 orthogonal planes (เหมือนเดิม) =====
mipXY = squeeze(max(mag3,[],3));
mipRZ = squeeze(max(mag3,[],2));
mipXZ = squeeze(max(mag3,[],1));

figure(2); set(gcf,'Name','Fig2: multi-pass MIP planes','Color','k','Position',[80 60 1250 420]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

a1=nexttile; imagesc(xScene,yScene,20*log10(mipXY+eps)); set(a1,'YDir','normal'); clim([-40 0]); colormap(a1,'gray');
hold on; plot(boxCr,boxRg,'c--','LineWidth',1.5); hold off;
xlabel('Cross-Range (m)','Color','w'); ylabel('Range (m)','Color','w'); title('MIP top-down','Color','w');
set(a1,'Color','k','XColor','w','YColor','w');

a2=nexttile; imagesc(zScene,yScene,20*log10(mipRZ+eps)); set(a2,'YDir','normal'); clim([-40 0]); colormap(a2,'gray');
hold on; xline(0,'c:'); xline(boxH,'c:'); hold off;
xlabel('Height (m)','Color','w'); ylabel('Range (m)','Color','w'); title('MIP range-height','Color','w');
set(a2,'Color','k','XColor','w','YColor','w');

a3=nexttile; imagesc(zScene,xScene,20*log10(mipXZ+eps)); set(a3,'YDir','normal'); clim([-40 0]); colormap(a3,'gray');
hold on; xline(0,'c:'); xline(boxH,'c:'); hold off;
xlabel('Height (m)','Color','w'); ylabel('Cross-Range (m)','Color','w'); title('MIP cross-height','Color','w');
set(a3,'Color','k','XColor','w','YColor','w');
cb2=colorbar; cb2.Color='w'; ylabel(cb2,'dB','Color','w'); cb2.Layout.Tile='east';
sgtitle('Fig 2 | W8 v2 — multi-pass 3D BP: MIP onto 3 planes','Color','w');

%% ===== FIGURE 3: isosurface ของ reconstruction — 4 มุมมอง  << ใหม่ (v2) =====
%  นี่คือ "กล่อง 3D จริง" ที่ tomography กู้กลับมา (ไม่ใช่จุดกล่องจริงระบายสี)
figure(3); set(gcf,'Name','Fig3: tomographic isosurface 4 views (v2)','Color','w','Position',[60 60 1150 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

isoLv  = [0.4, 0.15];                      % ~-8 dB (ทึบ) และ ~-16 dB (โปร่ง)
isoCol = [0.85 0.25 0.15; 1.00 0.65 0.20];
isoAl  = [0.85, 0.22];

viewAz = [-35, 145,  90, -35];
viewEl = [ 22,  22,   0,  85];
viewNm = {'Perspective (หน้ากล่อง)','Perspective (หลังกล่อง)', ...
          'Front view (มองตามแนว range)','Top-down'};

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
    daspect(ax,[1 1 1]);
    grid(ax,'on'); box(ax,'on');
    view(ax, viewAz(k), viewEl(k));
    camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
end
sgtitle(sprintf(['Fig 3 | W8 v2 — tomographic reconstruction, isosurface @ %.0f/%.0f dB + true box (blue)\n' ...
    '(Back/Bottom ไม่มีในภาพ = ถูกต้อง: เรดาร์ไม่เคยเห็น)'], ...
    20*log10(isoLv(1)), 20*log10(isoLv(2))));
rotate3d on;

%% ===== FIGURE 4: point cloud ของ voxel ที่แรง  << ใหม่ (v2) =====
pcThr = -18;                               % dB
idxPC = find(magSmdB >= pcThr);
figure(4); set(gcf,'Name','Fig4: voxel point cloud (v2)','Color','w','Position',[120 80 780 640]);
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
title(sprintf('Fig 4 | W8 v2 — voxels \\geq %d dB + true box', pcThr));
rotate3d on;

%% ===== FIGURE 5: ACQUISITION GEOMETRY (เหมือนเดิม) =====
figure(5); set(gcf,'Name','Fig5: acquisition geometry','Color','w','Position',[140 60 1180 520]);
cmap = parula(Npass);

subplot(1,2,1);
patch('Vertices',Vb,'Faces',Fb,'FaceColor',[0.6 0.8 1],'FaceAlpha',0.25, ...
      'EdgeColor',[0.1 0.3 0.9],'LineWidth',1.5); hold on;
for p = 1:Npass
    Hp = baseH + (p-1)*dzPass;
    plot3(ry, zeros(1,Naz), Hp*ones(1,Naz), '-', 'Color', cmap(p,:), 'LineWidth', 1);
end
plot3(ry(1),0,baseH,'go','MarkerFaceColor','g','MarkerSize',8);
plot3(ry(end),0,baseH,'rs','MarkerFaceColor','r','MarkerSize',8);
pm = round(Npass/2); Hpm = baseH+(pm-1)*dzPass;
for cc = 1:size(Vb,1)
    plot3([0 Vb(cc,1)], [0 Vb(cc,2)], [Hpm Vb(cc,3)], '-', 'Color',[1 0.6 0.1]);
end
plot3(0,0,Hpm,'k^','MarkerFaceColor','y','MarkerSize',9);
text(0,0,Hpm+30,'  radar (mid pass)','FontSize',8);
text(0,boxCenter(1),15,'  box','FontSize',8,'Color',[0.1 0.3 0.9]);
hold off; grid on; view(35,18);
xlabel('Cross (m)'); ylabel('Range (m)'); zlabel('Height (m)');
title('A | full scene: flight tracks (colored) + beams to box');
xlim([-220 220]); ylim([-50 1050]); zlim([0 560]);

subplot(1,2,2);
for p = 1:Npass
    Hp = baseH + (p-1)*dzPass;
    plot3(ry, zeros(1,Naz), Hp*ones(1,Naz), '-', 'Color', cmap(p,:), 'LineWidth', 1.5); hold on;
end
hold off; grid on; view(20,12);
xlabel('Cross (m)'); ylabel('Range (m)'); zlabel('Height (m)');
title(sprintf('B | zoom: %d passes stacked, %.0f m apart (elevation baseline = %.0f m)', ...
    Npass, dzPass, (Npass-1)*dzPass));
zlim([baseH-2 baseH+(Npass-1)*dzPass+2]); ylim([-30 30]);
cb5 = colorbar; ylabel(cb5,'pass height (m)'); colormap(gca,parula);
clim([baseH baseH+(Npass-1)*dzPass]);
sgtitle('Fig 5 | W8 v2 — how the radar acquires: repeat the cross sweep at stepped heights','FontWeight','bold');

%% ===== EVALUATION (ข้อมูลดิบ ไม่ใช่ตัว smooth) =====
fprintf('\n========== EVALUATION — W8 v2: Multi-pass Tomographic SAR ==========\n');
[e1r,e1c,e1z] = peakExtents(mag1, xScene,yScene,zScene);
[e3r,e3c,e3z] = peakExtents(mag3, xScene,yScene,zScene);
fprintf('\n[-3 dB extent through peak]   (box: L=%.0f  W=%.0f  H=%.0f m)\n', boxL,boxW,boxH);
fprintf('  %-14s %-10s %-10s %-10s\n','config','range','cross','height');
fprintf('  %-14s %-10.2f %-10.2f %-10.2f\n','single pass', e1r,e1c,e1z);
fprintf('  %-14s %-10.2f %-10.2f %-10.2f\n','multi pass',  e3r,e3c,e3z);
fprintf('\n  -> height collapses from ~%.0f m (smeared) to ~%.1f m (focused)\n', e1z, e3z);
fprintf('  Height resolution theory: %.2f m ; ambiguity span: %.1f m\n', dz_res, Ha_amb);
fprintf('\n[Note]\n');
fprintf('  v2 = display fix (smooth3 / daspect / 4 views / crop / point cloud)\n');
fprintf('  ฟิสิกส์เท่าเดิม: เห็นได้มากสุด 4 หน้า (Back/Bottom ถูกบังจริง)\n');

%% ===== SAVE FIGURES =====
figDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end
exportgraphics(figure(1), fullfile(figDir,'fig_w8v2_single_vs_multi.png'),'Resolution',150);
exportgraphics(figure(2), fullfile(figDir,'fig_w8v2_mip.png'),           'Resolution',150);
exportgraphics(figure(3), fullfile(figDir,'fig_w8v2_iso4views.png'),     'Resolution',150);
exportgraphics(figure(4), fullfile(figDir,'fig_w8v2_pointcloud.png'),    'Resolution',150);
exportgraphics(figure(5), fullfile(figDir,'fig_w8v2_geometry.png'),      'Resolution',150);
fprintf('\nFigures saved to %s\n', figDir);


%%%% ===== LOCAL FUNCTIONS =====

function [pos, amp] = buildBoxScatterers(bc, L, Wd, H, Nf, lambda, faceVisible, ampMode, surfaceModel)
% 6 หน้าของกล่อง, ตัดหน้าที่ถูกบัง (occlusion) ออก
    if nargin<9, surfaceModel='smooth'; end
    rng(1);
    hL=L/2; hW=Wd/2; hH=H/2; u=linspace(-1,1,Nf); pf=Nf^2;
    Afb=Wd*H; Alr=L*H; Atb=L*Wd;
    sPlate=4*pi/lambda^2*[Afb^2 Afb^2 Alr^2 Alr^2 Atb^2 Atb^2];
    faceRCS=sPlate/pf^2;
    pos=[]; amp=[];
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
    for f=1:6
        if ~faceVisible(f), continue; end
        pos=[pos, F{f}];
        if strcmpi(ampMode,'physics')
            a0 = sqrt(faceRCS(f))*ones(1,pf);
        else
            a0 = ones(1,pf);
        end
        if strcmpi(surfaceModel,'rough')
            a0 = a0 .* exp(1j*2*pi*rand(1,pf));
        end
        amp=[amp, a0];
    end
end


function S = simPassAnalytic(pos, amp, plat, rbins, rngres, lambda)
% analytic range-compressed data: sinc ที่ slant range + azimuth phase
    Naz=size(plat,2); Nr=numel(rbins); S=zeros(Nr,Naz);
    for k=1:size(pos,2)
        R = sqrt(sum((plat-pos(:,k)).^2, 1));
        S = S + amp(k) * sinc((rbins-R)/rngres) .* exp(-1j*4*pi*R/lambda);
    end
end


function [V,F] = makeBoxWire(bc, hl, hw, H)
% wireframe กล่องจริง (X=cross, Y=range, Z=height ; z=0..H) — สำหรับ patch (Fig5 geometry)
    cx=bc(2); cy=bc(1);
    V=[cx-hw,cy-hl,0; cx+hw,cy-hl,0; cx+hw,cy+hl,0; cx-hw,cy+hl,0; ...
       cx-hw,cy-hl,H; cx+hw,cy-hl,H; cx+hw,cy+hl,H; cx-hw,cy+hl,H];
    F=[1 2 3 4; 5 6 7 8; 1 2 6 5; 3 4 8 7; 1 4 8 5; 2 3 7 6];
end

function [Xe,Ye,Ze] = boxEdges(bc, hl, hw, H)
% 12 ขอบเป็น polyline เดียว (คั่น NaN) -> plot3 หมุนไม่บั๊ก
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
% -3 dB width ตามแต่ละแกน ผ่าน peak voxel
    [~,pk]=max(mag(:)); [iy,ix,iz]=ind2sub(size(mag),pk);
    rExt=ax3(yS, squeeze(mag(:,ix,iz)));
    cExt=ax3(xS, squeeze(mag(iy,:,iz)));
    zExt=ax3(zS, squeeze(mag(iy,ix,:)));
end

function w = ax3(ax, cut)
    ax=ax(:); m=abs(cut(:)); d=20*log10(m/max(m)+eps); a=ax(d>=-3);
    if isempty(a), w=0; else, w=max(a)-min(a); end
end
