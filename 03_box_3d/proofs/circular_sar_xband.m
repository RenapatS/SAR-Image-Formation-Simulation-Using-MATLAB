%% W7_CircularSAR_Xband — CIRCULAR SAR รุ่น X-BAND (9.6 GHz) ไว้เทียบ RESOLUTION
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  โคลนจาก W7_CircularSAR.m (C-band 4 GHz) เปลี่ยนแค่ "band" เพื่อตอบคำถาม:
%     "ถ้า sim ใช้ X-band (เหมือน GOTCHA จริง) resolution ดีขึ้นแค่ไหน?"
%
%  ตัวแปรที่คุมให้ "เท่าเดิม" (fair comparison):
%     - Bandwidth 50 MHz เท่าเดิม  -> range res = c/2BW = 3 m เท่าเดิม
%       ** range res ขึ้นกับ BANDWIDTH ไม่ใช่ band -> เปลี่ยน fc เฉย ๆ ไม่ช่วย range **
%     - Elevation baseline = 17.5 m เท่าเดิม -> ใช้สูตร height res เทียบตรง ๆ ได้
%
%  ตัวแปรที่เปลี่ยน:
%     - fc 4 GHz -> 9.6 GHz  (lambda 7.5 cm -> 3.1 cm)  = X-band แบบ GOTCHA
%       => cross-range และ height res ละเอียดขึ้น ~ lambda_C/lambda_X = 7.5/3.1 = 2.4 เท่า
%     - Npass 8 -> 15, dzPass 2.5 -> 1.25 m  (baseline (15-1)*1.25 = 17.5 m คงเดิม!)
%       ** ทำไมต้องซอยชั้นถี่ขึ้น: lambda สั้นลง -> unambiguous height เล็กลง
%          height ambiguity = lambda*r0/(2*dzPass); ถ้าใช้ dzPass 2.5 เดิมที่ X-band
%          จะได้ ~7 m ซึ่ง < ความสูงกล่อง 10 m -> ยอดกล่อง "wrap" (aliasing แนวสูง)
%          ลด dzPass เป็น 1.25 -> ambiguity ~14 m > 10 m -> ปลอดภัย
%          และคง baseline 17.5 m -> height RESOLUTION ยังเทียบกับ C-band ได้ตรง **
%
%  ค่าที่คาดไว้ (ทฤษฎี, r0 = sqrt(1000^2+500^2) = 1118 m, cos^2(th)):
%     |            | lambda  | range res | cross res(arc 15deg) | height res |
%     | C-band 4G  | 7.5 cm  |   3 m     |  ~0.14 m            |  ~3.0 m    |
%     | X-band 9.6G| 3.1 cm  |   3 m     |  ~0.06 m            |  ~1.25 m   |
%     (height res ใช้สูตร lambda*r0/(2*B*cos^2 th), th=atan(500/1000))
%  -> เลขจริงวัดจาก [-3 dB extent] ที่พิมพ์ท้าย script (รันแล้วเทียบกับไฟล์เดิม)
%
%  OPTIONAL: ถ้าอยากได้ range res 0.24 m แบบ GOTCHA ด้วย ให้ตั้ง rngres = 0.24
%     แต่ต้องเพิ่ม Nfreq (เช่น ~800) และ Nzp ให้ range window ยังกว้างพอ ไม่งั้น alias
%     (BW 624 MHz + Nfreq 64 -> window แค่ ~15 m < กล่อง 30 m) -> ช้าลงมาก
%
%  RUNTIME: ~10-16 นาที (Npass 15 = ~2 เท่าของไฟล์เดิม)

clear; clc; close all;

%% ===== RADAR PARAMETERS (X-BAND) =====
c   = physconst('LightSpeed');
fc  = 9.6e9;    lambda = c/fc;         % << X-band 9.6 GHz (เดิม 4e9)
rngres = 3;                            % m  (BW เท่าเดิม -> range res เท่าเดิม)
bw  = c/(2*rngres);                    % 50 MHz

Nfreq = 64;                            % จำนวน frequency samples (แบบ GOTCHA)
fgrid = fc + linspace(-bw/2, bw/2, Nfreq);
df    = fgrid(2)-fgrid(1);
wFreq = taylorwin(Nfreq, 4, -35).';    % Taylor window กด range sidelobe
Nzp   = 512;                           % zero-pad IFFT
drAxis = ((0:Nzp-1).' - Nzp/2) * c/(2*df*Nzp);   % differential range axis (m)
kc = 4*pi*fc/c;
fprintf('[X-BAND] fc %.1f GHz | lambda %.1f cm | BW %.0f MHz | range window ~%.0f m\n', ...
    fc/1e9, lambda*100, bw/1e6, c/(2*df));

%% ===== CIRCULAR GEOMETRY (GOTCHA-style) =====
Rg     = 1000;
baseH  = 500;
Npass  = 15;                           % << ซอยชั้นถี่ขึ้น (เดิม 8) กัน height ambiguity
dzPass = 1.25;                         % << (15-1)*1.25 = 17.5 m = baseline เดิม
arcDeg = 15;
azStep = 0.05;
arcs   = 0 : arcDeg : 360-arcDeg;
Narc   = numel(arcs);
NazArc = round(arcDeg/azStep);
fprintf('Circular: %d arcs x %d deg | %d passes x %.2f m (baseline %.1f m) | %d pulses/arc/pass\n', ...
    Narc, arcDeg, Npass, dzPass, (Npass-1)*dzPass, NazArc);

%% ===== BOX SCATTERERS =====
boxCenter = [1000; 0; 0];
boxL = 30; boxW = 20; boxH = 10;
Nf = 20;                               % 400 จุด/หน้า (ห่าง ~1-1.6 m)
surfaceModel = 'rough';

[scatPosAll, scatAmpAll, scatFace] = buildBoxScatterers( ...
    boxCenter, boxL, boxW, boxH, Nf, surfaceModel);
faceNormal = [ -1 0 0;  1 0 0;  0 -1 0;  0 1 0;  0 0 1;  0 0 -1];
faceLabels = {'Front','Back','Left','Right','Top','Bottom'};
fprintf('Scatterers (all faces): %d | surface: %s\n', size(scatPosAll,2), surfaceModel);

%% ===== VOXEL GRID =====
voxStep = 1.0;
xScene = (boxCenter(2)-25) : voxStep : (boxCenter(2)+25);
yScene = (boxCenter(1)-25) : voxStep : (boxCenter(1)+25);
zScene = -5               : voxStep : 20;
[Xg, Yg, Zg] = meshgrid(xScene, yScene, zScene);
sz   = size(Xg);
voxT = [Yg(:), Xg(:), Zg(:)];
fprintf('Voxel grid: %d x %d x %d = %d voxels\n', sz(1),sz(2),sz(3),numel(Xg));

%% ===== CIRCULAR TOMOGRAPHIC BP (multilook fusion) =====
wElev  = hamming(Npass);
Icomb  = zeros(sz);                % สะสม intensity (multilook)
magArc1 = [];
tAll = tic;

for iArc = 1:Narc
    azv = arcs(iArc) + (0:NazArc-1)*azStep;
    platXY = boxCenter(1:2) + Rg*[cosd(azv); sind(azv)];

    % occlusion ราย arc
    Hm = baseH + (Npass-1)/2*dzPass;
    losA = [platXY; Hm*ones(1,NazArc)] - boxCenter;
    visF = false(1,6);
    for f = 1:6, visF(f) = any(faceNormal(f,:)*losA > 0); end
    ampA = scatAmpAll .* visF(scatFace);
    keep = ampA ~= 0;
    posK = scatPosAll(:,keep);  ampK = ampA(keep);

    bpA = zeros(sz);
    for p = 1:Npass
        Hp   = baseH + (p-1)*dzPass;
        plat = [platXY; Hp*ones(1,NazArc)];
        [rc, rref] = simPassFreq(posK, ampK, plat, boxCenter, fgrid, wFreq, Nzp, c);
        bpA  = bpA + wElev(p) * bpOnePassDiff(rc, rref, plat, drAxis, voxT, kc, sz);
    end
    magA = abs(bpA);
    if iArc == 1, magArc1 = magA; end
    Icomb = Icomb + magA.^2;                    % multilook (เฉลี่ย intensity)
    fprintf('  arc %2d/%2d (az %3d deg) | faces: %s | %.0f s\n', ...
        iArc, Narc, arcs(iArc), strjoin(faceLabels(visF),','), toc(tAll));
end
magComb = sqrt(Icomb / Narc);
fprintf('Circular BP complete (%.1f s total).\n', toc(tAll));

magComb = magComb ./ max(magComb(:));
magArc1 = magArc1 ./ max(magArc1(:));

%% ===== THEORY =====
r0  = sqrt(Rg^2 + baseH^2);
th  = atan(baseH/Rg);  cos2 = cos(th)^2;
B_elev = (Npass-1)*dzPass;
fprintf('\n[Design X-band] elev baseline %.1f m -> height res ~%.2f m | ambiguity %.1f m\n', ...
    B_elev, lambda*r0/(2*B_elev*cos2), lambda*r0/(2*dzPass*cos2));
fprintf('[Ref  C-band ] lambda 7.5cm ที่ baseline เดิม -> height res ~%.2f m (เทียบ)\n', ...
    (c/4e9)*r0/(2*B_elev*cos2));

%% ===== POST-PROCESSING FOR DISPLAY =====
magSm = smooth3(magComb, 'gaussian', [5 5 5], 1.5);
magSm = magSm ./ max(magSm(:));
magSmdB = 20*log10(magSm + eps);

[Xe,Ye,Ze] = boxEdges(boxCenter, boxL/2, boxW/2, boxH);
boxCr = boxCenter(2)+boxW/2*[-1 1 1 -1 -1];
boxRg = boxCenter(1)+boxL/2*[-1 -1 1 1 -1];
xlims = [-18 18]; ylims = [982 1018]; zlims = [-3 15];

%% ===== FIGURE 1: 1 ARC vs FULL CIRCLE — top-down MIP =====
mipA1 = squeeze(max(magArc1,[],3));
mipMC = squeeze(max(magComb,[],3));

figure(1); set(gcf,'Name','Fig1 Xband: arc vs full circle','Color','k','Position',[80 60 950 460]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

axa = nexttile;
imagesc(xScene, yScene, 20*log10(mipA1+eps)); set(axa,'YDir','normal'); clim([-30 0]); colormap(axa,'gray');
hold on; plot(boxCr,boxRg,'c--','LineWidth',1.5); hold off;
xlabel('Cross (m)','Color','w'); ylabel('Range (m)','Color','w');
title(sprintf('single %d\\circ arc',arcDeg),'Color','w');
set(axa,'Color','k','XColor','w','YColor','w');

axb = nexttile;
imagesc(xScene, yScene, 20*log10(mipMC+eps)); set(axb,'YDir','normal'); clim([-30 0]); colormap(axb,'gray');
hold on; plot(boxCr,boxRg,'c--','LineWidth',1.5); hold off;
xlabel('Cross (m)','Color','w'); ylabel('Range (m)','Color','w');
title('full 360\circ (multilook)','Color','w');
set(axb,'Color','k','XColor','w','YColor','w');
cb = colorbar; cb.Color='w'; ylabel(cb,'dB','Color','w'); cb.Layout.Tile='east';
sgtitle('Fig 1 | X-band 9.6 GHz — Taylor + multilook','Color','w');

%% ===== FIGURE 2: isosurface 4 มุมมอง =====
figure(2); set(gcf,'Name','Fig2 Xband: full-box isosurface','Color','w','Position',[60 60 1150 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

isoLv  = [0.45, 0.2];
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
    plot3(ax, Xe, Ye, Ze, '-', 'Color',[0.1 0.3 0.9], 'LineWidth',1.5);
    hold(ax,'off');
    xlabel('Cross (m)'); ylabel('Range (m)'); zlabel('Height (m)');
    title(viewNm{k});
    xlim(xlims); ylim(ylims); zlim(zlims);
    daspect(ax,[1 1 1]); grid(ax,'on'); box(ax,'on');
    view(ax, viewAz(k), viewEl(k));
    camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
end
sgtitle(sprintf('Fig 2 | X-band circular SAR box (isosurface @ %.0f/%.0f dB)', ...
    20*log10(isoLv(1)), 20*log10(isoLv(2))));
rotate3d on;

%% ===== FIGURE 3: point cloud =====
pcThr = -10;    % << เดิม -16 ทำให้ point cloud โดนพื้น speckle (~-15 dB) ฟุ้ง;
                %    -10 ตัดพื้นออก เหลือเฉพาะกล่อง (ผนัง -4..-7 dB) ให้คมเท่า isosurface
idxPC = find(magSmdB >= pcThr);
figure(3); set(gcf,'Name','Fig3 Xband: voxel point cloud','Color','w','Position',[120 80 780 640]);
scatter3(Xg(idxPC), Yg(idxPC), Zg(idxPC), 12, magSmdB(idxPC), 'filled', ...
    'MarkerFaceAlpha', 0.4);
hold on;
plot3(Xe, Ye, Ze, '-', 'Color',[0.1 0.3 0.9], 'LineWidth',1.8);
hold off;
colormap(jet); clim([pcThr 0]);
cb3 = colorbar; ylabel(cb3,'dB');
xlabel('Cross (m)'); ylabel('Range (m)'); zlabel('Height (m)');
xlim(xlims); ylim(ylims); zlim(zlims);
daspect([1 1 1]); grid on; box on; view(-35, 22);
title(sprintf('Fig 3 | X-band — voxels \\geq %d dB + true box', pcThr));
rotate3d on;

%% ===== EVALUATION =====
fprintf('\n========== EVALUATION — X-band Circular SAR ==========\n');
faceCtr = [boxCenter(1)-boxL/2, boxCenter(2),        boxH/2;
           boxCenter(1)+boxL/2, boxCenter(2),        boxH/2;
           boxCenter(1),        boxCenter(2)-boxW/2, boxH/2;
           boxCenter(1),        boxCenter(2)+boxW/2, boxH/2;
           boxCenter(1),        boxCenter(2),        boxH;
           boxCenter(1),        boxCenter(2),        0     ];
fprintf('\n[Recovered level at each face centre]  (0 dB = image peak)\n');
for f = 1:6
    vC = interp3(Xg,Yg,Zg, magSm, faceCtr(f,2), faceCtr(f,1), faceCtr(f,3), 'linear', 0);
    fprintf('  %-7s : %7.1f dB\n', faceLabels{f}, 20*log10(vC+eps));
end

vWall = interp3(Xg,Yg,Zg, magSm, 0, boxCenter(1)-boxL/2, boxH/2, 'linear', 0);
vIn   = interp3(Xg,Yg,Zg, magSm, 0, boxCenter(1),        boxH/2, 'linear', 0);
vOut  = interp3(Xg,Yg,Zg, magSm, 20, boxCenter(1),       boxH/2, 'linear', 0);
fprintf('\n[Contrast]  wall %.1f dB | inside box %.1f dB | outside %.1f dB\n', ...
    20*log10(vWall+eps), 20*log10(vIn+eps), 20*log10(vOut+eps));
fprintf('  -> wall-to-interior contrast = %.1f dB\n', 20*log10(vWall/(vIn+eps)+eps));

[e3r,e3c,e3z] = peakExtents(magComb, xScene,yScene,zScene);
fprintf('\n[-3 dB extent through peak — X-band]  range %.2f | cross %.2f | height %.2f m\n', e3r,e3c,e3z);
fprintf('   (เทียบไฟล์เดิม C-band ~ range 3 | cross ~0.14 | height ~3.0 m)\n');
fprintf('   -> range แทบไม่เปลี่ยน (BW เท่ากัน) | cross+height คมขึ้น ~2.4x (lambda สั้นลง)\n');

%% ===== SAVE FIGURES =====
figDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end
exportgraphics(figure(1), fullfile(figDir,'fig_xband_arc_vs_full.png'),'Resolution',150);
exportgraphics(figure(2), fullfile(figDir,'fig_xband_iso4views.png'),  'Resolution',150);
exportgraphics(figure(3), fullfile(figDir,'fig_xband_pointcloud.png'), 'Resolution',150);
fprintf('\nFigures saved to %s\n', figDir);


%%%% ===== LOCAL FUNCTIONS =====

function [pos, amp, faceIdx] = buildBoxScatterers(bc, L, Wd, H, Nf, surfaceModel)
% จุดครบ 6 หน้า (mask ราย arc ข้างนอก)
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
        amp = amp .* exp(1j*2*pi*rand(1,6*pf));
    end
end


function [rc, rref] = simPassFreq(pos, amp, plat, sceneCtr, fgrid, wFreq, Nzp, c)
% frequency-domain forward model (โครงเดียวกับ GOTCHA phase history)
    Naz  = size(plat,2);
    rref = sqrt(sum((plat - sceneCtr).^2, 1));               % 1 x Naz
    dR   = sqrt( (pos(1,:)'-plat(1,:)).^2 + (pos(2,:)'-plat(2,:)).^2 ...
               + (pos(3,:)'-plat(3,:)).^2 ) - rref;          % Nscat x Naz
    Nf = numel(fgrid);
    S  = zeros(Nf, Naz);
    for n = 1:Nf
        S(n,:) = wFreq(n) * (amp * exp(-1j*4*pi*fgrid(n)/c * dR));
    end
    rc = fftshift(ifft(S, Nzp, 1), 1);
end


function img = bpOnePassDiff(rc, rref, plat, drAxis, voxT, kc, sz)
% 3D BP บนแกน differential range
    Nzp = numel(drAxis); ddr = drAxis(2)-drAxis(1); d1 = drAxis(1);
    Nvox = size(voxT,1);
    acc  = zeros(Nvox,1);
    for a = 1:size(plat,2)
        d  = voxT - plat(:,a).';
        dv = sqrt(sum(d.^2, 2)) - rref(a);
        fi = (dv - d1)/ddr + 1;
        i0 = floor(fi);
        v  = (i0 >= 1) & (i0 < Nzp);
        w  = fi - i0;
        Sa = rc(:,a);
        iv = zeros(Nvox,1);
        iv(v) = (1-w(v)).*Sa(i0(v)) + w(v).*Sa(i0(v)+1);
        acc = acc + iv .* exp(1j*kc*dv);
    end
    img = reshape(acc, sz);
end


function [Xe,Ye,Ze] = boxEdges(bc, hl, hw, H)
% 12 ขอบของกล่องเป็น polyline เดียว (คั่นด้วย NaN)
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
