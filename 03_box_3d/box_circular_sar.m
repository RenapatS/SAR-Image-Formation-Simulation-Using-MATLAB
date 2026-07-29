%% CIRCULAR SAR — บินวนรอบกล่องแบบ GOTCHA เพื่อลดความฟุ้ง (sidelobe + speckle floor)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  ต่อจาก box_multiaspect_4dir.m ซึ่งได้กล่องแล้วแต่ฟุ้ง/พื้นหลังเทา สาเหตุ+วิธีแก้:
%
%  (1) RANGE SIDELOBE: เดิม forward model เป็น sinc ตรง ๆ (sidelobe -13 dB)
%      -> ฟุ้งตามแนว range รอบผนังแรง ๆ
%      แก้: เปลี่ยนเป็น frequency-domain forward model (แบบเดียวกับ GOTCHA
%      จริง: fp(f,pulse) = Σ a·exp(-j4πfR/c)) แล้วใส่ TAYLOR WINDOW ก่อน
%      IFFT -> sidelobe ลงไป ~-35 dB   [ข้อดี: pipeline นี้ = ของจริง 100%]
%
%  (2) NOISE FLOOR จากการรวม arc ด้วย MAX: max ของ speckle 24 ชุดดันพื้นขึ้น
%      แก้: รวมแบบ MULTILOOK (เฉลี่ย intensity: sqrt(mean(|I|^2))) —
%      ผนังอยู่ครบ (เห็นในหลาย arc) แต่ speckle เฉลี่ยลง -> คอนทราสต์ดีขึ้น
%
%  (3) ผนังแตกเป็นปุ่ม: จุด scatterer ห่าง ~4 m แต่ sub-aperture resolution
%      ละเอียดกว่ามาก -> BP แยกเห็นเป็นรายจุด
%      แก้: เพิ่มความหนาแน่นจุดเป็น 20x20/หน้า (ห่าง 1-1.6 m) + smooth แรงขึ้น
%
%  RUNTIME: ~5-8 นาที (จุดเยอะขึ้น แต่ sim เขียนใหม่แบบ vectorized ต่อ freq)

clear; clc; close all;

%% ===== RADAR PARAMETERS =====
c   = physconst('LightSpeed');
fc  = 4e9;      lambda = c/fc;
rngres = 3;                        % m
bw  = c/(2*rngres);                % 50 MHz

Nfreq = 64;                        % จำนวน frequency samples (แบบ GOTCHA)
fgrid = fc + linspace(-bw/2, bw/2, Nfreq);
df    = fgrid(2)-fgrid(1);
wFreq = taylorwin(Nfreq, 4, -35).';     % << (1) Taylor window กด range sidelobe
Nzp   = 512;                            % zero-pad IFFT
drAxis = ((0:Nzp-1).' - Nzp/2) * c/(2*df*Nzp);   % differential range axis (m)
kc = 4*pi*fc/c;
fprintf('Freq-domain model: %d freqs | BW %.0f MHz | range window ~%.0f m | Taylor -35 dB\n', ...
    Nfreq, bw/1e6, c/(2*df));

%% ===== CIRCULAR GEOMETRY (GOTCHA-style) =====
Rg     = 1000;
baseH  = 500;
Npass  = 8;
dzPass = 2.5;
arcDeg = 15;
azStep = 0.05;
arcs   = 0 : arcDeg : 360-arcDeg;
Narc   = numel(arcs);
NazArc = round(arcDeg/azStep);
fprintf('Circular: %d arcs x %d deg | %d passes | %d pulses/arc/pass\n', ...
    Narc, arcDeg, Npass, NazArc);

%% ===== BOX SCATTERERS =====
boxCenter = [1000; 0; 0];
boxL = 30; boxW = 20; boxH = 10;
Nf = 20;                           % << (3) 400 จุด/หน้า (ห่าง ~1-1.6 m)
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
Icomb  = zeros(sz);                % << (2) สะสม intensity (multilook)
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
fprintf('\n[Design] elev baseline %.1f m -> height res ~%.1f m | ambiguity %.0f m\n', ...
    B_elev, lambda*r0/(2*B_elev*cos2), lambda*r0/(2*dzPass*cos2));

%% ===== POST-PROCESSING FOR DISPLAY =====
magSm = smooth3(magComb, 'gaussian', [5 5 5], 1.5);   % แรงขึ้นเล็กน้อย (เชื่อมจุดผนัง)
magSm = magSm ./ max(magSm(:));
magSmdB = 20*log10(magSm + eps);

[Xe,Ye,Ze] = boxEdges(boxCenter, boxL/2, boxW/2, boxH);   % กรอบกล่องเป็นเส้น (หมุนไม่บั๊ก)
boxCr = boxCenter(2)+boxW/2*[-1 1 1 -1 -1];
boxRg = boxCenter(1)+boxL/2*[-1 -1 1 1 -1];
xlims = [-18 18]; ylims = [982 1018]; zlims = [-3 15];

%% ===== FIGURE 1: 1 ARC vs FULL CIRCLE — top-down MIP =====
mipA1 = squeeze(max(magArc1,[],3));
mipMC = squeeze(max(magComb,[],3));

figure(1); set(gcf,'Name','Fig1: arc vs full circle (v2)','Color','k','Position',[80 60 950 460]);
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
sgtitle('Fig 1 | W9 v2 — Taylor window + multilook: cleaner walls, lower floor','Color','w');

%% ===== FIGURE 2: isosurface 4 มุมมอง =====
figure(2); set(gcf,'Name','Fig2: full-box isosurface (v2)','Color','w','Position',[60 60 1150 850]);
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
    plot3(ax, Xe, Ye, Ze, '-', 'Color',[0.1 0.3 0.9], 'LineWidth',1.5);   % กรอบกล่อง (เส้น)
    hold(ax,'off');
    xlabel('Cross (m)'); ylabel('Range (m)'); zlabel('Height (m)');
    title(viewNm{k});
    xlim(xlims); ylim(ylims); zlim(zlims);
    daspect(ax,[1 1 1]); grid(ax,'on'); box(ax,'on');
    view(ax, viewAz(k), viewEl(k));
    camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
end
sgtitle(sprintf('Fig 2 | W9 v2 — circular SAR box (isosurface @ %.0f/%.0f dB)', ...
    20*log10(isoLv(1)), 20*log10(isoLv(2))));
rotate3d on;

%% ===== FIGURE 3: point cloud =====
pcThr = -10;
idxPC = find(magSmdB >= pcThr);
figure(3); set(gcf,'Name','Fig3: voxel point cloud (v2)','Color','w','Position',[120 80 780 640]);
scatter3(Xg(idxPC), Yg(idxPC), Zg(idxPC), 12, magSmdB(idxPC), 'filled', ...
    'MarkerFaceAlpha', 0.4);
hold on;
plot3(Xe, Ye, Ze, '-', 'Color',[0.1 0.3 0.9], 'LineWidth',1.8);   % กรอบกล่อง (เส้น)
hold off;
colormap(jet); clim([pcThr 0]);
cb3 = colorbar; ylabel(cb3,'dB');
xlabel('Cross (m)'); ylabel('Range (m)'); zlabel('Height (m)');
xlim(xlims); ylim(ylims); zlim(zlims);
daspect([1 1 1]); grid on; box on; view(-35, 22);
title(sprintf('Fig 3 | W9 v2 — voxels \\geq %d dB + true box', pcThr));
rotate3d on;

%% ===== EVALUATION =====
fprintf('\n========== EVALUATION — W9 v2: Circular SAR (windowed + multilook) ==========\n');
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

% wall-to-clutter contrast: ผนังหน้า vs กลางกล่อง (ข้างในว่าง) vs นอกกล่อง
vWall = interp3(Xg,Yg,Zg, magSm, 0, boxCenter(1)-boxL/2, boxH/2, 'linear', 0);
vIn   = interp3(Xg,Yg,Zg, magSm, 0, boxCenter(1),        boxH/2, 'linear', 0);
vOut  = interp3(Xg,Yg,Zg, magSm, 20, boxCenter(1),       boxH/2, 'linear', 0);
fprintf('\n[Contrast]  wall %.1f dB | inside box %.1f dB | outside %.1f dB\n', ...
    20*log10(vWall+eps), 20*log10(vIn+eps), 20*log10(vOut+eps));
fprintf('  -> wall-to-interior contrast = %.1f dB (v1 แบบ max+sinc ต่ำกว่านี้มาก)\n', ...
    20*log10(vWall/(vIn+eps)+eps));

[e3r,e3c,e3z] = peakExtents(magComb, xScene,yScene,zScene);
fprintf('\n[-3 dB extent through peak]  range %.2f | cross %.2f | height %.2f m\n', e3r,e3c,e3z);
fprintf('\n[Note] Fixes v2: Taylor -35 dB (range sidelobe), multilook แทน max\n');
fprintf('       (speckle floor), จุดถี่ขึ้น 400/หน้า (ผนังต่อเนื่อง)\n');

%% ===== SAVE FIGURES =====
figDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end
exportgraphics(figure(1), fullfile(figDir,'fig_w9v2_arc_vs_full.png'),'Resolution',150);
exportgraphics(figure(2), fullfile(figDir,'fig_w9v2_iso4views.png'),  'Resolution',150);
exportgraphics(figure(3), fullfile(figDir,'fig_w9v2_pointcloud.png'), 'Resolution',150);
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
% frequency-domain forward model (โครงเดียวกับ GOTCHA phase history):
%   fp(f,pulse) = w(f) · Σ_k a_k · exp(-j4πf R_k /c)  (demodulate ด้วย rref)
%   คืน rc = range-compressed profiles (Nzp x Naz) บนแกน differential range
    Naz  = size(plat,2);
    rref = sqrt(sum((plat - sceneCtr).^2, 1));               % 1 x Naz
    dR   = sqrt( (pos(1,:)'-plat(1,:)).^2 + (pos(2,:)'-plat(2,:)).^2 ...
               + (pos(3,:)'-plat(3,:)).^2 ) - rref;          % Nscat x Naz
    Nf = numel(fgrid);
    S  = zeros(Nf, Naz);
    for n = 1:Nf                                             % vectorized ต่อ freq
        S(n,:) = wFreq(n) * (amp * exp(-1j*4*pi*fgrid(n)/c * dR));
    end
    rc = fftshift(ifft(S, Nzp, 1), 1);
end


function img = bpOnePassDiff(rc, rref, plat, drAxis, voxT, kc, sz)
% 3D BP บนแกน differential range (แบบ gotcha_BP_reference)
    Nzp = numel(drAxis); ddr = drAxis(2)-drAxis(1); d1 = drAxis(1);
    Nvox = size(voxT,1);
    acc  = zeros(Nvox,1);
    for a = 1:size(plat,2)
        d  = voxT - plat(:,a).';
        dv = sqrt(sum(d.^2, 2)) - rref(a);       % differential range ของ voxel
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
% 12 ขอบของกล่องเป็น polyline เดียว (คั่นด้วย NaN) -> วาดด้วย plot3 หมุนไม่บั๊ก
%   แกน: X=cross, Y=range, Z=height ; hl=range half, hw=cross half
    cx=bc(2); cy=bc(1);
    C=[cx-hw,cy-hl,0; cx+hw,cy-hl,0; cx+hw,cy+hl,0; cx-hw,cy+hl,0; ...   % ล่าง 1-4
       cx-hw,cy-hl,H; cx+hw,cy-hl,H; cx+hw,cy+hl,H; cx-hw,cy+hl,H];      % บน  5-8
    E=[1 2;2 3;3 4;4 1; 5 6;6 7;7 8;8 5; 1 5;2 6;3 7;4 8];              % 12 ขอบ
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
