%% W7 v2 — 3D BACKPROJECTION (voxel grid) + IMPROVED VISUALIZATION
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  แก้จาก W7_3DBackprojection.m : ภาพ isosurface เดิมเป็น "แท่งลายทางแนวตั้ง"
%  ดูไม่ออกว่าเป็นกล่อง สาเหตุ + วิธีแก้ในไฟล์นี้:
%
%  (1) ลายทาง = carrier fringe ตามแกน z (คาบ ~lambda/2/sin(elev) ~ 8 cm)
%      แต่ voxel 1 m -> aliasing เป็นลายโมเร่
%      แก้: smooth3 (Gaussian, เน้นแกน z) กับ |image| = incoherent averaging
%      เพื่อการแสดงผล (ค่า evaluation ยังใช้ข้อมูลดิบเหมือนเดิม)
%  (2) แกนภาพไม่ equal -> กล่องสัดส่วนเพี้ยน
%      แก้: daspect([1 1 1]) ทุก axes 3D
%  (3) มุมมองเดียว + isosurface ชั้นเดียวทึบ -> เห็นแค่บางด้าน
%      แก้: 4 มุมมอง (perspective/หลังกล่อง/ด้านหน้า/top-down),
%      isosurface 2 ชั้นโปร่งแสง (-6 dB ทึบ, -12 dB โปร่ง),
%      กล่องจริงวาดเป็นหน้าโปร่งแสง + เพิ่ม Fig 4 point cloud สีตาม dB
%
%  หมายเหตุฟิสิกส์: single pass ไม่มี elevation aperture -> ความสูงยัง
%  ไม่ resolve (พลังงาน smear ตาม z) ภาพจะเป็น "กำแพง/แผ่น" ตรงตำแหน่ง
%  หน้ากล่องที่มองเห็น ไม่ใช่กล่องปิด 6 ด้าน — อันนั้นต้องรอ W8 (multi-pass)
%
%  Section 1-4 (config / box / echo / BP) เหมือน W7 เดิมทุกประการ

clear; clc; close all;

%% ===== RADAR CONFIGURATION (เดียวกับ W6/W7) =====

c  = physconst('LightSpeed');
fc = 4e9;
rangeResolution = 3;
bw = c / (2*rangeResolution);

prf      = 1000;
aperture = 4;
tpd      = 3e-6;
fs       = 120e6;

waveform = phased.LinearFMWaveform( ...
    'SampleRate', fs, 'PulseWidth', tpd, ...
    'PRF', prf, 'SweepBandwidth', bw);

speed          = 100;
flightDuration = 4;

slowTime  = 1/prf;
numpulses = flightDuration/slowTime + 1;

maxRange          = 2500;
truncrangesamples = ceil((2*maxRange/c)*fs);

antenna     = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);
antennaGain = aperture2gain(aperture, c/fc);
transmitter = phased.Transmitter('PeakPower', 50e3, 'Gain', antennaGain);
radiator    = phased.Radiator('Sensor', antenna, 'OperatingFrequency', fc, 'PropagationSpeed', c);
collector   = phased.Collector('Sensor', antenna, 'PropagationSpeed', c, 'OperatingFrequency', fc);
receiver    = phased.ReceiverPreamp('SampleRate', fs, 'NoiseFigure', 30);
channel     = phased.FreeSpace('PropagationSpeed', c, 'OperatingFrequency', fc, ...
                               'SampleRate', fs, 'TwoWayPropagation', true);

%% ===== 3D BOX DEFINITION =====

boxCenter = [1000; 0; 0];        % [range; cross; height]
boxL = 30;  halfL = boxL/2;      % along range
boxW = 20;  halfW = boxW/2;      % along cross
boxH = 10;  halfH = boxH/2;      % along height (box spans z = 0..10)

%% ===== PHYSICS-BASED RCS (flat plate, sub-patch) =====
lambda = c / fc;
N_face = 6;
ptsPerFace = N_face^2;

A_front_back = boxW * boxH;
A_left_right = boxL * boxH;
A_top_bottom = boxL * boxW;

splate_fb = 4*pi * A_front_back^2 / lambda^2;
splate_lr = 4*pi * A_left_right^2 / lambda^2;
splate_tb = 4*pi * A_top_bottom^2 / lambda^2;
sigmaPlate = [splate_fb, splate_fb, splate_lr, splate_lr, splate_tb, splate_tb];
faceRCS    = sigmaPlate / ptsPerFace^2;

faceLabels = {'Front','Back','Left','Right','Top','Bottom'};

%% ===== SCATTERER POSITIONS (6 faces) =====
u = linspace(-1, 1, N_face);
scatPos = [];  scatRCS = [];

[V, W] = meshgrid(u*halfW, u*halfH);
% Front (x = center - halfL)
xF = (boxCenter(1)-halfL)*ones(size(V));
scatPos = [scatPos, [xF(:)'; (boxCenter(2)+V(:))'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(1)*ones(1,N_face^2)];
% Back
xF = (boxCenter(1)+halfL)*ones(size(V));
scatPos = [scatPos, [xF(:)'; (boxCenter(2)+V(:))'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(2)*ones(1,N_face^2)];

[U, W] = meshgrid(u*halfL, u*halfH);
% Left
yF = (boxCenter(2)-halfW)*ones(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; yF(:)'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(3)*ones(1,N_face^2)];
% Right
yF = (boxCenter(2)+halfW)*ones(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; yF(:)'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(4)*ones(1,N_face^2)];

[U, V] = meshgrid(u*halfL, u*halfW);
% Top (z = boxH)
zF = (boxCenter(3)+boxH)*ones(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; (boxCenter(2)+V(:))'; zF(:)']];
scatRCS = [scatRCS, faceRCS(5)*ones(1,N_face^2)];
% Bottom (z = 0)
zF = zeros(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; (boxCenter(2)+V(:))'; zF(:)']];
scatRCS = [scatRCS, faceRCS(6)*ones(1,N_face^2)];

numScatterers = size(scatPos, 2);

%% ===== PLATFORM HISTORY =====
radarpos0 = [0; -200; 500];
radarvel0 = [0; speed; 0];
timeAxis  = ((1:numpulses)-1) / prf;
radarPosHistory = radarpos0 + radarvel0 * timeAxis;

%% ===== STATIC OCCLUSION =====
faceNormal = [ -1 0 0;  1 0 0;  0 -1 0;  0 1 0;  0 0 1;  0 0 -1];
losAll = radarPosHistory - boxCenter;
faceVisible = true(1,6);
for f = 1:6
    faceVisible(f) = any(faceNormal(f,:) * losAll > 0);
end
for f = 1:6
    idx = (f-1)*ptsPerFace + (1:ptsPerFace);
    scatRCS(idx) = scatRCS(idx) * faceVisible(f);
end
fprintf('Occlusion ON | visible: %s | hidden: %s\n', ...
    strjoin(faceLabels(faceVisible), ', '), strjoin(faceLabels(~faceVisible), ', '));

%% ===== SIMULATE RAW ECHO =====
fprintf('Simulating raw echo (%d pulses, %d scatterers)...\n', numpulses, numScatterers);
rxsig = zeros(truncrangesamples, numpulses);
boxPlatform = phased.Platform('InitialPosition', scatPos, 'Velocity', zeros(size(scatPos)));
boxTarget   = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', scatRCS);
radarPlatform = phased.Platform('InitialPosition', radarpos0, 'Velocity', radarvel0);

for ii = 1:numpulses
    [radarpos, radarvel] = radarPlatform(slowTime);
    [tpos, tvel]         = boxPlatform(slowTime);
    [~, targetAngle]     = rangeangle(tpos, radarpos);
    sig = waveform();
    sig = sig(1:truncrangesamples);
    sig = transmitter(sig);
    targetAngle(1,:) = zeros(1, numScatterers);
    sig = radiator(sig, targetAngle);
    sig = channel(sig, radarpos, tpos, radarvel, tvel);
    sig = boxTarget(sig);
    sig = collector(sig, targetAngle);
    rxsig(:,ii) = receiver(sig);
end
fprintf('Raw echo done.\n');

%% ===== RANGE COMPRESSION =====
refPulse     = waveform();
chirpSamples = round(tpd * fs);
refChirp     = refPulse(1:chirpSamples);
matchedFilter = conj(flipud(refChirp));
mfLen = length(matchedFilter);

rxsigRC = zeros(size(rxsig,1)+mfLen-1, numpulses);
for ii = 1:numpulses
    rxsigRC(:,ii) = conv(rxsig(:,ii), matchedFilter, 'full');
end
fprintf('Range compression done.\n');

%% ===== VOXEL GRID =====
gridMode = 'coarse';            % 'coarse' | 'fine'
switch gridMode
    case 'fine',   voxStep = 0.5;
    otherwise,     voxStep = 1.0;
end

xScene = (boxCenter(2)-30) : voxStep : (boxCenter(2)+30);   % cross  [-30..30]
yScene = (boxCenter(1)-20) : voxStep : (boxCenter(1)+20);   % range  [980..1020]
zScene = -5               : voxStep : 25;                   % height [-5..25]

[Xg, Yg, Zg] = meshgrid(xScene, yScene, zScene);
Ny = numel(yScene);  Nx = numel(xScene);  Nz = numel(zScene);
fprintf('Voxel grid (%s): %d x %d x %d = %d voxels (step %.2f m)\n', ...
    gridMode, Ny, Nx, Nz, Ny*Nx*Nz, voxStep);

%% ===== 3D BACKPROJECTION =====
fprintf('Running 3D backprojection...\n');
tic;
bp3 = backproject3D(rxsigRC, radarPosHistory, Xg, Yg, Zg, c, fs, fc, mfLen);
tBP = toc;
fprintf('3D backprojection complete (%.1f s).\n', tBP);

bp3Mag = abs(bp3);
bp3Mag = bp3Mag ./ max(bp3Mag(:));
bp3dB  = 20*log10(bp3Mag + eps);

%% ===== POST-PROCESSING FOR DISPLAY  << ใหม่ (v2) =====
%  ลบ carrier fringe ที่ alias ตามแกน z ด้วย incoherent (magnitude) smoothing
%  kernel [3 3 9]: เบาในระนาบ range-cross, หนักตามแกน z (แกนที่มีลายทาง)
bpSm = smooth3(bp3Mag, 'gaussian', [3 3 9], 2);
bpSm = bpSm ./ max(bpSm(:));
bpSmdB = 20*log10(bpSm + eps);

% true box geometry (ใช้ overlay ทุกรูป)
[Vb, Fb] = makeBoxWire(boxCenter, halfL, halfW, boxH);
[Xe, Ye, Ze] = boxEdges(boxCenter, halfL, halfW, boxH);   % กรอบเส้น (หมุนไม่บั๊ก)
boxCr = boxCenter(2) + halfW*[-1 1 1 -1 -1];
boxRg = boxCenter(1) + halfL*[-1 -1 1 1 -1];

%% ===== FIGURE 1: three MIP planes (เหมือนเดิม) =====
mipXY = squeeze(max(bp3Mag, [], 3));
mipRZ = squeeze(max(bp3Mag, [], 2));
mipXZ = squeeze(max(bp3Mag, [], 1));

figure(1); set(gcf,'Name','Fig1: 3D BP — MIP planes','Color','k','Position',[60 80 1250 420]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
imagesc(xScene, yScene, 20*log10(mipXY+eps)); set(ax1,'YDir','normal'); clim([-40 0]); colormap(ax1,'jet');
hold on; plot(boxCr, boxRg, 'w--', 'LineWidth', 1.5); hold off;
xlabel('Cross-Range (m)','Color','w'); ylabel('Range (m)','Color','w');
title('MIP top-down (collapse height)','Color','w');
set(ax1,'Color','k','XColor','w','YColor','w');

ax2 = nexttile;
imagesc(zScene, yScene, 20*log10(mipRZ+eps)); set(ax2,'YDir','normal'); clim([-40 0]); colormap(ax2,'jet');
hold on;
yline(boxCenter(1)-halfL,'w--','LineWidth',1); yline(boxCenter(1)+halfL,'w--','LineWidth',1);
xline(0,'w:','LineWidth',1); xline(boxH,'w:','LineWidth',1);
hold off;
xlabel('Height (m)','Color','w'); ylabel('Range (m)','Color','w');
title('MIP range-height (collapse cross)','Color','w');
set(ax2,'Color','k','XColor','w','YColor','w');

ax3 = nexttile;
imagesc(zScene, xScene, 20*log10(mipXZ+eps)); set(ax3,'YDir','normal'); clim([-40 0]); colormap(ax3,'jet');
hold on;
yline(-halfW,'w--','LineWidth',1); yline(halfW,'w--','LineWidth',1);
xline(0,'w:','LineWidth',1); xline(boxH,'w:','LineWidth',1);
hold off;
xlabel('Height (m)','Color','w'); ylabel('Cross-Range (m)','Color','w');
title('MIP cross-height (collapse range)','Color','w');
set(ax3,'Color','k','XColor','w','YColor','w');
cb = colorbar; cb.Color='w'; ylabel(cb,'dB','Color','w'); cb.Layout.Tile='east';

sgtitle('Fig 1 | 3D BP — MIP onto 3 orthogonal planes (single pass) [W7 v2]','Color','w');

%% ===== FIGURE 2: isosurface 4 มุมมอง  << แก้ใหม่ทั้งหมด (v2) =====
%  - ใช้ volume ที่ smooth แล้ว -> ผิวเรียบ ไม่เป็นลายทาง
%  - isosurface 2 ชั้น: -6 dB (ทึบ) + -12 dB (โปร่งแสง) เห็นทั้ง core และขอบ
%  - daspect([1 1 1]) -> สัดส่วนจริง
%  - 4 มุมมองในรูปเดียว + หมุนได้ (rotate3d)
figure(2); set(gcf,'Name','Fig2: 3D BP — isosurface 4 views (v2)','Color','w','Position',[60 60 1150 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

isoLv  = [0.5, 0.25];                      % -6 dB, -12 dB ของ peak (หลัง smooth)
isoCol = [0.85 0.25 0.15; 1.00 0.65 0.20];
isoAl  = [0.85, 0.25];

viewAz = [-35, 145,  90, -35];
viewEl = [ 22,  22,   0,  85];
viewNm = {'Perspective (หน้ากล่อง)','Perspective (หลังกล่อง)', ...
          'Front view (มองตามแนว range)','Top-down'};

for k = 1:4
    ax = nexttile; hold(ax,'on');
    for L = 1:2
        p = patch(ax, isosurface(Xg, Yg, Zg, bpSm, isoLv(L)));
        isonormals(Xg, Yg, Zg, bpSm, p);
        set(p,'FaceColor',isoCol(L,:),'EdgeColor','none','FaceAlpha',isoAl(L));
    end
    % true box: วาดเป็นเส้น 12 ขอบ (หมุนไม่บั๊ก — เลี่ยง depth-sort ของหน้าโปร่งแสง)
    plot3(ax, Xe, Ye, Ze, '-', 'Color',[0.1 0.3 0.9], 'LineWidth',1.5);
    hold(ax,'off');
    xlabel('Cross (m)'); ylabel('Range (m)'); zlabel('Height (m)');
    title(viewNm{k});
    xlim([min(xScene) max(xScene)]); ylim([min(yScene) max(yScene)]); zlim([min(zScene) max(zScene)]);
    daspect(ax,[1 1 1]);                   % << สำคัญ: สัดส่วนจริง
    grid(ax,'on'); box(ax,'on');
    view(ax, viewAz(k), viewEl(k));
    camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
end
sgtitle(sprintf('Fig 2 | 3D BP isosurface (smoothed) @ %.0f/%.0f dB + true box (blue) [W7 v2]', ...
    20*log10(isoLv(1)), 20*log10(isoLv(2))));
rotate3d on;

%% ===== FIGURE 3: point cloud สีตามความแรง  << ใหม่ (v2) =====
%  อีกวิธีที่อ่านง่าย: โชว์เฉพาะ voxel ที่แรงกว่า threshold เป็นจุดสี dB
pcThr = -15;                               % dB
idxPC = find(bpSmdB >= pcThr);
figure(3); set(gcf,'Name','Fig3: 3D BP — point cloud (v2)','Color','w','Position',[120 80 780 640]);
scatter3(Xg(idxPC), Yg(idxPC), Zg(idxPC), 10, bpSmdB(idxPC), 'filled', ...
    'MarkerFaceAlpha', 0.35);
hold on;
plot3(Xe, Ye, Ze, '-', 'Color',[0.1 0.3 0.9], 'LineWidth',1.8);   % กรอบเส้น
hold off;
colormap(jet); clim([pcThr 0]);
cb3 = colorbar; ylabel(cb3,'dB');
xlabel('Cross (m)'); ylabel('Range (m)'); zlabel('Height (m)');
xlim([min(xScene) max(xScene)]); ylim([min(yScene) max(yScene)]); zlim([min(zScene) max(zScene)]);
daspect([1 1 1]); grid on; box on; view(-35, 22);
title(sprintf('Fig 3 | Voxels \\geq %d dB + true box [W7 v2]', pcThr));
rotate3d on;

%% ===== EVALUATION (ใช้ข้อมูลดิบ ไม่ใช่ตัว smooth) =====
fprintf('\n========== EVALUATION — W7 v2: 3D Backprojection (single pass) ==========\n');

[~, pk] = max(bp3Mag(:));
[pyi, pxi, pzi] = ind2sub([Ny Nx Nz], pk);
fprintf('\n[Peak voxel]\n');
fprintf('  location  : range=%.1f m, cross=%.1f m, height=%.1f m\n', ...
    yScene(pyi), xScene(pxi), zScene(pzi));
fprintf('  box centre: range=%.0f m, cross=%.0f m, box height 0..%.0f m\n', ...
    boxCenter(1), boxCenter(2), boxH);

rCut = 20*log10(squeeze(bp3Mag(:, pxi, pzi)) + eps);  rCut = rCut - max(rCut);
cCut = 20*log10(squeeze(bp3Mag(pyi, :, pzi)) + eps);  cCut = cCut - max(cCut);
zCut = 20*log10(squeeze(bp3Mag(pyi, pxi, :)) + eps);  zCut = zCut - max(zCut);
rExt = axisExtent(yScene, rCut, -3);
cExt = axisExtent(xScene, cCut, -3);
zExt = axisExtent(zScene, zCut, -3);

fprintf('\n[-3 dB extent through peak]\n');
fprintf('  Range  : %.2f m   (box length %.0f m)\n', rExt, boxL);
fprintf('  Cross  : %.2f m   (box width  %.0f m)\n', cExt, boxW);
fprintf('  Height : %.2f m   (box height %.0f m)\n', zExt, boxH);
fprintf('\n[Note]\n');
fprintf('  v2 = display fix เท่านั้น (smoothing/aspect/views) ฟิสิกส์เท่าเดิม:\n');
fprintf('  single pass ไม่มี elevation aperture -> ความสูงยัง smear ตาม z\n');
fprintf('  ภาพจึงเป็น "แผ่น/กำแพง" ตรงหน้ากล่องที่มองเห็น ไม่ใช่กล่องปิด\n');
fprintf('  -> resolve ความสูง 10 m ต้องใช้ multi-pass tomography (W8).\n');

%% ===== SAVE FIGURES =====
figDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end
exportgraphics(figure(1), fullfile(figDir,'fig_w7v2_mip.png'),        'Resolution',150);
exportgraphics(figure(2), fullfile(figDir,'fig_w7v2_iso4views.png'),  'Resolution',150);
exportgraphics(figure(3), fullfile(figDir,'fig_w7v2_pointcloud.png'), 'Resolution',150);
fprintf('\nFigures saved to %s\n', figDir);


%%%% ===== LOCAL FUNCTIONS =====

function bp3 = backproject3D(rc, rph, Xg, Yg, Zg, c, fs, fc, mfLen)
% 3D coherent backprojection onto a voxel grid.
    sz  = size(Xg);
    vox = [Yg(:)'; Xg(:)'; Zg(:)'];
    bp3 = zeros(sz);
    np  = size(rph, 2);
    Nrc = size(rc, 1);
    for ii = 1:np
        rp = rph(:, ii);
        dv = vox - rp;
        sr = sqrt(sum(dv.^2, 1));
        si = 2*sr/c*fs + mfLen;
        vm = (si >= 1) & (si <= Nrc);
        iv = zeros(1, size(vox,2));
        if any(vm)
            iv(vm) = interp1(1:Nrc, rc(:,ii), si(vm), 'linear', 0);
        end
        bp3 = bp3 + reshape(iv .* exp(1j*4*pi*fc/c .* sr), sz);
    end
end


function [V, F] = makeBoxWire(bc, hl, hw, H)
% 8 vertices / 6 faces of the true box (X=cross, Y=range, Z=height).
    cx = bc(2); cy = bc(1);
    V = [cx-hw, cy-hl, 0;
         cx+hw, cy-hl, 0;
         cx+hw, cy+hl, 0;
         cx-hw, cy+hl, 0;
         cx-hw, cy-hl, H;
         cx+hw, cy-hl, H;
         cx+hw, cy+hl, H;
         cx-hw, cy+hl, H];
    F = [1 2 3 4;
         5 6 7 8;
         1 2 6 5;
         3 4 8 7;
         1 4 8 5;
         2 3 7 6];
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


function ext = axisExtent(ax, pdB, lvl)
% width of the region >= lvl (dB) around the peak of a 1-D cut.
    ax = ax(:); pdB = pdB(:);
    a = ax(pdB >= lvl);
    if isempty(a); ext = 0; else; ext = max(a) - min(a); end
end
