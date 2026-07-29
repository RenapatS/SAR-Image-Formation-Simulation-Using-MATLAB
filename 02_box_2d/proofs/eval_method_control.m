%% SAR Backprojection — CONTROL: Week 5 RCS (manual) + Week 6 NEW evaluation
%  หมายเหตุ: บล็อก RADAR CONFIGURATION สืบทอดมาจาก 01_point_target/point3_bp_2d.m
%  ซึ่งดัดแปลงจากตัวอย่าง Stripmap SAR ของ MathWorks (ดู attribution ในไฟล์นั้น)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%
%  จุดประสงค์: พิสูจน์ว่าตัวเลข W5 ที่ "ดูเพี้ยน" เกิดจากวิธีวัดผล ไม่ใช่ RCS
%    - ใช้ RCS แบบ Week 5 (ตั้งมือ): faceRCS = [2.0, 0.5, 1.0, 1.0, 0.8, 0.3]
%    - แต่ใช้ evaluation แบบ Week 6: point-target calibration + box extent
%
%  ผลที่คาดไว้:
%    - Calibration (resolution/PSLR) = เท่ากับ try104 เป๊ะ
%      (calibration ใช้ point target RCS=1 แยกต่างหาก ไม่ขึ้นกับ RCS กล่อง)
%    - Box extent / localization / dynamic range = กลับไปเท่าค่า W5 เดิม
%      (~7.88 m, loc err ~17.62 m, ~96.2 dB) เพราะขึ้นกับน้ำหนัก RCS
%  => ยืนยันว่าค่า W5 เดิม = target EXTENT ไม่ใช่ resolution

clear; clc; close all;

%% ===== RADAR CONFIGURATION =====

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

radarPlatform = phased.Platform( ...
    'InitialPosition', [0; -200; 500], 'Velocity', [0; speed; 0]);

slowTime   = 1/prf;
numpulses  = flightDuration/slowTime + 1;

maxRange          = 2500;
truncrangesamples = ceil((2*maxRange/c)*fs);
fastTime          = (0:1/fs:(truncrangesamples-1)/fs);

antenna     = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);
antennaGain = aperture2gain(aperture, c/fc);
transmitter = phased.Transmitter('PeakPower', 50e3, 'Gain', antennaGain);
radiator    = phased.Radiator('Sensor', antenna, 'OperatingFrequency', fc, 'PropagationSpeed', c);
collector   = phased.Collector('Sensor', antenna, 'PropagationSpeed', c, 'OperatingFrequency', fc);
receiver    = phased.ReceiverPreamp('SampleRate', fs, 'NoiseFigure', 30);
channel     = phased.FreeSpace('PropagationSpeed', c, 'OperatingFrequency', fc, ...
                               'SampleRate', fs, 'TwoWayPropagation', true);

%% ===== 3D BOX DEFINITION =====

boxCenter = [1000; 0; 0];
boxL = 30;   halfL = boxL/2;
boxW = 20;   halfW = boxW/2;
boxH = 10;   halfH = boxH/2;

%% ===== RCS: WEEK 5 MANUAL VALUES (control) =====
%  ใช้ค่า RCS แบบ Week 5 (try103) ตรงๆ ไม่ใช้สูตรฟิสิกส์
%  เพื่อแยกผลของ "RCS" ออกจากผลของ "วิธี evaluation"

lambda_calc = c / fc;
ptsPerFace  = 36;

% [Front, Back, Left, Right, Top, Bottom] — Week 5 hand-picked values
faceRCS    = [2.0, 0.5, 1.0, 1.0, 0.8, 0.3];
sigmaPlate = faceRCS;            % no plate formula here; keep name for eval compatibility

fprintf('===== CONTROL: Week 5 manual RCS + Week 6 evaluation =====\n');
fprintf('  faceRCS (manual) = [%.1f %.1f %.1f %.1f %.1f %.1f] m^2\n', faceRCS);
fprintf('  Evaluation       = point-target calibration + box extent (same as try104)\n');
fprintf('  Expectation: calibration == try104; box extent/loc == old W5 numbers\n\n');

%% ===== SCATTERER POSITIONS =====

N_face = 6;
u = linspace(-1, 1, N_face);

scatPos = [];
scatRCS = [];
faceColors = [1 0.2 0.2; 0.2 0.4 1; 0.2 0.8 0.2; 0.8 0.2 0.8; 0.2 0.8 0.8; 1 0.8 0.2];
faceLabels = {'Front','Back','Left','Right','Top','Bottom'};

[V, W] = meshgrid(u*halfW, u*halfH);

% Face 1: FRONT  (x = center − halfL, facing radar approaching from −y)
xF = (boxCenter(1)-halfL)*ones(size(V));
scatPos = [scatPos, [xF(:)'; (boxCenter(2)+V(:))'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(1)*ones(1,N_face^2)];

% Face 2: BACK
xF = (boxCenter(1)+halfL)*ones(size(V));
scatPos = [scatPos, [xF(:)'; (boxCenter(2)+V(:))'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(2)*ones(1,N_face^2)];

[U, W] = meshgrid(u*halfL, u*halfH);

% Face 3: LEFT
yF = (boxCenter(2)-halfW)*ones(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; yF(:)'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(3)*ones(1,N_face^2)];

% Face 4: RIGHT
yF = (boxCenter(2)+halfW)*ones(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; yF(:)'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(4)*ones(1,N_face^2)];

[U, V] = meshgrid(u*halfL, u*halfW);

% Face 5: TOP
zF = (boxCenter(3)+boxH)*ones(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; (boxCenter(2)+V(:))'; zF(:)']];
scatRCS = [scatRCS, faceRCS(5)*ones(1,N_face^2)];

% Face 6: BOTTOM
zF = zeros(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; (boxCenter(2)+V(:))'; zF(:)']];
scatRCS = [scatRCS, faceRCS(6)*ones(1,N_face^2)];

numScatterers = size(scatPos, 2);
ptsPerFace    = N_face^2;
fprintf('Box: %d scatterers total\n', numScatterers);

%% ===== PLATFORM HISTORY =====

radarpos0 = [0; -200; 500];
radarvel0 = [0; speed; 0];
pulseIndex = 1:numpulses;
timeAxis   = (pulseIndex-1) / prf;
radarPosHistory = radarpos0 + radarvel0 * timeAxis;

%% ===== SIMULATE RAW ECHO =====

fprintf('Simulating raw echo...\n');
rxsig = zeros(truncrangesamples, numpulses);

scatVel      = zeros(size(scatPos));
boxPlatform2 = phased.Platform('InitialPosition', scatPos, 'Velocity', scatVel);
boxTarget    = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', scatRCS);

radarPlatform2 = phased.Platform('InitialPosition', [0;-200;500], 'Velocity', [0;speed;0]);

for ii = 1:numpulses
    [radarpos, radarvel] = radarPlatform2(slowTime);
    [tpos, tvel]         = boxPlatform2(slowTime);
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

rangeAxis = fastTime * c / 2;

%% ===== RANGE COMPRESSION =====

refPulse      = waveform();
chirpSamples  = round(tpd * fs);
refChirp      = refPulse(1:chirpSamples);
matchedFilter = conj(flipud(refChirp));

rxsigRC = zeros(size(rxsig,1)+length(matchedFilter)-1, numpulses);
for ii = 1:numpulses
    rxsigRC(:,ii) = conv(rxsig(:,ii), matchedFilter, 'full');
end

rangeAxisRC = ((0:size(rxsigRC,1)-1).' - (length(matchedFilter)-1)) / fs * c/2;
rxsigRCmag  = abs(rxsigRC);
rxsigRCmag  = rxsigRCmag ./ max(rxsigRCmag(:));

%% ===== BACKPROJECTION =====

fprintf('Running Backprojection...\n');

xScene = linspace(-100, 100, 401);
yScene = linspace(850,  1200, 401);
[xGrid, yGrid] = meshgrid(xScene, yScene);
Ny = length(yScene);  Nx = length(xScene);
pixelPos = [yGrid(:)'; xGrid(:)'; zeros(1,numel(xGrid))];

bpImage = zeros(Ny, Nx);
for ii = 1:numpulses
    rp         = radarPosHistory(:, ii);
    diffVec    = pixelPos - rp;
    slantRange = sqrt(sum(diffVec.^2, 1));
    sampleIdx  = 2 .* slantRange ./ c .* fs + length(matchedFilter);
    validMask  = (sampleIdx >= 1) & (sampleIdx <= size(rxsigRC,1));
    interpVals = zeros(1, numel(xGrid));
    if any(validMask)
        interpVals(validMask) = interp1(1:size(rxsigRC,1), rxsigRC(:,ii), ...
                                        sampleIdx(validMask), 'linear', 0);
    end
    phaseCorr = exp(1j*4*pi*fc/c .* slantRange);
    bpImage   = bpImage + reshape(interpVals .* phaseCorr, Ny, Nx);
end
fprintf('Backprojection complete.\n');

bpImageMag = abs(bpImage);
bpImageMag = bpImageMag ./ max(bpImageMag(:));
bpImagedB  = 20*log10(bpImageMag + eps);

% box outline
boxX_cr = boxCenter(2) + halfW*[-1 1 1 -1 -1];
boxY_r  = boxCenter(1) + halfL*[-1 -1 1 1 -1];

%% ===== POINT-TARGET CALIBRATION — TRUE system resolution / PSLR =====
%  resolution = ความสามารถแยกสองจุดที่อยู่ใกล้กัน -> ต้องวัดจาก impulse
%  response (IRF) ของ "จุดเดียว" ไม่ใช่จากกล่อง
%  (กล่องยาว 30 m ทำให้ -3 dB width = ขนาดเป้า ไม่ใช่ resolution ระบบ)
%  ใช้ forward model + BP ตัวเดียวกับกล่อง แต่ noise-free (ไม่ผ่าน receiver)
%  และ backproject บน grid ละเอียดรอบจุด เพื่อ sample IRF ให้พอ

fprintf('\nRunning point-target calibration...\n');

calPos       = [boxCenter(1); boxCenter(2); 0];   % single scatterer at scene center
calTarget    = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', 1);
calTx        = phased.Transmitter('PeakPower', 50e3, 'Gain', antennaGain);
calRadiator  = phased.Radiator('Sensor', antenna, 'OperatingFrequency', fc, 'PropagationSpeed', c);
calCollector = phased.Collector('Sensor', antenna, 'PropagationSpeed', c, 'OperatingFrequency', fc);
calChannel   = phased.FreeSpace('PropagationSpeed', c, 'OperatingFrequency', fc, ...
                                'SampleRate', fs, 'TwoWayPropagation', true);

rxsigCal = zeros(truncrangesamples, numpulses);
for ii = 1:numpulses
    rpos = radarPosHistory(:, ii);
    [~, tAng] = rangeangle(calPos, rpos);
    s = waveform(); s = s(1:truncrangesamples);
    s = calTx(s);
    tAng(1) = 0;                                   % mirror box model: boresight in az
    s = calRadiator(s, tAng);
    s = calChannel(s, rpos, calPos, [0;speed;0], [0;0;0]);
    s = calTarget(s);
    s = calCollector(s, tAng);
    rxsigCal(:, ii) = s;                           % no receiver -> clean, noise-free IRF
end

% range compression (reuse matched filter)
rxsigCalRC = zeros(size(rxsigCal,1)+length(matchedFilter)-1, numpulses);
for ii = 1:numpulses
    rxsigCalRC(:,ii) = conv(rxsigCal(:,ii), matchedFilter, 'full');
end

% backprojection on a FINE grid centered on the point (resolve sub-metre IRF)
xCal = linspace(-2, 2, 401);                              % cross-range, 0.01 m step
yCal = linspace(boxCenter(1)-10, boxCenter(1)+10, 401);   % range, 0.05 m step
[xGc, yGc] = meshgrid(xCal, yCal);
Nyc = numel(yCal);  Nxc = numel(xCal);
pixC = [yGc(:)'; xGc(:)'; zeros(1,numel(xGc))];
bpCal = zeros(Nyc, Nxc);
for ii = 1:numpulses
    rp = radarPosHistory(:, ii);
    dv = pixC - rp;
    sr = sqrt(sum(dv.^2, 1));
    si = 2 .* sr ./ c .* fs + length(matchedFilter);
    vm = (si >= 1) & (si <= size(rxsigCalRC,1));
    iv = zeros(1, numel(xGc));
    if any(vm)
        iv(vm) = interp1(1:size(rxsigCalRC,1), rxsigCalRC(:,ii), si(vm), 'linear', 0);
    end
    bpCal = bpCal + reshape(iv .* exp(1j*4*pi*fc/c .* sr), Nyc, Nxc);
end
bpCalMag = abs(bpCal);  bpCalMag = bpCalMag ./ max(bpCalMag(:));

% IRF cuts through the peak
[~, pkLin] = max(bpCalMag(:));
[pkR, pkC] = ind2sub([Nyc Nxc], pkLin);
calRangeCutdB = 20*log10(bpCalMag(:, pkC) + eps);  calRangeCutdB = calRangeCutdB - max(calRangeCutdB);
calCrossCutdB = 20*log10(bpCalMag(pkR, :) + eps);  calCrossCutdB = calCrossCutdB - max(calCrossCutdB);

% theoretical resolution (for comparison)
lambda      = c / fc;
apertureLen = speed * flightDuration;
rangeRes_th = c / (2 * bw);
crRes_th    = lambda * boxCenter(1) / (2 * apertureLen);

[calRangeRes, calPSLR_r] = irfMetrics(yCal(:), calRangeCutdB(:));
[calCrossRes, calPSLR_c] = irfMetrics(xCal(:), calCrossCutdB(:));

fprintf('Calibration done: range res = %.2f m, cross res = %.3f m, PSLR(r/c) = %.1f / %.1f dB\n', ...
    calRangeRes, calCrossRes, calPSLR_r, calPSLR_c);

%% ===== FIGURE 1: 3D Box Geometry =====

figure(1);
set(gcf,'Name','Fig1: 3D Box Geometry');

[Verts, Faces] = makeBoxPatch(boxCenter(2), boxCenter(1), boxCenter(3), halfW, halfL, boxH);
patch('Vertices', Verts, 'Faces', Faces, ...
      'FaceColor', [0.7 0.85 1.0], 'FaceAlpha', 0.15, ...
      'EdgeColor', [0.3 0.3 0.8], 'LineWidth', 1.5);
hold on;
for f = 1:6
    idx = (f-1)*ptsPerFace + (1:ptsPerFace);
    % marker size scaled by plate dBsm (physically meaningful face RCS)
    markerSz = 10*log10(sigmaPlate(f)) - 60;   % plate dBsm offset -> marker size
    scatter3(scatPos(2,idx), scatPos(1,idx), scatPos(3,idx), ...
        max(markerSz, 10), ...
        repmat(faceColors(f,:), ptsPerFace, 1), 'filled', ...
        'DisplayName', sprintf('%s (\\sigma_{plate}=%.1f dBsm)', faceLabels{f}, 10*log10(sigmaPlate(f))));
end
xlabel('Cross-Range (m)'); ylabel('Range (m)'); zlabel('Height (m)');
title('Fig 1 | 3D Box Target — W5 manual RCS (control)');
legend('Box (transparent)', faceLabels{:}, 'Location','northeast');
grid on; view(40, 25);
xlim([-60 60]); ylim([920 1080]); zlim([-5 35]);
hold off;

%% ===== FIGURE 2: SAR Acquisition Scene 3D =====

figure(2);
set(gcf,'Name','Fig2: SAR Acquisition Scene 3D');

bpImg2dB = max(bpImagedB, -30);
surf(xScene, yScene, zeros(Ny,Nx), bpImg2dB, 'EdgeAlpha',0, 'FaceAlpha',0.85);
colormap(gca,'jet'); hold on;

[Verts2, Faces2] = makeBoxPatch(boxCenter(2),boxCenter(1),boxCenter(3),halfW,halfL,boxH);
patch('Vertices',Verts2,'Faces',Faces2, ...
    'FaceColor',[0.3 0.6 1],'FaceAlpha',0.3,'EdgeColor',[1 0.9 0.1],'LineWidth',3, ...
    'DisplayName','Box');
scatter3(scatPos(2,:), scatPos(1,:), scatPos(3,:), ...
    15, 'w', 'filled', 'MarkerFaceAlpha',0.5, 'HandleVisibility','off');

rph2 = radarPosHistory;
plot3(rph2(2,1:50:end), rph2(1,1:50:end), rph2(3,1:50:end), ...
    'b-', 'LineWidth', 2.5, 'DisplayName','Radar path');
scatter3(rph2(2,1),   rph2(1,1),   rph2(3,1),   200, 'gs','filled','DisplayName','Start');
scatter3(rph2(2,end), rph2(1,end), rph2(3,end), 200, 'rs','filled','DisplayName','End');

midP2 = round(numpulses/2);
rp2   = rph2(:, midP2);
fpX2  = boxCenter(2) + halfW*2*[-1 1 1 -1];
fpY2  = boxCenter(1) + halfL*2*[-1 -1 1 1];
for cr = 1:4
    hb = plot3([rp2(2) fpX2(cr)],[rp2(1) fpY2(cr)],[rp2(3) 1],'c-','LineWidth',1.2);
    try hb.Color(4)=0.45; catch; hb.Color=[0 0.85 0.85]; end
end

cb2 = colorbar; clim([-30 0]); ylabel(cb2,'dB');
xlabel('Cross-Range (m)'); ylabel('Range (m)'); zlabel('Height (m)');
title('Fig 2 | SAR Scene — Radar Path + Beam + Box + BP Image (dB) [W6]');
legend({'BP Image','Box','Radar path','Start','End'},'Location','northwest');
grid on;
xlim([-150 150]); ylim([-280 1100]); zlim([-5 560]);
view(-20, 28); camlight left; lighting gouraud;
hold off;

%% ===== FIGURE 3: Raw Echo — imagesc (linear amplitude) =====

figure(3);
set(gcf,'Name','Fig3: Raw Echo imagesc (linear)');

inRange2  = (rangeAxis >= 1060) & (rangeAxis <= 1170);
pulseAxis = 1:numpulses;
radarY    = radarpos0(2) + radarvel0(2)*(pulseAxis-1)/prf;
expectedRange = sqrt((boxCenter(1)-radarpos0(1))^2 + (boxCenter(2)-radarY).^2 + radarpos0(3)^2);

imagesc(1:numpulses, rangeAxis(inRange2), abs(rxsig(inRange2, :)));
set(gca,'YDir','normal');
colormap(gca,'parula'); colorbar;
xlabel('Pulse Index (slow time)'); ylabel('Slant Range (m)');
title('Fig 3 | Raw Echo — Range × Pulse (linear amplitude) [W6]');
hold on;
plot(pulseAxis, expectedRange, 'r--', 'LineWidth', 1.5, 'DisplayName','Expected slant range');
ylim([1060 1170]);
legend('Location','northeast'); grid on;
hold off;

%% ===== FIGURE 4: Raw Echo — imagesc (dB scale) =====

figure(4);
set(gcf,'Name','Fig4: Raw Echo imagesc (dB)');

rawdB = 20*log10(abs(rxsig(inRange2, :)) + eps);
rawdB = rawdB - max(rawdB(:));
imagesc(1:numpulses, rangeAxis(inRange2), rawdB);
set(gca,'YDir','normal');
clim([-40 0]);
colormap(gca,'jet'); cb4 = colorbar; ylabel(cb4,'dB');
xlabel('Pulse Index (slow time)'); ylabel('Slant Range (m)');
title('Fig 4 | Raw Echo — Range × Pulse (dB, clip -40 dB) [W6]');
hold on;
plot(pulseAxis, expectedRange, 'w--', 'LineWidth', 1.5, 'DisplayName','Expected slant range');
ylim([1060 1170]);
legend('Location','northeast'); grid on;
hold off;

%% ===== FIGURE 5: Range Compressed — imagesc (dB) =====

figure(5);
set(gcf,'Name','Fig5: Range Compressed');

inR5    = (rangeAxisRC >= 1060) & (rangeAxisRC <= 1170);
step5   = 10;
pulses5 = 1:step5:numpulses;
imagesc(pulses5, rangeAxisRC(inR5), rxsigRCmag(inR5, pulses5));
set(gca,'YDir','normal');
colormap(gca,'jet'); cb5 = colorbar; ylabel(cb5,'Normalized Amp');
xlabel('Pulse Index'); ylabel('Slant Range (m)');
title('Fig 5 | Range-Compressed Data [W6]');
ylim([1060 1170]); grid on;

%% ===== FIGURE 6: BP Image — 3D Surface + 2D Views =====
%  NOTE FIX: Range Profile (bottom-right) — normalize to peak 0 dB + set xlim

figure(6);
set(gcf,'Name','Fig6: BP Image 3D+2D');
set(gcf,'Color','k');

bpClip = max(bpImagedB, -40);

tl = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');

% Left: 3D Surface (span 2 rows)
ax1 = nexttile([2 1]);
surf(xScene, yScene, bpClip, 'EdgeAlpha', 0);
shading interp; colormap(ax1,'jet'); clim([-40 0]);
hold on;
contour3(xScene, yScene, bpClip, 12, 'k', 'LineWidth', 0.4);
plot3(boxX_cr, boxY_r, 2*ones(1,5), 'y--', 'LineWidth', 2);
xlabel('Cross-Range (m)','Color','w'); ylabel('Range (m)','Color','w'); zlabel('dB','Color','w');
title('3D Surface (dB)', 'Color','w');
colorbar; view(-40, 30); grid on; camlight left; lighting gouraud;
set(ax1,'Color','k','XColor','w','YColor','w','ZColor','w'); hold off;

% Top-right: Top-down imagesc
ax2 = nexttile;
imagesc(xScene, yScene, bpClip); axis xy;
colormap(ax2,'jet'); clim([-40 0]);
hold on;
plot(boxX_cr, boxY_r, 'y--', 'LineWidth', 1.8);
hold off;
xlabel('Cross-Range (m)','Color','w'); ylabel('Range (m)','Color','w');
title('Top-down View', 'Color','w');
colorbar; set(ax2,'Color','k','XColor','w','YColor','w');

% Bottom-right: Range Profile — FIX: normalize + proper xlim
ax3 = nexttile;
[~, xMid] = min(abs(xScene - boxCenter(2)));
rangeProf    = bpImageMag(:, xMid);
rangeProfDB  = 20*log10(rangeProf + eps);
rangeProfDB  = rangeProfDB - max(rangeProfDB);   % normalize to 0 dB peak
plot(rangeProfDB, yScene, 'c-', 'LineWidth', 1.8); hold on;
xline(-3,  'r--', 'LineWidth', 1.2, 'DisplayName','-3 dB');
xline(-10, 'm--', 'LineWidth', 1.0, 'DisplayName','-10 dB');
% X axis = dB, Y axis = Range -> use yline for the face markers
yline(boxCenter(1)-halfL, 'y--', 'LineWidth',1,'DisplayName','Front face');
yline(boxCenter(1)+halfL, 'g--', 'LineWidth',1,'DisplayName','Back face');
xlabel('Normalized (dB)','Color','w'); ylabel('Range (m)','Color','w');
xlim([-40 5]);   % FIX: กำหนด dB scale ให้ถูกต้อง
ylim([850 1200]);
title('Range Profile (cross=0) — normalized', 'Color','w');
grid on; set(ax3,'Color','k','XColor','w','YColor','w'); hold off;

sgtitle('Fig 6 | BP Image [W6 — Physics RCS]', 'Color','w', 'FontSize',11);

%% ===== FIGURE 7: Cross-Range Profiles — 3D Filled Slices =====

figure(7);
set(gcf,'Name','Fig7: Cross-Range Profiles 3D');

rangeSlices  = [boxCenter(1)-halfL, boxCenter(1), boxCenter(1)+halfL];
sliceLabels  = {'Front face','Box center','Back face'};
sliceColors3 = {[1 0.2 0.2], [0.1 0.1 0.1], [0.2 0.7 0.2]};

hold on;
for k = 1:3
    [~, rIdx] = min(abs(yScene - rangeSlices(k)));
    prof = bpImageMag(rIdx, :);
    fill3(xScene([1:end end 1]), ...
          rangeSlices(k)*ones(1,Nx+2), ...
          [prof 0 0], ...
          sliceColors3{k}, 'FaceAlpha', 0.3, 'EdgeColor', sliceColors3{k}, ...
          'DisplayName', sliceLabels{k});
end
plot3([boxX_cr; boxX_cr], [boxY_r; boxY_r], ...
      [zeros(1,5); 0.8*ones(1,5)], 'k--', 'LineWidth', 1);
xlabel('Cross-Range (m)'); ylabel('Ground Range (m)'); zlabel('Amplitude');
title('Fig 7 | Cross-Range Profiles — 3D Filled Slices [W6]');
legend show; grid on; view(-30, 25);
hold off;

%% ===== FIGURE 8: SAR Image 2D — Top-down (dB) =====

figure(8);
set(gcf,'Name','Fig8: SAR Image 2D (dB)');

imagesc(xScene, yScene, bpImagedB);
set(gca,'YDir','normal');
clim([-40 0]);
colormap(gca,'jet');
cb8 = colorbar; ylabel(cb8,'dB');
hold on;
plot(boxX_cr, boxY_r, 'w--', 'LineWidth', 2, 'DisplayName','Box outline');
plot([0 0], [yScene(1) yScene(end)], 'w:', 'LineWidth', 1, 'HandleVisibility','off');
plot([xScene(1) xScene(end)], [boxCenter(1) boxCenter(1)], 'w:', 'LineWidth', 1, 'HandleVisibility','off');
xlabel('Cross-Range (m)'); ylabel('Ground Range (m)');
title('Fig 8 | SAR BP Image — Top-down 2D (dB) [W6]');
legend('Box outline','Location','northeast');
axis equal tight; grid on;
hold off;

%% ===== FIGURE 9: Range Profile 2D — dB =====

figure(9);
set(gcf,'Name','Fig9: Range Profile 2D (dB)');

[~, xZeroIdx9] = min(abs(xScene));
rProf9   = bpImageMag(:, xZeroIdx9);
rProf9dB = 20*log10(rProf9 + eps);
rProf9dB = rProf9dB - max(rProf9dB);

plot(yScene, rProf9dB, 'b-', 'LineWidth', 1.8); hold on;
yline(-3,  'r--', 'LineWidth', 1.5, 'DisplayName','-3 dB');
yline(-10, 'm--', 'LineWidth', 1,   'DisplayName','-10 dB');

faceR9 = [boxCenter(1)-halfL, boxCenter(1), boxCenter(1)+halfL];
fClr9  = {'r','k','g'};
fNm9   = {'Front (985m)','Center (1000m)','Back (1015m)'};
for k = 1:3
    xline(faceR9(k), fClr9{k}, fNm9{k}, 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom');
end

above3r9 = yScene(rProf9dB >= -3);
if ~isempty(above3r9)
    w3r9 = max(above3r9) - min(above3r9);
    xline(min(above3r9), 'r:', 'LineWidth', 1, 'HandleVisibility','off');
    xline(max(above3r9), 'r:', 'LineWidth', 1, 'HandleVisibility','off');
    title(sprintf('Fig 9 | Range Profile (Cross-Range=0) | -3dB = %.1f m  (theory ≈ %.0f m) [W6]', ...
        w3r9, rangeResolution));
else
    title('Fig 9 | Range Profile (Cross-Range=0) [W6]');
end

xlim([850 1200]); ylim([-40 5]);
xlabel('Ground Range (m)'); ylabel('Normalized (dB)');
legend('Location','northeast'); grid on;
hold off;

%% ===== FIGURE 10: Cross-Range Profile 2D — dB =====

figure(10);
set(gcf,'Name','Fig10: Cross-Range Profile 2D (dB)');

[~, rCenterIdx10] = min(abs(yScene - boxCenter(1)));
crProf10   = bpImageMag(rCenterIdx10, :);
crProf10dB = 20*log10(crProf10 + eps);
crProf10dB = crProf10dB - max(crProf10dB);

plot(xScene, crProf10dB, 'b-', 'LineWidth', 1.8); hold on;
yline(-3,  'r--', 'LineWidth', 1.5, 'DisplayName','-3 dB');
yline(-10, 'm--', 'LineWidth', 1,   'DisplayName','-10 dB');
xline(-halfW, 'k--', 'LineWidth', 1.5, 'DisplayName','Box edges (±10m)');
xline(+halfW, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');

above3cr10 = xScene(crProf10dB >= -3);
if ~isempty(above3cr10)
    w3cr10 = max(above3cr10) - min(above3cr10);
    xline(min(above3cr10), 'r:', 'LineWidth', 1, 'HandleVisibility','off');
    xline(max(above3cr10), 'r:', 'LineWidth', 1, 'HandleVisibility','off');
    lambda_ev = c/fc;
    apertureLen_ev = speed * flightDuration;
    crRes_th10 = lambda_ev * boxCenter(1) / (2 * apertureLen_ev);
    title(sprintf('Fig 10 | Cross-Range Profile | -3dB = %.1f m  (theory ≈ %.2f m) [W6]', ...
        w3cr10, crRes_th10));
else
    title('Fig 10 | Cross-Range Profile (Range=1000m) [W6]');
end

xlim([-100 100]); ylim([-40 5]);
xlabel('Cross-Range (m)'); ylabel('Normalized (dB)');
legend('Location','northeast'); grid on;
hold off;

%% ===== FIGURE 12: Point-Target Calibration — IRF (true resolution) =====

figure(12);
set(gcf,'Name','Fig12: Point-Target Calibration IRF');

subplot(1,2,1);
plot(yCal, calRangeCutdB, 'b-', 'LineWidth', 1.8); hold on;
yline(-3, 'r--', 'LineWidth', 1.2, 'DisplayName','-3 dB');
yline(calPSLR_r, 'm:', 'LineWidth', 1.2, 'DisplayName', sprintf('PSLR %.1f dB', calPSLR_r));
xlabel('Range (m)'); ylabel('Normalized (dB)');
title(sprintf('Range IRF | -3dB = %.2f m  (theory %.2f m)', calRangeRes, rangeRes_th));
xlim([boxCenter(1)-10 boxCenter(1)+10]); ylim([-40 3]); grid on;
legend('Location','northeast'); hold off;

subplot(1,2,2);
plot(xCal, calCrossCutdB, 'b-', 'LineWidth', 1.8); hold on;
yline(-3, 'r--', 'LineWidth', 1.2, 'DisplayName','-3 dB');
yline(calPSLR_c, 'm:', 'LineWidth', 1.2, 'DisplayName', sprintf('PSLR %.1f dB', calPSLR_c));
xlabel('Cross-Range (m)'); ylabel('Normalized (dB)');
title(sprintf('Cross IRF | -3dB = %.3f m  (theory %.4f m)', calCrossRes, crRes_th));
xlim([-2 2]); ylim([-40 3]); grid on;
legend('Location','northeast'); hold off;

sgtitle('Fig 12 | Point-Target Calibration — true system IRF (single scatterer, noise-free)');

%% ===== SAVE FIGURES =====

figDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'figure');
if ~exist(figDir, 'dir'); mkdir(figDir); end
for fn = [1:10, 12]   % Fig 11 (eval text) is built later and saved after creation
    fh = figure(fn);
    exportgraphics(fh, fullfile(figDir, sprintf('fig%02d.png', fn)), 'Resolution', 150);
end
fprintf('Figures saved to %s\n', figDir);

%% ===== EVALUATION =====

fprintf('\n========== EVALUATION — CONTROL: W5 manual RCS + new eval ==========\n');

% --- Week 5 (try103) reference numbers, for delta comparison ---
W5_rangeRes = 7.88;  W5_crossRes = 20.00;
W5_pslr_r   = -3.4;  W5_pslr_cr  = -10.2;
W5_locErrR  = 17.62; W5_locErrC  = 0.00;
W5_snr      = 96.2;

lambda_e    = c / fc;
apertureLen = speed * flightDuration;

fprintf('\n[Radar Parameters]\n');
fprintf('  Center frequency : %.1f GHz\n', fc/1e9);
fprintf('  Bandwidth        : %.1f MHz\n', bw/1e6);
fprintf('  Wavelength       : %.4f m\n', lambda_e);
fprintf('  PRF              : %d Hz  |  Aperture: %.0f m\n', prf, apertureLen);
fprintf('  Platform speed   : %.0f m/s  |  Height: %.0f m\n', speed, radarpos0(3));

fprintf('\n[Box Target]\n');
fprintf('  Center: (range=%.0f, cross=%.0f, height=%.0f) m\n', boxCenter(1), boxCenter(2), boxCenter(3));
fprintf('  Size  : %.0f × %.0f × %.0f m  (L × W × H)\n', boxL, boxW, boxH);
fprintf('  Scatterers: %d total (%d/face)\n', numScatterers, ptsPerFace);

fprintf('\n[RCS — Week 5 manual values (this control run)]\n');
fprintf('  %-12s  %-12s  %-10s\n', 'Face', 'RCS (m^2)', '(dBsm)');
for f = 1:6
    fprintf('  %-12s  %-12.1f  %-10.1f\n', faceLabels{f}, faceRCS(f), 10*log10(faceRCS(f)));
end

fprintf('\n[Theoretical Resolution]\n');
rangeRes_th = c / (2 * bw);
crRes_th    = lambda_e * boxCenter(1) / (2 * apertureLen);
fprintf('  Range      : %.2f m\n', rangeRes_th);
fprintf('  Cross-range: %.4f m (%.2f cm)\n', crRes_th, crRes_th*100);

fprintf('\n[Point-Target Calibration — TRUE system resolution & PSLR]\n');
fprintf('  (from a single scatterer''s impulse response, NOT the box)\n');
fprintf('  Range  res (-3dB): %.2f m   (theory %.2f m)\n', calRangeRes, rangeRes_th);
fprintf('  Cross  res (-3dB): %.3f m   (theory %.4f m)\n', calCrossRes, crRes_th);
fprintf('  PSLR  range      : %.1f dB  (ideal < -13 dB)\n', calPSLR_r);
fprintf('  PSLR  cross      : %.1f dB\n', calPSLR_c);
fprintf('  NOTE: W5 reported res/PSLR were measured on the box -> they were\n');
fprintf('        target EXTENT, not true resolution. These calibration values\n');
fprintf('        are the correct system metrics and stay constant in Step 2.\n');

fprintf('\n[Box Reconstruction — EXTENT of bright region, NOT resolution]\n');
fprintf('  (box is %.0f x %.0f m, so -3/-6 dB widths track target SIZE)\n', boxL, boxW);
[~, xZeroE] = min(abs(xScene));
rProfE   = bpImageMag(:, xZeroE);
rProfEdB = 20*log10(rProfE + eps);  rProfEdB = rProfEdB - max(rProfEdB);
above3rE = yScene(rProfEdB >= -3);
above6rE = yScene(rProfEdB >= -6);
if ~isempty(above3rE)
    fprintf('  Range -3 dB extent: %.2f m  (box length = %.0f m)\n', max(above3rE)-min(above3rE), boxL);
end
if ~isempty(above6rE)
    fprintf('  Range -6 dB extent: %.2f m\n', max(above6rE)-min(above6rE));
end
[~, rCtrE] = min(abs(yScene - boxCenter(1)));
crProfE   = bpImageMag(rCtrE, :);
crProfEdB = 20*log10(crProfE + eps);  crProfEdB = crProfEdB - max(crProfEdB);
above3crE = xScene(crProfEdB >= -3);
if ~isempty(above3crE)
    fprintf('  Cross -3 dB extent: %.2f m  (box width = %.0f m)\n', max(above3crE)-min(above3crE), boxW);
end

fprintf('\n[Dynamic Range (noise-free)]\n');
bpPeakE   = max(bpImageMag(:));
cornerMaskE = (abs(xGrid - boxCenter(2)) > halfW*5) & (abs(yGrid - boxCenter(1)) > halfL*5);
noiseFloorE = mean(bpImageMag(cornerMaskE));
snr_dBE = 20*log10(bpPeakE / noiseFloorE);
fprintf('  Peak / background: %.1f dB  (dynamic range; NOT true SNR -- no noise\n', snr_dBE);
fprintf('                     yet, real SNR comes in Step 2.  W5: 96.2 dB)\n');
fprintf('  Background floor : %.6f  (normalized)\n', noiseFloorE);

fprintf('\n[Box Localization]\n');
brightMaskE = bpImageMag > 0.5;
if sum(brightMaskE(:)) > 0
    xCentE = mean(xGrid(brightMaskE));
    yCentE = mean(yGrid(brightMaskE));
    fprintf('  True center  : (cross=%.0f, range=%.0f) m\n', boxCenter(2), boxCenter(1));
    fprintf('  BP centroid  : (cross=%.1f, range=%.1f) m\n', xCentE, yCentE);
    fprintf('  Location err : cross=%.2f m,  range=%.2f m\n', ...
        abs(xCentE-boxCenter(2)), abs(yCentE-boxCenter(1)));
end

fprintf('\n[Control Summary — W5 RCS measured with the NEW eval]\n');
fprintf('  - Calibration (resolution/PSLR) should MATCH try104 exactly,\n');
fprintf('    because it uses a point target RCS=1, independent of box RCS.\n');
fprintf('  - Box extent/localization/dyn-range should reproduce the OLD W5\n');
fprintf('    numbers (range ext ~7.88 m, loc err ~17.62 m, ~96.2 dB).\n');
fprintf('  => Confirms the W5 "resolution" values were target EXTENT, and that\n');
fprintf('     the W5 vs W6 difference is purely RCS weighting, not the system.\n');

%% ===== FIGURE 11: Evaluation Text (Command Window output) =====

% Build text string from computed variables
w5vals_fig = [2.0, 0.5, 1.0, 1.0, 0.8, 0.3];

evalLines = {};
evalLines{end+1} = '====== CONTROL — W5 manual RCS + Week 6 evaluation ======';
evalLines{end+1} = '';
evalLines{end+1} = '[Radar Parameters]';
evalLines{end+1} = sprintf('  Center freq : %.1f GHz   BW: %.1f MHz', fc/1e9, bw/1e6);
evalLines{end+1} = sprintf('  Wavelength  : %.4f m      PRF: %d Hz', lambda_e, prf);
evalLines{end+1} = sprintf('  Speed: %.0f m/s   Height: %.0f m   Aperture: %.0f m', speed, radarpos0(3), apertureLen);
evalLines{end+1} = '';
evalLines{end+1} = '[Box Target]';
evalLines{end+1} = sprintf('  Center: range=%.0f m, cross=%.0f m', boxCenter(1), boxCenter(2));
evalLines{end+1} = sprintf('  Size  : %.0f x %.0f x %.0f m   Scatterers: %d', boxL, boxW, boxH, numScatterers);
evalLines{end+1} = '';
evalLines{end+1} = '[RCS — Week 5 manual values]';
evalLines{end+1} = sprintf('  %-8s  %-12s  %s', 'Face','RCS (m^2)','(dBsm)');
evalLines{end+1} = repmat('-', 1, 40);
for f = 1:6
    evalLines{end+1} = sprintf('  %-8s  %-12.1f  %.1f', ...
        faceLabels{f}, faceRCS(f), 10*log10(faceRCS(f)));
end
evalLines{end+1} = '';
evalLines{end+1} = '[Theoretical Resolution]';
evalLines{end+1} = sprintf('  Range      : %.2f m  (c/2B)', rangeRes_th);
evalLines{end+1} = sprintf('  Cross-range: %.4f m = %.2f cm  (lambdaR/2L)', crRes_th, crRes_th*100);
evalLines{end+1} = '';
evalLines{end+1} = '[Calibration — TRUE resolution (single point, not box)]';
evalLines{end+1} = sprintf('  Range res : %.2f m   (theory %.2f m)', calRangeRes, rangeRes_th);
evalLines{end+1} = sprintf('  Cross res : %.3f m   (theory %.4f m)', calCrossRes, crRes_th);
evalLines{end+1} = sprintf('  PSLR range: %.1f dB   PSLR cross: %.1f dB', calPSLR_r, calPSLR_c);
evalLines{end+1} = '  (W5 res/PSLR were measured on the box = target EXTENT, invalid)';
evalLines{end+1} = '';
evalLines{end+1} = '[Box Reconstruction — bright-region EXTENT, not resolution]';
if ~isempty(above3rE)
    evalLines{end+1} = sprintf('  Range -3dB extent : %.2f m   (box L = %.0f m)', ...
        max(above3rE)-min(above3rE), boxL);
end
if ~isempty(above3crE)
    evalLines{end+1} = sprintf('  Cross -3dB extent : %.2f m   (box W = %.0f m)', ...
        max(above3crE)-min(above3crE), boxW);
end
evalLines{end+1} = '';
evalLines{end+1} = '[Dynamic Range (noise-free)]';
evalLines{end+1} = sprintf('  Peak/background : %.1f dB   (NOT true SNR; W5: 96.2 dB)', snr_dBE);
evalLines{end+1} = sprintf('  Background floor: %.6f  (normalized)', noiseFloorE);
evalLines{end+1} = '';
evalLines{end+1} = '[Box Localization]';
if exist('xCentE','var')
    evalLines{end+1} = sprintf('  True center : cross=%.0f m, range=%.0f m', boxCenter(2), boxCenter(1));
    evalLines{end+1} = sprintf('  BP centroid : cross=%.1f m, range=%.1f m', xCentE, yCentE);
    evalLines{end+1} = sprintf('  Location err: cross=%.2f m,  range=%.2f m  [W5: 17.62 m]', ...
        abs(xCentE-boxCenter(2)), abs(yCentE-boxCenter(1)));
end
evalLines{end+1} = '';
evalLines{end+1} = '[Control Summary — W5 RCS, new eval]';
evalLines{end+1} = sprintf('  Calibration res : range %.2f m, cross %.3f m', calRangeRes, calCrossRes);
evalLines{end+1} = '    -> should equal try104 (point target, RCS-independent)';
if ~isempty(above3rE)
    mR = max(above3rE)-min(above3rE);
    evalLines{end+1} = sprintf('  Box range extent: %.2f m   (old W5: %.2f m)', mR, W5_rangeRes);
end
if exist('snr_dBE','var')
    evalLines{end+1} = sprintf('  Peak/background : %.1f dB   (old W5: %.1f dB)', snr_dBE, W5_snr);
end
if exist('yCentE','var')
    leR = abs(yCentE-boxCenter(1));
    evalLines{end+1} = sprintf('  Loc err (range) : %.2f m   (old W5: %.2f m)', leR, W5_locErrR);
end
evalLines{end+1} = '  If box metrics match old W5 but calibration matches try104,';
evalLines{end+1} = '  it proves: W5 "resolution" was target EXTENT; true resolution';
evalLines{end+1} = '  is set by BW/aperture and is identical for W5 and W6 RCS.';
evalLines{end+1} = repmat('=', 1, 54);

% Display as Figure 11
figure(11);
set(gcf, 'Name', 'Fig11: Evaluation Text', 'Color', 'k', ...
    'Position', [100 60 770 860]);
ax11 = axes('Position', [0 0 1 1], 'Visible', 'off', ...
    'Color', 'k', 'XColor', 'k', 'YColor', 'k');
hold(ax11, 'on');

nLines = numel(evalLines);
yStep  = 1 / (nLines + 2);
for li = 1:nLines
    line = evalLines{li};
    % Color coding
    if startsWith(line, '===') || startsWith(line, '---')
        clr = [0.6 0.6 0.6];   % grey separator
        fw  = 'normal';
    elseif startsWith(line, '[')
        clr = [0.3 0.85 1.0];  % cyan header
        fw  = 'bold';
    elseif contains(line, 'W5:') || contains(line, '[W5')
        clr = [1.0 0.75 0.2];  % orange for W5 comparison
        fw  = 'normal';
    elseif startsWith(line, '  +')
        clr = [0.3 1.0 0.4];   % green for improvement
        fw  = 'normal';
    elseif startsWith(line, '  -')
        clr = [1.0 0.4 0.4];   % red for worse
        fw  = 'normal';
    elseif startsWith(line, '  =')
        clr = [0.8 0.8 0.8];   % white for neutral
        fw  = 'normal';
    else
        clr = [0.92 0.92 0.92];
        fw  = 'normal';
    end
    text(0.02, 1 - li*yStep, line, ...
        'Units', 'normalized', 'Color', clr, ...
        'FontName', 'Courier', 'FontSize', 8, ...
        'FontWeight', fw, 'Interpreter', 'none', ...
        'Parent', ax11);
end
hold(ax11, 'off');

% export Fig 11 now that it exists (save loop above only covered 1..10)
exportgraphics(figure(11), fullfile(figDir, 'fig11.png'), 'Resolution', 150);

fprintf('\nDone. Figures 1–12 saved (11 = eval text, 12 = calibration IRF).\n');
fprintf('=============================================================\n');

%% ===== LOCAL FUNCTIONS =====

function [Verts, Faces] = makeBoxPatch(cx, cy, cz, hl, hw, hh)
    Verts = [cx-hl, cy-hw, cz;
             cx+hl, cy-hw, cz;
             cx+hl, cy+hw, cz;
             cx-hl, cy+hw, cz;
             cx-hl, cy-hw, cz+hh;
             cx+hl, cy-hw, cz+hh;
             cx+hl, cy+hw, cz+hh;
             cx-hl, cy+hw, cz+hh];
    Faces = [1 2 3 4;
             5 6 7 8;
             1 2 6 5;
             3 4 8 7;
             1 4 8 5;
             2 3 7 6];
end

function [res3dB, pslr] = irfMetrics(ax, pdB)
% irfMetrics — measure -3 dB width and PSLR from a normalized impulse response.
%   ax   : axis values (column vector)
%   pdB  : profile in dB, peak normalized to 0 dB (column vector)
%   res3dB : -3 dB main-lobe width (interpolated)
%   pslr   : peak side-lobe ratio (dB), max level outside the main lobe
    ax  = ax(:);  pdB = pdB(:);
    [~, ip] = max(pdB);                       % main-lobe peak (0 dB)
    xl = interpCross(ax, pdB, ip, -3, -1);    % left  -3 dB crossing
    xr = interpCross(ax, pdB, ip, -3, +1);    % right -3 dB crossing
    res3dB = xr - xl;
    hw  = max(ax(ip) - xl, xr - ax(ip));      % main-lobe half width
    out = abs(ax - ax(ip)) > 2*hw;            % region outside the main lobe
    if any(out)
        pslr = max(pdB(out));                 % highest side lobe (dB, < 0)
    else
        pslr = NaN;
    end
end

function xc = interpCross(ax, pdB, ip, lvl, dir)
% walk from the peak (index ip) in direction dir (+1/-1) to the first sample
% below lvl, then linearly interpolate the exact crossing on the axis.
    n = numel(pdB);
    i = ip;
    while (i+dir) >= 1 && (i+dir) <= n && pdB(i+dir) > lvl
        i = i + dir;
    end
    j = i + dir;
    if j < 1 || j > n
        xc = ax(i);                           % no crossing in window -> clamp
    else
        t  = (lvl - pdB(i)) / (pdB(j) - pdB(i));
        xc = ax(i) + t*(ax(j) - ax(i));
    end
end
