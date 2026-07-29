%% SAR Backprojection — 3D Box Target Simulation
%  หมายเหตุ: บล็อก RADAR CONFIGURATION สืบทอดมาจาก 01_point_target/point3_bp_2d.m
%  ซึ่งดัดแปลงจากตัวอย่าง Stripmap SAR ของ MathWorks (ดู attribution ในไฟล์นั้น)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  ต่อยอดจาก try101 โดยแทน 3 point targets ด้วยกล่องสี่เหลี่ยม 3D
%  Box ถูก discretize เป็น Point Scatterers บน 6 หน้า (faces)
%  เป้าหมาย: ดูว่า SAR image ของ 3D shape หน้าตาเป็นยังไง
%             เปรียบเทียบกับ Point Target
clear; clc; close all;

%% ===== RADAR CONFIGURATION (เหมือนเดิมทุกอย่าง) =====

c  = physconst('LightSpeed');
fc = 4e9;               % C-band carrier frequency
rangeResolution = 3;    % m
bw = c / (2*rangeResolution);  % LFM bandwidth

prf   = 1000;           % Pulse Repetition Frequency (Hz)
aperture = 4;           % antenna aperture (m)
tpd   = 3e-6;           % pulse duration (s)
fs    = 120e6;          % ADC sampling frequency (Hz)

waveform = phased.LinearFMWaveform( ...
    'SampleRate',    fs, ...
    'PulseWidth',    tpd, ...
    'PRF',           prf, ...
    'SweepBandwidth', bw);

% Platform
speed         = 100;        % radar speed (m/s) along y-axis
flightDuration = 4;         % seconds

radarPlatform = phased.Platform( ...
    'InitialPosition', [0; -200; 500], ...
    'Velocity',        [0; speed; 0]);

slowTime   = 1/prf;
numpulses  = flightDuration/slowTime + 1;

% Range sampling
maxRange          = 2500;
truncrangesamples = ceil((2*maxRange/c)*fs);
fastTime          = (0:1/fs:(truncrangesamples-1)/fs);

% Antenna / RF chain
antenna     = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);
antennaGain = aperture2gain(aperture, c/fc);

transmitter = phased.Transmitter('PeakPower', 50e3, 'Gain', antennaGain);

radiator = phased.Radiator('Sensor', antenna, ...
    'OperatingFrequency', fc, ...
    'PropagationSpeed',   c);

collector = phased.Collector('Sensor', antenna, ...
    'PropagationSpeed',   c, ...
    'OperatingFrequency', fc);

receiver = phased.ReceiverPreamp('SampleRate', fs, 'NoiseFigure', 30);

channel = phased.FreeSpace( ...
    'PropagationSpeed',  c, ...
    'OperatingFrequency', fc, ...
    'SampleRate',        fs, ...
    'TwoWayPropagation', true);

%% ===== 3D BOX DEFINITION =====
%
%  กล่องสี่เหลี่ยม กว้าง x ยาว x สูง = boxW x boxL x boxH (เมตร)
%  วางอยู่บนพื้น (z = 0) ตรงกลางที่ตำแหน่ง boxCenter
%
%  แต่ละหน้าถูก discretize ด้วย N_face x N_face จุด
%  RCS ของแต่ละ scatterer ถูกปรับตาม face orientation
%  (หน้าที่หันหา radar มี RCS สูงกว่า)

boxCenter = [1000; 0; 0];   % [range; cross-range; height]  (m)
boxL = 30;   % ความยาวตามแนว range       (m)
boxW = 20;   % ความกว้างตามแนว cross-range (m)
boxH = 10;   % ความสูง                    (m)

N_face = 6;  % จำนวน scatterer ต่อด้านในแต่ละมิติ (N x N ต่อหน้า)

% ─── สร้าง scatterer grid บน 6 หน้า ───────────────────────────────────
% แต่ละหน้าเก็บเป็น 3×Npts (x=range, y=cross-range, z=height)

scatPos = [];   % 3 × totalScatterers
scatRCS = [];   % 1 × totalScatterers  (per-scatterer RCS ในหน่วย m²)

halfL = boxL/2;
halfW = boxW/2;
halfH = boxH/2;

% --- กำหนด grid ของแต่ละหน้า ---
u = linspace(-1, 1, N_face);   % normalized grid -1 to 1

% Face 1: หน้าด้านใกล้ radar (range = boxCenter(1) - halfL) — FRONT
%         หันหา radar มาก → RCS สูง
[V, W] = meshgrid(u*halfW, u*halfH);
xF = (boxCenter(1) - halfL) * ones(size(V));
yF = boxCenter(2) + V;
zF = boxCenter(3) + halfH + W;  % z from 0 to boxH
scatPos = [scatPos, [xF(:)'; yF(:)'; zF(:)']];
scatRCS = [scatRCS, 2.0 * ones(1, N_face^2)];

% Face 2: หน้าด้านไกล radar (range = boxCenter(1) + halfL) — BACK
%         หันหนีจาก radar → RCS ต่ำ
xF = (boxCenter(1) + halfL) * ones(size(V));
yF = boxCenter(2) + V;
zF = boxCenter(3) + halfH + W;
scatPos = [scatPos, [xF(:)'; yF(:)'; zF(:)']];
scatRCS = [scatRCS, 0.5 * ones(1, N_face^2)];

% Face 3: หน้าด้านซ้าย (cross-range = boxCenter(2) - halfW) — LEFT
[U, W] = meshgrid(u*halfL, u*halfH);
xF = boxCenter(1) + U;
yF = (boxCenter(2) - halfW) * ones(size(U));
zF = boxCenter(3) + halfH + W;
scatPos = [scatPos, [xF(:)'; yF(:)'; zF(:)']];
scatRCS = [scatRCS, 1.0 * ones(1, N_face^2)];

% Face 4: หน้าด้านขวา (cross-range = boxCenter(2) + halfW) — RIGHT
xF = boxCenter(1) + U;
yF = (boxCenter(2) + halfW) * ones(size(U));
zF = boxCenter(3) + halfH + W;
scatPos = [scatPos, [xF(:)'; yF(:)'; zF(:)']];
scatRCS = [scatRCS, 1.0 * ones(1, N_face^2)];

% Face 5: หน้าบน (z = boxCenter(3) + boxH) — TOP
[U, V] = meshgrid(u*halfL, u*halfW);
xF = boxCenter(1) + U;
yF = boxCenter(2) + V;
zF = (boxCenter(3) + boxH) * ones(size(U));
scatPos = [scatPos, [xF(:)'; yF(:)'; zF(:)']];
scatRCS = [scatRCS, 0.8 * ones(1, N_face^2)];

% Face 6: พื้นกล่อง (z = 0) — BOTTOM (usually shadowed)
xF = boxCenter(1) + U;
yF = boxCenter(2) + V;
zF = zeros(size(U));
scatPos = [scatPos, [xF(:)'; yF(:)'; zF(:)']];
scatRCS = [scatRCS, 0.3 * ones(1, N_face^2)];

numScatterers = size(scatPos, 2);
fprintf('Box discretized into %d point scatterers (%d faces × %d² pts)\n', ...
    numScatterers, 6, N_face);

%% ===== FIGURE 1: Ground Truth — Box Geometry =====

figure(1);
% วาด outline ของกล่อง (top-down view: cross-range vs range)
boxX = boxCenter(1) + halfL * [-1 1 1 -1 -1];   % range corners
boxY = boxCenter(2) + halfW * [-1 -1 1 1 -1];    % cross-range corners
plot(boxY, boxX, 'k-', 'LineWidth', 2);
hold on;

% วาด scatterer positions แยกตามหน้า
faceColors = {'r','b','g','m','c','y'};
faceLabels = {'Front','Back','Left','Right','Top','Bottom'};
ptsPerFace = N_face^2;
for f = 1:6
    idx = (f-1)*ptsPerFace + (1:ptsPerFace);
    scatter(scatPos(2,idx), scatPos(1,idx), 20, faceColors{f}, 'filled', ...
        'DisplayName', faceLabels{f});
end

legend('Box outline', faceLabels{:}, 'Location', 'best');
xlabel('Cross-Range (m)');
ylabel('Range (m)');
title(sprintf('Ground Truth: 3D Box Target (%dm × %dm × %dm)', boxL, boxW, boxH));
grid on;
axis equal;
hold off;

%% ===== SIMULATE PLATFORM POSITION HISTORY =====

radarpos0 = [0; -200; 500];
radarvel0 = [0; speed; 0];
pulseIndex = 1:numpulses;
timeAxis   = (pulseIndex - 1) / prf;
radarPosHistory = radarpos0 + radarvel0 * timeAxis;  % 3 × numpulses

%% ===== SIMULATE RAW ECHO SIGNAL =====
%
%  ต่างจากเดิม: target ไม่ใช่ phased.Platform เดียว
%  แต่เป็น scatterers หลายร้อยจุด → ต้อง loop ทีละ scatterer แล้ว sum
%  เพื่อใช้ phased toolbox ได้ เราส่ง scatterers ทีละ batch

fprintf('Simulating raw echo for Box target (%d scatterers)...\n', numScatterers);

rxsig = zeros(truncrangesamples, numpulses);
refangle = zeros(1, numScatterers);

% สร้าง Platform และ RadarTarget สำหรับ scatterers ทั้งหมดพร้อมกัน
scatVel = zeros(size(scatPos));   % scatterers หยุดนิ่ง
boxPlatform  = phased.Platform('InitialPosition', scatPos, 'Velocity', scatVel);
boxTarget    = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', scatRCS);

% Reset radar platform
radarPlatform2 = phased.Platform( ...
    'InitialPosition', [0; -200; 500], ...
    'Velocity',        [0; speed; 0]);

for ii = 1:numpulses
    [radarpos, radarvel] = radarPlatform2(slowTime);
    [tpos, tvel]         = boxPlatform(slowTime);

    [targetRange, targetAngle] = rangeangle(tpos, radarpos);

    sig = waveform();
    sig = sig(1:truncrangesamples);
    sig = transmitter(sig);

    % บังคับ azimuth = 0 (broadside beam ไม่ tilt)
    targetAngle(1,:) = zeros(1, numScatterers);

    sig = radiator(sig, targetAngle);
    sig = channel(sig, radarpos, tpos, radarvel, tvel);
    sig = boxTarget(sig);
    sig = collector(sig, targetAngle);
    rxsig(:,ii) = receiver(sig);
end

fprintf('Raw echo simulation complete.\n');

%% ===== FIGURE 2 & 3: Raw Echo Data =====

figure(2);
imagesc(abs(rxsig));
title('SAR Raw Echo — Box Target (Magnitude)');
xlabel('Slow-time / Pulse index');
ylabel('Fast-time / Range samples');
colorbar;

rangeAxis = fastTime * c / 2;

figure(3);
imagesc(1:numpulses, rangeAxis, abs(rxsig));
set(gca, 'YDir', 'normal');
title('SAR Raw Echo with Range Axis — Box Target');
xlabel('Slow-time / Pulse index');
ylabel('Slant range (m)');
colorbar;
ylim([800 1300]);

%% ===== RANGE COMPRESSION =====

refPulse    = waveform();
chirpSamples = round(tpd * fs);
refChirp    = refPulse(1:chirpSamples);
matchedFilter = conj(flipud(refChirp));

rxsigRC = zeros(size(rxsig,1) + length(matchedFilter) - 1, numpulses);

for ii = 1:numpulses
    rxsigRC(:,ii) = conv(rxsig(:,ii), matchedFilter, 'full');
end

rangeAxisRC = ((0:size(rxsigRC,1)-1).' - (length(matchedFilter)-1)) / fs * c/2;

rxsigRCmag = abs(rxsigRC);
rxsigRCmag = rxsigRCmag ./ max(rxsigRCmag(:));

figure(4);
imagesc(1:numpulses, rangeAxisRC, rxsigRCmag);
set(gca, 'YDir', 'normal');
title('Range Compressed SAR Data — Box Target');
xlabel('Slow-time / Pulse index');
ylabel('Slant range (m)');
colorbar;
ylim([850 1200]);

%% ===== BACKPROJECTION =====

fprintf('Running Backprojection for Box target...\n');

% Scene grid — ขยายให้ครอบ box ทั้งหมด
xScene = linspace(-100, 100, 401);   % cross-range (m)
yScene = linspace(850, 1200, 401);   % ground range (m)

[xGrid, yGrid] = meshgrid(xScene, yScene);
Ny = length(yScene);
Nx = length(xScene);

% 3D pixel positions (height=0 สำหรับ ground plane)
pixelPos = [yGrid(:)'; xGrid(:)'; zeros(1, numel(xGrid))];

bpImage = zeros(Ny, Nx);

for ii = 1:numpulses
    rp = radarPosHistory(:, ii);

    diffVec    = pixelPos - rp;
    slantRange = sqrt(sum(diffVec.^2, 1));

    sampleIdx = 2 .* slantRange ./ c .* fs + length(matchedFilter);
    validMask = (sampleIdx >= 1) & (sampleIdx <= size(rxsigRC,1));

    interpVals = zeros(1, numel(xGrid));
    if any(validMask)
        interpVals(validMask) = interp1( ...
            1:size(rxsigRC,1), ...
            rxsigRC(:,ii), ...
            sampleIdx(validMask), ...
            'linear', 0);
    end

    phaseCorr    = exp(1j * 4*pi*fc/c .* slantRange);
    contribution = reshape(interpVals .* phaseCorr, Ny, Nx);
    bpImage      = bpImage + contribution;
end

fprintf('Backprojection complete.\n');

bpImageMag = abs(bpImage);
bpImageMag = bpImageMag ./ max(bpImageMag(:));
bpImagedB  = 20 * log10(bpImageMag + eps);

%% ===== FIGURE 5: BP Image (Linear Scale) =====

figure(5);
imagesc(xScene, yScene, bpImageMag);
set(gca, 'YDir', 'normal');
xlabel('Cross-range (m)');
ylabel('Ground range (m)');
title('SAR BP Image — Box Target (Linear Scale)');
colorbar;
axis equal tight;

% วาด box outline ทับ
hold on;
boxX_cr = boxCenter(2) + halfW * [-1 1 1 -1 -1];  % cross-range corners
boxY_r  = boxCenter(1) + halfL * [-1 -1 1 1 -1];  % range corners
plot(boxX_cr, boxY_r, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Box outline (truth)');
legend('Location','northeast');
hold off;

%% ===== FIGURE 6: BP Image (dB Scale) =====

figure(6);
imagesc(xScene, yScene, bpImagedB, [-40 0]);
set(gca, 'YDir', 'normal');
colormap('jet');
xlabel('Cross-range (m)');
ylabel('Ground range (m)');
title('SAR BP Image — Box Target (dB Scale, clipped −40 dB)');
colorbar;
axis equal tight;

hold on;
plot(boxX_cr, boxY_r, 'w--', 'LineWidth', 1.5);
legend('Box outline (truth)', 'Location','northeast');
hold off;

%% ===== FIGURE 7: Range Profile (เปรียบเทียบ Front / Center / Back) =====
% ตัด cross-section ที่ x=0 เพื่อดู range response ของ box

[~, xZeroIdx] = min(abs(xScene));
rangeProfile   = bpImageMag(:, xZeroIdx);

figure(7);
plot(yScene, rangeProfile, 'b-', 'LineWidth', 1.5);
hold on;
xline(boxCenter(1) - halfL, '--r', 'Front face',  'LabelVerticalAlignment','bottom');
xline(boxCenter(1),          '--k', 'Box center',  'LabelVerticalAlignment','bottom');
xline(boxCenter(1) + halfL, '--g', 'Back face',   'LabelVerticalAlignment','bottom');
hold off;
grid on;
xlabel('Ground range (m)');
ylabel('Normalized amplitude');
title('BP Range Profile at Cross-Range = 0 m');
xlim([boxCenter(1)-halfL-20, boxCenter(1)+halfL+20]);

%% ===== FIGURE 8: Cross-Range Profile (Front / Center / Back faces) =====

figure(8);
hold on;
rangeSlices = [boxCenter(1)-halfL, boxCenter(1), boxCenter(1)+halfL];
sliceLabels = {'Front face (near)','Box center','Back face (far)'};
sliceColors = {'r','k','g'};

for k = 1:3
    [~, rIdx] = min(abs(yScene - rangeSlices(k)));
    plot(xScene, bpImageMag(rIdx,:), sliceColors{k}, ...
        'LineWidth', 1.5, 'DisplayName', sliceLabels{k});
end
hold off;
grid on;
xlabel('Cross-range (m)');
ylabel('Normalized amplitude');
title('BP Cross-Range Profile at Front / Center / Back of Box');
legend show;
xlim([-60 60]);

%% ===== FIGURE 9: 3D Visualization ของ Box Scatterer Geometry =====

figure(9);
scatter3(scatPos(2,:), scatPos(1,:), scatPos(3,:), ...
    30, scatRCS, 'filled');
colorbar;
xlabel('Cross-range (m)');
ylabel('Range (m)');
zlabel('Height (m)');
title(sprintf('3D Box Scatterer Positions (colored by RCS)  [%d pts]', numScatterers));
grid on;
view(35, 25);

%% ===== EVALUATION =====

fprintf('\n========== EVALUATION — Box Target ==========\n');
fprintf('Box size: %.0fm (range) x %.0fm (cross-range) x %.0fm (height)\n', ...
    boxL, boxW, boxH);
fprintf('Box center: range=%.0fm, cross-range=%.0fm\n', ...
    boxCenter(1), boxCenter(2));
fprintf('Total scatterers: %d\n', numScatterers);

% Range extent of bright region in BP image
rangeProfileDB = 20*log10(rangeProfile + eps);
rangeProfileDB = rangeProfileDB - max(rangeProfileDB);
above6dB = yScene(rangeProfileDB >= -6);  % -6 dB extent
if ~isempty(above6dB)
    fprintf('\n--- BP Image Range Extent (−6 dB) ---\n');
    fprintf('  From %.1fm to %.1fm (width = %.1fm)\n', ...
        min(above6dB), max(above6dB), max(above6dB)-min(above6dB));
    fprintf('  Box true length = %.0fm\n', boxL);
end

% Cross-range extent at front face
[~, rFront] = min(abs(yScene - (boxCenter(1)-halfL)));
crProfileFront = bpImageMag(rFront,:);
crProfileFrontDB = 20*log10(crProfileFront + eps);
crProfileFrontDB = crProfileFrontDB - max(crProfileFrontDB);
above6dB_cr = xScene(crProfileFrontDB >= -6);
if ~isempty(above6dB_cr)
    fprintf('\n--- BP Image Cross-Range Extent at Front Face (−6 dB) ---\n');
    fprintf('  From %.1fm to %.1fm (width = %.1fm)\n', ...
        min(above6dB_cr), max(above6dB_cr), max(above6dB_cr)-min(above6dB_cr));
    fprintf('  Box true width = %.0fm\n', boxW);
end

fprintf('\n==============================================\n');
fprintf('Done. Check Figures 1–9.\n');