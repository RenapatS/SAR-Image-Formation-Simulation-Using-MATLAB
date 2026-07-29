%% Configuration
%  -------------------------------------------------------------------
%  ATTRIBUTION / ที่มาของโค้ด
%  -------------------------------------------------------------------
%  ส่วน RADAR CONFIGURATION, SCENE CONFIGURATION และลูปสร้าง RAW ECHO
%  (ตั้งแต่ต้นไฟล์ถึงท้ายลูป "Simulate Raw Echo Signal")
%  ดัดแปลงมาจากตัวอย่างของ MathWorks:
%
%     "Stripmap Synthetic Aperture Radar (SAR) Image Formation"
%     Radar Toolbox Documentation, The MathWorks, Inc.
%     https://www.mathworks.com/help/radar/ug/stripmap-synthetic-aperture-radar-sar-image-formation.html
%     Portions Copyright The MathWorks, Inc.
%
%  ส่วนที่เขียนขึ้นเองในไฟล์นี้ (ไม่ได้มาจากตัวอย่างข้างต้น):
%     - matched filter สำหรับ range compression (conj(flipud(chirp)) + conv)
%       ตัวอย่างของ MathWorks ใช้ phased.RangeResponse แทน
%     - backprojection แบบ exact ที่คำนวณ slant range ต่อพิกเซลแล้วแก้เฟส
%       ตัวอย่างของ MathWorks ใช้ approximate BP + hanning window คนละอัลกอริทึม
%     - Fig 4, 5, 7, 8, 13 และส่วน EVALUATION ทั้งหมด
%
%  Reference เชิงแนวคิดเพิ่มเติม (ไม่ได้คัดลอกโค้ด):
%     Fahnemann C., Rother N., Blume H., "Interactive synthetic aperture radar
%     simulator generating and visualizing realistic FMCW data,"
%     Int. Conf. on Radar Systems (RADAR 2022), pp. 725-730. doi:10.1049/ICP.2023.1281
%  -------------------------------------------------------------------
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
clear; clc; close all;

% ---- โฟลเดอร์เก็บรูป (สร้างให้เองถ้ายังไม่มี) ----
figDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figure');
if ~exist(figDir,'dir'), mkdir(figDir); end
% helper: เซฟ figure ปัจจุบันเป็น png
saveFig = @(name) exportgraphics(gcf, fullfile(figDir,name), 'Resolution', 150);

% Radar Config

c = physconst('LightSpeed');
fc = 4e9; % C-band

rangeResolution = 3;
crossRangeResolution = 3;

bw = c/(2*rangeResolution); % LFM sweep bandwidth (Hz)

prf = 1000; % Pulse Repetition Frequency
aperture = 4;  % (m)
tpd = 3*10^-6;  % pulse duration(s)
fs = 120*10^6;  % sampling frequency

waveform = phased.LinearFMWaveform('SampleRate',fs, ...
    'PulseWidth', tpd, ...
    'PRF', prf, ...
    'SweepBandwidth', bw);

% Platform Config

speed = 100;  
flightDuration = 4; % ระยะเวลาที่ radar เก็บข้อมูล

radarPlatform  = phased.Platform( ...
    'InitialPosition', [0;-200;500], ...
    'Velocity', [0; speed; 0]);

slowTime = 1/prf;
numpulses = flightDuration/slowTime +1;

% Range Sampling Config

maxRange = 2500;
truncrangesamples = ceil((2*maxRange/c)*fs);
fastTime = (0:1/fs:(truncrangesamples-1)/fs);
Rc = 1000; % reference range

% antenna
antenna = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);
antennaGain = aperture2gain(aperture,c/fc);

%เพิ่ม power/gain ให้ signal
transmitter = phased.Transmitter('PeakPower', 50e3, 'Gain', antennaGain); 
%ยิง signal ออกจาก antenna ไปตามมุมที่กำหนด
radiator = phased.Radiator('Sensor', antenna,...
    'OperatingFrequency', fc, ...
    'PropagationSpeed', c);
    
% antenna รับ echo ที่กลับมาจากทิศของ target
collector = phased.Collector('Sensor', antenna, ...
    'PropagationSpeed', c,...
    'OperatingFrequency', fc);
% รับ echo ที่กลับมา แล้วเก็บลง raw data matrix
receiver = phased.ReceiverPreamp('SampleRate', fs, 'NoiseFigure', 30);

%จำลองการเดินทางของคลื่นในอากาศ/อวกาศว่าง ระหว่าง radar กับ target
channel = phased.FreeSpace('PropagationSpeed', c, ...
    'OperatingFrequency', fc,...
    'SampleRate', fs,...
    'TwoWayPropagation', true);



%[[next]]

targetpos= [
    800,0,0;
    1000,0,0; 
    1300,0,0]'; 

% transpose = //index start at 1,1
% 800    1000    1300
% 0       0       0
% 0       0       0

targetvel = [ % velocity = 0 m/s
    0,0,0;
    0,0,0; 
    0,0,0]';

target = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', [1,1,1]); %RCS = radar cross section(ไว้ใช้ตอนสะท้อน)
pointTargets = phased.Platform('InitialPosition', targetpos,'Velocity',targetvel);

% the ground truth based on the target locations.
figure(1);
h = axes;
plot(targetpos(2,1),targetpos(1,1),'*g');
hold on;
plot(targetpos(2,2),targetpos(1,2),'*r');
hold on;plot(targetpos(2,3),targetpos(1,3),'*b');
hold off;

set(h,'Ydir','reverse');
xlim([-10 10]);
ylim([700 1500]);

figure(1);
title('Ground Truth');
ylabel('Range');
xlabel('Cross-Range');
saveFig('fig_w13_01_ground_truth.png');
saveFig('fig_pipe2_scene_setup.png');        % stage 2 ของ pipeline (ใช้ในสไลด์)

%% Simulate Raw Echo Signal
refangle = zeros(1,size(targetpos,2));
rxsig = zeros(truncrangesamples,numpulses);

for ii = 1:numpulses
    % Update radar platform and target position
    [radarpos, radarvel] = radarPlatform(slowTime);
    [targetpos,targetvel] = pointTargets(slowTime);

    % Get the range and angle to the point targets
    [targetRange, targetAngle] = rangeangle(targetpos, radarpos);

    % Generate the LFM pulse (chirp)
    sig = waveform();
    % Use only the pulse length that will cover the targets.
    sig = sig(1:truncrangesamples);

    % Transmit the pulse
    sig = transmitter(sig);

    % Define no tilting of beam in azimuth direction
    targetAngle(1,:) = refangle;

    % Radiate the pulse towards the targets
    sig = radiator(sig, targetAngle);

    % Propagate the pulse to the point targets in free space
    sig = channel(sig, radarpos, targetpos, radarvel, targetvel);

    % Reflect the pulse off the targets
    sig = target(sig);

    % Collect the reflected pulses at the antenna
    sig = collector(sig, targetAngle);

    % Receive the signal  
    rxsig(:,ii) = receiver(sig);

end

%% Visualize Raw Echo Data
figure(2);
imagesc(real(rxsig));
title('SAR Raw Data')
xlabel('Cross-Range Samples')
ylabel('Range Samples')
saveFig('fig_w13_02_rawdata_real.png');

figure(3);
imagesc(abs(rxsig));
title('SAR Raw Echo Data - Magnitude');
xlabel('Slow-time / Pulse index / Cross-range samples');
ylabel('Fast-time / Range samples');
colorbar;
saveFig('fig_w13_03_rawdata_magnitude.png');

rangeAxis = fastTime * c / 2;

figure(4);
imagesc(1:numpulses, rangeAxis, abs(rxsig));
set(gca, 'YDir', 'normal');
title('SAR Raw Echo Data with Range Axis');
xlabel('Slow-time / Pulse index');
ylabel('Slant range (m)');
colorbar;
saveFig('fig_w13_04_rawdata_range_axis.png');

%% Figure 5: Actual raw echo + theoretical slant range curves

% Save original geometry manually
targetpos0 = [
    800,0,0;
    1000,0,0; 
    1300,0,0]';

radarpos0 = [0; -200; 500];
radarvel0 = [0; speed; 0];

pulseIndex = 1:numpulses;
timeAxis = (pulseIndex - 1) / prf;

% radar position history, size = 3 x numpulses
radarPosHistory = radarpos0 + radarvel0 * timeAxis;

targetRangeHistory = zeros(size(targetpos0,2), numpulses);

for k = 1:size(targetpos0,2)
    diffPos = radarPosHistory - targetpos0(:,k);
    targetRangeHistory(k,:) = sqrt(sum(diffPos.^2, 1));
end

figure(5);

imagesc(pulseIndex, rangeAxis, abs(rxsig));
set(gca, 'YDir', 'normal');
xlabel('Slow-time / Pulse index');
ylabel('Slant range (m)');
title('Actual Raw Echo Data with Theoretical Slant Range Curves');
colorbar;
hold on;

plot(pulseIndex, targetRangeHistory(1,:), 'w', 'LineWidth', 2);
plot(pulseIndex, targetRangeHistory(2,:), 'r', 'LineWidth', 2);
plot(pulseIndex, targetRangeHistory(3,:), 'g', 'LineWidth', 2);

ylim([800 1800]);

% หมายเหตุ: legend ไม่นับ object ที่มาจาก imagesc จึงใส่เฉพาะ 3 เส้น
legend('Target 1 theoretical range', ...
    'Target 2 theoretical range', ...
    'Target 3 theoretical range');

hold off;
saveFig('fig_w13_05_rawecho_vs_theory.png');
saveFig('fig_pipe3_raw_echo.png');           % stage 3 ของ pipeline (ใช้ในสไลด์)

%% Figure 13: Radar configuration — chirp + geometry  (stage 1 ของ pipeline)
%  ขั้นแรกของ pipeline เดิมไม่มีรูป เพิ่มไว้เพื่อใช้อธิบายในสไลด์
%  ซ้าย  = LFM chirp ความถี่ทันทีไล่ขึ้นเป็นเส้นตรงตลอดพัลส์
%  ขวา   = เรขาคณิตการบิน 3 มิติ (เส้นทางบิน + เป้า 3 จุด)

figure(13); clf; set(gcf,'Color','w','Position',[60 60 1040 400]);

subplot(1,2,1);
tCh = (0:1/fs:tpd-1/fs);
plot(tCh*1e6, (bw/tpd*tCh)/1e6, 'LineWidth', 2.5, 'Color', [0.17 0.43 0.39]);
grid on; box on;
xlabel('Time within pulse (\mus)');
ylabel('Instantaneous frequency (MHz)');
title(sprintf('LFM chirp  |  B = %.1f MHz , T_p = %.0f \\mus', bw/1e6, tpd*1e6));
xlim([0 tpd*1e6]); ylim([0 bw/1e6*1.05]);
text(0.06*tpd*1e6, 0.80*bw/1e6, sprintf('\\delta r = c/2B = %.2f m', c/(2*bw)), ...
    'BackgroundColor',[0.95 0.97 0.98], 'EdgeColor',[0.72 0.78 0.84], 'Margin',5);

subplot(1,2,2);
plot3(radarPosHistory(2,:), radarPosHistory(1,:), radarPosHistory(3,:), ...
      'LineWidth', 3, 'Color', [0.13 0.19 0.35]); hold on;
plot3(targetpos0(2,:), targetpos0(1,:), targetpos0(3,:), 'r*', ...
      'MarkerSize', 11, 'LineWidth', 2);
midIdx = round(numpulses/2);
for k = 1:size(targetpos0,2)
    plot3([radarPosHistory(2,midIdx) targetpos0(2,k)], ...
          [radarPosHistory(1,midIdx) targetpos0(1,k)], ...
          [radarPosHistory(3,midIdx) targetpos0(3,k)], ':', ...
          'Color',[0.45 0.62 0.76], 'LineWidth', 1.2);
end
hold off; grid on; box on;
xlabel('Cross-range (m)'); ylabel('Ground range (m)'); zlabel('Height (m)');
title(sprintf('Geometry  |  h = %d m , v = %d m/s , aperture = %d m', ...
      round(radarPosHistory(3,1)), speed, round(speed*flightDuration)));
legend({'Flight path','Point targets'}, 'Location','northeast');
view(-37, 22);

saveFig('fig_w13_13_radar_config.png');
saveFig('fig_pipe1_radar_config.png');       % stage 1 ของ pipeline (ใช้ในสไลด์)

%% Range Compression using Matched Filtering
% Range compression = correlate received echo with transmitted LFM chirp
% This compresses each long chirp echo into a narrow peak in range direction.

% Recreate one transmitted baseband chirp for matched filter reference
refPulse = waveform();
chirpSamples = round(tpd * fs);
refChirp = refPulse(1:chirpSamples);

% Matched filter impulse response h(t) = conjugate time-reversed transmitted chirp
matchedFilter = conj(flipud(refChirp));

% Apply matched filter along fast-time/range dimension for every slow-time pulse
rxsigRC = zeros(size(rxsig,1) + length(matchedFilter) - 1, size(rxsig,2));

for ii = 1:numpulses
    rxsigRC(:,ii) = conv(rxsig(:,ii), matchedFilter, 'full');
end

% Range axis after full convolution.
% Subtract matched-filter delay so the peak corresponds to physical slant range.
rangeAxisRC = ((0:size(rxsigRC,1)-1).' - (length(matchedFilter)-1)) / fs * c/2;

% Normalize for easier visualization
rxsigRCmag = abs(rxsigRC);
rxsigRCmag = rxsigRCmag ./ max(rxsigRCmag(:));

figure(6);
imagesc(1:numpulses, rangeAxisRC, rxsigRCmag);
set(gca, 'YDir', 'normal');
xlabel('Slow-time / Pulse index');
ylabel('Slant range (m)');
title('Range Compressed SAR Data');
colorbar;
ylim([800 1800]);
saveFig('fig_w13_06_range_compressed.png');

%% Figure 7: Range-compressed data + theoretical slant range curves
figure(7);
imagesc(1:numpulses, rangeAxisRC, rxsigRCmag);
set(gca, 'YDir', 'normal');
xlabel('Slow-time / Pulse index');
ylabel('Slant range (m)');
title('Range Compressed Data with Theoretical Slant Range Curves');
colorbar;
hold on;

plot(pulseIndex, targetRangeHistory(1,:), 'w', 'LineWidth', 2);
plot(pulseIndex, targetRangeHistory(2,:), 'r', 'LineWidth', 2);
plot(pulseIndex, targetRangeHistory(3,:), 'g', 'LineWidth', 2);

ylim([800 1800]);
% หมายเหตุ: legend ไม่นับ object ที่มาจาก imagesc จึงใส่เฉพาะ 3 เส้น
legend('Target 1 theoretical range', ...
    'Target 2 theoretical range', ...
    'Target 3 theoretical range');
hold off;
saveFig('fig_w13_07_rangecompressed_vs_theory.png');

%% Optional: Compare one pulse before and after range compression
midPulse = round(numpulses/2);

figure(8);
plot(rangeAxis, abs(rxsig(:,midPulse)) ./ max(abs(rxsig(:,midPulse))));
hold on;
plot(rangeAxisRC, rxsigRCmag(:,midPulse));
hold off;
grid on;
xlabel('Slant range (m)');
ylabel('Normalized magnitude');
title('One Pulse Before vs After Range Compression');
legend('Raw echo', 'Range compressed echo');
xlim([700 1800]);
saveFig('fig_w13_08_before_after_rc.png');
saveFig('fig_pipe4_range_compression.png');  % stage 4 ของ pipeline (ใช้ในสไลด์)


%% Backprojection (BP) SAR Focusing
%  Algorithm overview
%  For every pixel p = (x_p, y_p) in the scene grid:
%    1. Compute the instantaneous slant range R(n) = ||radarPos(n) - p|| 
%       for each pulse n.
%    2. Convert R(n) to the corresponding fast-time sample index in the
%       range-compressed data.
%    3. Interpolate rxsigRC at that sample index.
%    4. Apply a phase correction to remove the carrier round-trip phase:
%         phase_corr(n) = exp(+j * 4*pi*fc/c * R(n))
%    5. Coherently sum (integrate) the phase-corrected, interpolated
%       samples across all pulses => focused pixel amplitude.
%
%  The result is a 2-D complex SAR image in (cross-range x, range y)
%  coordinates with both range and azimuth compression applied.
% ======================================================================

fprintf('Running Backprojection... this may take a moment.\n');

%% Scene grid definition

xScene   = linspace(-200, 200, 401);   % cross-range axis  (m)
yScene   = linspace(700,  1500, 401);  % ground range axis (m)

[xGrid, yGrid] = meshgrid(xScene, yScene);   % size: (Ny x Nx)

% 3-D scene pixel positions (height assumed = 0 for ground-plane targets)
% pixelPos: 3 x (Ny*Nx)
pixelPos = [yGrid(:)'; xGrid(:)'; zeros(1, numel(xGrid))];

%% Backprojection loop
Ny = length(yScene);
Nx = length(xScene);
bpImage = zeros(Ny, Nx);  % complex focused SAR image

for ii = 1:numpulses
    
    rp = radarPosHistory(:, ii); % Radar position at pulse ii  (3 x 1)

    diffVec    = pixelPos - rp;                   % 3 x Npix
    slantRange = sqrt(sum(diffVec.^2, 1));         % 1 x Npix

    % Convert slant range to sample index in rangeAxisRC
    sampleIdx = 2 .* slantRange ./ c .* fs + length(matchedFilter);  % 1 x Npix
    
    % Clamp indices to valid range
    validMask = (sampleIdx >= 1) & (sampleIdx <= size(rxsigRC, 1));

    % Interpolate range-compressed pulse at fractional sample positions
    interpVals = zeros(1, numel(xGrid));
    if any(validMask)
        interpVals(validMask) = interp1( ...
            1:size(rxsigRC,1), ...          % sample indices (integer grid)
            rxsigRC(:, ii), ...             % complex RC data for this pulse
            sampleIdx(validMask), ...       % query positions
            'linear', 0);                   % zero outside bounds
    end

    % Phase correction: remove two-way carrier phase accumulated over slantRange
    phaseCorr = exp(1j * 4 * pi * fc / c .* slantRange);

    % Coherent accumulation into image grid
    contribution = reshape(interpVals .* phaseCorr, Ny, Nx);
    bpImage = bpImage + contribution;

end

fprintf('Backprojection complete.\n');

% Normalize
bpImageMag = abs(bpImage);
bpImageMag = bpImageMag ./ max(bpImageMag(:));

% dB scale for dynamic range display
bpImagedB = 20 * log10(bpImageMag + eps);

% ------------------------------------------------------------------
% Figure 9: BP image (linear scale)
% ------------------------------------------------------------------
figure(9);
imagesc(xScene, yScene, bpImageMag);
set(gca, 'YDir', 'normal');
xlabel('Cross-range (m)');
ylabel('Ground range (m)');
title('SAR Backprojection Image (Linear Scale)');
colorbar;
axis equal tight;

% Mark ground truth target positions
hold on;
plot(targetpos0(2,:), targetpos0(1,:), 'r*', 'MarkerSize', 12, 'LineWidth', 2);
legend('Ground truth targets');   % imagesc ไม่ถูกนับใน legend
hold off;
saveFig('fig_w13_09_bp_linear.png');

% ------------------------------------------------------------------
% Figure 10: BP image (dB scale, clipped at -40 dB)
% ------------------------------------------------------------------
figure(10);
imagesc(xScene, yScene, bpImagedB, [-40 0]);
set(gca, 'YDir', 'normal');
colormap('jet');
xlabel('Cross-range (m)');
ylabel('Ground range (m)');
title('SAR Backprojection Image (dB Scale, clipped at -40 dB)');
colorbar;
axis equal tight;

hold on;
plot(targetpos0(2,:), targetpos0(1,:), 'w*', 'MarkerSize', 12, 'LineWidth', 2);
legend('Ground truth targets');   % imagesc ไม่ถูกนับใน legend
hold off;
saveFig('fig_w13_10_bp_db.png');
saveFig('fig_pipe5_backprojection.png');     % stage 5 ของ pipeline (ใช้ในสไลด์)

% ------------------------------------------------------------------
% Figure 11: Range profile through each focused target
%   Slice the BP image along the cross-range=0 axis to see
%   range resolution of each reconstructed target.
% ------------------------------------------------------------------
[~, xZeroIdx] = min(abs(xScene));  % index closest to x = 0

figure(11);
plot(yScene, bpImageMag(:, xZeroIdx));
grid on;
xlabel('Ground range (m)');
ylabel('Normalized amplitude');
title('BP Image Range Profile at Cross-Range = 0 m');
hold on;
for k = 1:size(targetpos0, 2)
    xline(targetpos0(1,k), '--r', sprintf('T%d (%.0f m)', k, targetpos0(1,k)));
end
hold off;
xlim([700 1500]);
saveFig('fig_w13_11_range_profile.png');

% ------------------------------------------------------------------
% Figure 12: Cross-range profile at each target's range bin
%   Shows azimuth (along-track) focusing quality per target.
% ------------------------------------------------------------------
figure(12);
hold on;
colors = {'b','r','g'};
for k = 1:size(targetpos0, 2)
    [~, rIdx] = min(abs(yScene - targetpos0(1,k)));
    plot(xScene, bpImageMag(rIdx, :), colors{k}, ...
         'DisplayName', sprintf('Target %d (range = %.0f m)', k, targetpos0(1,k)));
end
hold off;
grid on;
xlabel('Cross-range (m)');
ylabel('Normalized amplitude');
title('BP Image Cross-Range Profile at Each Target Range');
legend show;
xlim([-100 100]);
saveFig('fig_w13_12_crossrange_profile.png');

%% ===== EVALUATION =====
% measure resolution from BP image

fprintf('========== EVALUATION ==========\n');

% ─────────────────────────────────────────────────────────────
% 1. RANGE RESOLUTION  (วัดจาก Fig 11 — range profile)

%    ตัด BP image ตามแนวตั้ง (range direction) ที่ cross-range = 0
%    แล้ววัดความกว้างของ peak ที่ระดับ half-power (3-dB width)
%    หาร 2 เพราะ sum(above) นับทั้งซ้ายและขวาของ peak รวมกัน
% ─────────────────────────────────────────────────────────────

[~, xZeroIdx] = min(abs(xScene));        % index ที่ x ใกล้ 0 ที่สุด
rangeProfile   = bpImageMag(:, xZeroIdx); % 1-D profile ตามแนว range

halfMax_r      = max(rangeProfile) / 2;
above_r        = rangeProfile >= halfMax_r;
measured_range_res = sum(above_r) * (yScene(2)-yScene(1)) / 2;

theory_range_res = c / (2 * bw);

fprintf('\n--- Range Resolution ---\n');
fprintf('  Theoretical : %.2f m\n', theory_range_res);
fprintf('  Measured    : %.2f m\n', measured_range_res);

% ─────────────────────────────────────────────────────────────
% 2. CROSS-RANGE RESOLUTION  (วัดจาก Fig 12 — cross-range profile)
%
%    ตัด BP image ตามแนวนอน (cross-range direction) ที่ range ของแต่ละ target
%    ใช้ spline interpolation ก่อนวัด เพราะ peak แคบกว่า grid step มาก
%    (true resolution ~0.1 m แต่ grid step = 1 m)
% ─────────────────────────────────────────────────────────────

targets_range = [800, 1000, 1300];   % range ของ T1, T2, T3 (m)
L = speed * flightDuration;          % aperture length = 100 × 4 = 400 m
lambda = c / fc;                     % wavelength = 0.075 m

fprintf('\n--- Cross-Range Resolution ---\n');
fprintf('  (interpolated to 0.01 m grid before measuring)\n');

for k = 1:3
    % ดึง cross-range profile ที่ range ของ target k
    [~, rIdx]   = min(abs(yScene - targets_range(k)));
    crProfile   = bpImageMag(rIdx, :);

    % หา peak position จริง (อาจไม่ตรง x=0 พอดีถ้า grid ไม่ตรง)
    [~, pkIdx]  = max(crProfile);
    xCenter     = xScene(pkIdx);

    % interpolate รอบ peak ±10 m ด้วย step 0.01 m
    xFine  = linspace(xCenter - 10, xCenter + 10, 2001);
    crFine = interp1(xScene, crProfile, xFine, 'spline');

    % วัด 3-dB width แล้วหาร 2
    halfMax_cr   = max(crFine) / 2;
    above_cr     = crFine >= halfMax_cr;
    measured_cr  = sum(above_cr) * (xFine(2)-xFine(1)) / 2;

    % Theoretical: δCR = λR / (2L)
    theory_cr = (lambda * targets_range(k)) / (2 * L);

    fprintf('  T%d (R=%4dm) | measured = %.3f m | theoretical = %.3f m\n', ...
        k, targets_range(k), measured_cr, theory_cr);
end

% ─────────────────────────────────────────────────────────────
% 3. PEAK LOCATION ACCURACY
%
%    ตรวจสอบว่า BP วาง target ถูกตำแหน่งไหม
%    โดยหา (x, y) ของ pixel ที่สว่างที่สุดในบริเวณแต่ละ target
%    แล้วเทียบกับ ground truth
% ─────────────────────────────────────────────────────────────

fprintf('\n--- Peak Location Accuracy ---\n');
fprintf('  (ground truth vs measured peak position)\n');

for k = 1:3
    R_true = targets_range(k);   % ground truth range
    x_true = 0;                  % ground truth cross-range (targets อยู่ที่ x=0)

    % หา peak ใน window รอบ target จริง ±50 m range, ±20 m cross-range
    yMask = abs(yScene - R_true) <= 50;
    xMask = abs(xScene - x_true) <= 20;
    subImg = bpImageMag(yMask, :);
    subImg = subImg(:, xMask);

    [~, flatIdx]  = max(subImg(:));
    [rSub, cSub]  = ind2sub(size(subImg), flatIdx);

    yVec = yScene(yMask);
    xVec = xScene(xMask);
    r_measured = yVec(rSub);
    x_measured = xVec(cSub);

    fprintf('  T%d | truth=(%4dm, 0m) | peak=(%4.0fm, %.0fm) | error=(%.0fm, %.0fm)\n', ...
        k, R_true, r_measured, x_measured, ...
        abs(r_measured - R_true), abs(x_measured - x_true));
end

% ─────────────────────────────────────────────────────────────
% ** ข้อจำกัดของตัวเลขในส่วน EVALUATION นี้ — อ่านก่อนนำไปอ้างอิง **
%
% 1) กริดหยาบเกินกว่าจะวัดได้แม่น
%    - แกน range  ห่าง 2.00 m แต่ mainlobe กว้าง ~1.5-3 m -> ค่าถูกปัดเป็นขั้นละ 1 m
%    - แกน cross  ห่าง 1.00 m แต่ mainlobe กว้าง ~0.09 m -> เล็กกว่าช่องกริด 10 เท่า
%      ค่า ~0.58 m ที่ได้คือความกว้างของ kernel ที่ใช้ interpolate ไม่ใช่ของเรดาร์
%    - peak location error ที่รายงานเป็น 0 m ก็เพราะกริดมองไม่เห็น error ระดับ cm
%
% 2) วัดคนละนิยามกับทฤษฎี
%    โค้ดข้างบนวัดที่ max/2 ของ "แอมพลิจูด" = -6 dB เชิงกำลัง
%    แต่นิยามมาตรฐานของ resolution คือความกว้างที่ -3 dB (ครึ่งกำลัง)
%    ซึ่งบนแกนแอมพลิจูดอยู่ที่ 0.7071 * peak
%
% ตัวเลขที่ถูกต้องวัดไว้แล้วใน  W1_3_SaveFigs_and_ExactEval.m
% (ยิง BP ซ้ำบนเส้นตัด 1 มิติ step 0.01 m / 0.002 m) ผลที่ได้:
%
%   range res (-3 dB) : T1 1.509 | T2 1.737 | T3 2.148  m   (ทฤษฎี c/2B = 3.000)
%       แคบกว่าทฤษฎีเพราะ BP รวมข้ามมุมกวาด 17-28 องศา ซึ่งขยาย support
%       ในแนว range ด้วย  ยิ่งไกลมุมยิ่งแคบ ค่ายิ่งลู่เข้า 3 m ตามตำรา
%
%   cross res (-3 dB) : T1 0.0785 | T2 0.0938 | T3 0.1167 m
%       เทียบทฤษฎี lambda*R/2L = 0.0749 | 0.0937 | 0.1218  -> ตรงภายใน 5%
%
%   peak error        : range <= 0.07 m , cross <= 0.002 m
% ─────────────────────────────────────────────────────────────

fprintf('\n>> รูปทั้งหมดเซฟไว้ที่: %s\n', figDir);
fprintf('>> หมายเหตุ: ตัวเลข resolution ข้างบนถูกจำกัดด้วยความละเอียดของกริด\n');
fprintf('   ค่าที่วัดแบบเป๊ะอยู่ใน W1_3_SaveFigs_and_ExactEval.m (ดูคอมเมนต์ท้ายไฟล์นี้)\n');

fprintf('\n=================================\n');