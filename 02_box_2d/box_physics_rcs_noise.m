%% COMBINED A — Core pipeline + Noise + Corner Localization
%  หมายเหตุ: บล็อก RADAR CONFIGURATION สืบทอดมาจาก 01_point_target/point3_bp_2d.m
%  ซึ่งดัดแปลงจากตัวอย่าง Stripmap SAR ของ MathWorks (ดู attribution ในไฟล์นั้น)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  ไฟล์รวมโค้ดจริง (ก๊อปจากไฟล์ย่อยมาต่อกัน) — แต่ละ section มี clear เอง
%  ทำงานอิสระต่อกัน; local function ทั้งหมดรวมไว้ท้ายไฟล์ (ชื่อไม่ชนกัน)
%  ที่มา: try104_PhysicsRCS.m, try105_AWGN_Noise.m, try109_CornerLocalization.m

%% ================= SECTION 1: CORE PIPELINE (try104)  (from try104_PhysicsRCS.m) =================
%% SAR Backprojection — Week 6 | Step 1: Physics-based RCS
%
%  เปลี่ยนจาก RCS ตั้งมือ → คำนวณจาก Flat Plate formula:
%      σ = 4πA²/λ²
%  โดย A = พื้นที่หน้าแต่ละด้าน, λ = c/fc
%
%  เปรียบเทียบกับ try103 (Week 5):
%    Week 5: faceRCS = [2.0, 0.5, 1.0, 1.0, 0.8, 0.3]  (ตั้งมือ)
%    Week 6: faceRCS = คำนวณจากสูตรฟิสิกส์ → ค่าจะใหญ่กว่ามาก (~10^6–10^8 m²)
%
%  Figures: 1–10 (pipeline order เหมือน try103)
%  Note fix: Range Profile (Fig 6 bottom-right) แก้สเกล dB axis

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

%% ===== STEP 1: PHYSICS-BASED RCS (Flat Plate Formula, sub-patch) =====
%
%  Plate broadside max:   sigma_plate = 4*pi*A^2 / lambda^2   (normal incidence ONLY)
%
%  Sub-patch discretization (each face = ptsPerFace point scatterers):
%    each scatterer represents a sub-area  a = A / ptsPerFace
%       sigma_i = 4*pi*a^2 / lambda^2 = sigma_plate / ptsPerFace^2
%    -> coherent in-phase sum at broadside reconstructs the full plate:
%       ( ptsPerFace * sqrt(sigma_i) )^2 = sigma_plate     (exact)
%
%  Face areas:
%    Front / Back  (YZ): A = boxW*boxH = 20*10 = 200 m^2
%    Left  / Right (XZ): A = boxL*boxH = 30*10 = 300 m^2
%    Top   / Bottom(XY): A = boxL*boxW = 30*20 = 600 m^2
%
%  LIMITATIONS (note for report):
%   - 4*pi*A^2/lambda^2 is the BROADSIDE MAXIMUM; this side/down-looking geometry
%     is never at normal incidence, so these are upper-bound RCS values.
%   - No aspect/occlusion model: the Bottom face (faces the ground) and the far
%     faces still return full RCS even though they are physically hidden.
%   - In a noise-free, self-normalized image the ABSOLUTE RCS scale cancels out;
%     only the RELATIVE face weighting changes the picture. Absolute level starts
%     to matter in Step 2 (AWGN), which is why the sub-patch level is used here.

lambda_calc = c / fc;            % 0.0749 m at 4 GHz
ptsPerFace  = 36;                % = N_face^2 scatterers per face (N_face = 6)

A_front_back = boxW * boxH;      % 200 m^2
A_left_right = boxL * boxH;      % 300 m^2
A_top_bottom = boxL * boxW;      % 600 m^2

% full-plate broadside RCS (reporting / physical reference)
splate_front_back = 4*pi * A_front_back^2 / lambda_calc^2;
splate_left_right = 4*pi * A_left_right^2 / lambda_calc^2;
splate_top_bottom = 4*pi * A_top_bottom^2 / lambda_calc^2;
sigmaPlate = [splate_front_back, splate_front_back, ...
              splate_left_right, splate_left_right, ...
              splate_top_bottom, splate_top_bottom];

% per-scatterer sub-patch RCS used in the simulation  (sigma_plate / ptsPerFace^2)
faceRCS = sigmaPlate / ptsPerFace^2;

fprintf('===== STEP 1: Physics-based RCS (flat plate, sub-patch) =====\n');
fprintf('  lambda = %.4f m  (fc = %.1f GHz)\n', lambda_calc, fc/1e9);
fprintf('  Scatterers/face = %d  ->  sub-area a = A/%d,  sigma_i = sigma_plate/%d^2\n\n', ...
    ptsPerFace, ptsPerFace, ptsPerFace);
fprintf('  %-12s %7s %16s %16s\n', 'Face', 'A(m^2)', 'plate sigma(m^2)', 'per-scat sigma');
fprintf('  %-12s %7.0f %16.3e %12.3e (%.1f dBsm)\n', ...
    'Front/Back', A_front_back, splate_front_back, faceRCS(1), 10*log10(faceRCS(1)));
fprintf('  %-12s %7.0f %16.3e %12.3e (%.1f dBsm)\n', ...
    'Left/Right', A_left_right, splate_left_right, faceRCS(3), 10*log10(faceRCS(3)));
fprintf('  %-12s %7.0f %16.3e %12.3e (%.1f dBsm)\n', ...
    'Top/Bottom', A_top_bottom, splate_top_bottom, faceRCS(5), 10*log10(faceRCS(5)));
fprintf('\n  Plate broadside dBsm : Front/Back=%.1f, Left/Right=%.1f, Top/Bottom=%.1f\n', ...
    10*log10(splate_front_back), 10*log10(splate_left_right), 10*log10(splate_top_bottom));
fprintf('  Week 5 (manual) RCS  : [2.0 0.5 1.0 1.0 0.8 0.3] m^2\n');
fprintf('  Relative weighting   : Front/Back : Left/Right : Top/Bottom = 1 : 2.25 : 9\n\n');

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

%% ===== STEP 1.5: STATIC OCCLUSION (face visibility) =====
%  ปัญหา: โมเดลเดิมให้ทุกหน้าคืน RCS เต็ม ทั้งที่หน้า Back (ด้านไกล) และ
%  Bottom (หันลงพื้น) เรดาร์มองไม่เห็น  -> ตัดหน้าที่ normal ไม่เคยหันเข้า
%  เรดาร์ตลอด aperture ออก (set RCS = 0)
%
%  เป็น static approximation: เช็คว่าหน้านั้น "เคยถูกเห็น" ระหว่างบินไหม
%  (per-pulse visibility + angular taper sigma(theta) = งานละเอียดขั้นต่อไป)
%  หมายเหตุ: occlusion ไม่กระทบ calibration (ใช้ point target แยก)

useOcclusion = true;   % << ปิดเป็น false เพื่อเทียบกับโมเดลเดิม

faceNormal = [ -1  0  0;    % Front  (-range, หันเข้าเรดาร์)
                1  0  0;    % Back   (+range)
                0 -1  0;    % Left   (-cross)
                0  1  0;    % Right  (+cross)
                0  0  1;    % Top    (+height)
                0  0 -1];   % Bottom (-height)

losAll = radarPosHistory - boxCenter;     % 3 x numpulses : เวกเตอร์ box -> radar
faceVisible = true(1,6);
for f = 1:6
    d = faceNormal(f,:) * losAll;         % 1 x numpulses dot products
    faceVisible(f) = any(d > 0);          % เห็นอย่างน้อย 1 pulse
end

if useOcclusion
    for f = 1:6
        idx = (f-1)*ptsPerFace + (1:ptsPerFace);
        scatRCS(idx) = scatRCS(idx) * faceVisible(f);   % หน้าที่ถูกบัง -> 0
    end
    fprintf('Occlusion ON  | visible: %s | hidden: %s\n', ...
        strjoin(faceLabels(faceVisible), ', '), strjoin(faceLabels(~faceVisible), ', '));
else
    fprintf('Occlusion OFF | all 6 faces illuminated (original model)\n');
end

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

%% ===== RANGE COMPRESSION (+ optional window for sidelobe control) =====
%  range sidelobe (PSLR) มาจากสเปกตรัมสี่เหลี่ยมของ LFM -> taper ตัว reference
%  chirp ก่อนทำ matched filter จะลด sidelobe (แลกกับ mainlobe กว้างขึ้น)
%  เพราะ LFM กวาดความถี่ตามเวลา (time ~ frequency) การ taper ตาม sample = taper สเปกตรัม

% FINDING (W6): ทดสอบแล้ว window ไม่ลด range PSLR ใน BP geometry นี้
%   (range-cut sidelobe ~-10 dB มาจาก 2D PSF coupling/depression angle
%    ไม่ใช่ matched-filter sidelobe) -> ตั้ง false ไว้ (ได้ res ดีกว่า, PSLR เท่าเดิม)
%   เปิด true เพื่อทดลองได้ แต่จะทำให้ resolution กว้างขึ้นโดยไม่ได้ PSLR กลับมา
useWindow = false;         % << true = ใส่ window (พิสูจน์แล้วไม่ช่วย geometry นี้)
winType   = 'taylor';      % 'taylor' | 'hamming'

refPulse      = waveform();
chirpSamples  = round(tpd * fs);
refChirp      = refPulse(1:chirpSamples);

if useWindow
    switch lower(winType)
        case 'hamming'
            rngWin = hamming(chirpSamples);
        case 'taylor'
            rngWin = taylorwin(chirpSamples, 4, -35);   % nbar=4, sidelobe -35 dB
        otherwise
            rngWin = ones(chirpSamples,1);
    end
    fprintf('Range window  : %s (sidelobe control ON)\n', winType);
else
    rngWin = ones(chirpSamples,1);
    fprintf('Range window  : none (rectangular)\n');
end

matchedFilter = conj(flipud(refChirp .* rngWin(:)));

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
legNames = cell(1,6);
for f = 1:6
    idx = (f-1)*ptsPerFace + (1:ptsPerFace);
    % marker size scaled by plate dBsm (physically meaningful face RCS)
    markerSz = 10*log10(sigmaPlate(f)) - 60;   % plate dBsm offset -> marker size
    if useOcclusion && ~faceVisible(f)
        mAlpha = 0.12;  vis = ' [hidden]';      % faded = occluded face
    else
        mAlpha = 0.9;   vis = '';
    end
    scatter3(scatPos(2,idx), scatPos(1,idx), scatPos(3,idx), ...
        max(markerSz, 10), ...
        repmat(faceColors(f,:), ptsPerFace, 1), 'filled', ...
        'MarkerFaceAlpha', mAlpha, ...
        'DisplayName', sprintf('%s (\\sigma_{plate}=%.1f dBsm)%s', faceLabels{f}, 10*log10(sigmaPlate(f)), vis));
    legNames{f} = sprintf('%s%s', faceLabels{f}, vis);
end
xlabel('Cross-Range (m)'); ylabel('Range (m)'); zlabel('Height (m)');
title('Fig 1 | 3D Box Target — Physics RCS + occlusion (W6)');
legend('Box (transparent)', legNames{:}, 'Location','northeast');
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

if useWindow
    winStr = sprintf('%s window', winType);
else
    winStr = 'no window';
end
sgtitle(sprintf('Fig 12 | Point-Target Calibration — true system IRF (single scatterer, noise-free, %s)', winStr));

%% ===== SAVE FIGURES =====

figDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figure');
if ~exist(figDir, 'dir'); mkdir(figDir); end
for fn = [1:10, 12]   % Fig 11 (eval text) is built later and saved after creation
    fh = figure(fn);
    exportgraphics(fh, fullfile(figDir, sprintf('fig%02d.png', fn)), 'Resolution', 150);
end
fprintf('Figures saved to %s\n', figDir);

%% ===== EVALUATION =====

fprintf('\n========== EVALUATION — Week 6 Step 1: Physics-based RCS ==========\n');

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
if useOcclusion
    fprintf('  Occlusion : ON  -> visible faces: %s  (hidden: %s)\n', ...
        strjoin(faceLabels(faceVisible), ', '), strjoin(faceLabels(~faceVisible), ', '));
else
    fprintf('  Occlusion : OFF -> all 6 faces illuminated\n');
end

fprintf('\n[RCS Comparison — W5 vs W6]\n');
fprintf('  W6 = full-plate broadside sigma; simulation uses sub-patch (sigma/%d^2 per scatterer)\n', ptsPerFace);
fprintf('  %-12s  %-12s  %-18s  %-10s\n', 'Face', 'W5 (m^2)', 'W6 plate (m^2)', 'W6 (dBsm)');
w5vals = [2.0, 0.5, 1.0, 1.0, 0.8, 0.3];
for f = 1:6
    fprintf('  %-12s  %-12.1f  %-18.3e  %-10.1f\n', ...
        faceLabels{f}, w5vals(f), sigmaPlate(f), 10*log10(sigmaPlate(f)));
end

fprintf('\n[Theoretical Resolution]\n');
rangeRes_th = c / (2 * bw);
crRes_th    = lambda_e * boxCenter(1) / (2 * apertureLen);
fprintf('  Range      : %.2f m\n', rangeRes_th);
fprintf('  Cross-range: %.4f m (%.2f cm)\n', crRes_th, crRes_th*100);

fprintf('\n[Point-Target Calibration — TRUE system resolution & PSLR]\n');
fprintf('  (from a single scatterer''s impulse response, NOT the box)\n');
if useWindow
    fprintf('  Range window     : %s ON  (lower range PSLR, wider main lobe)\n', winType);
else
    fprintf('  Range window     : none (rectangular)\n');
end
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

fprintf('\n[Step 1 / 1.5 Summary]\n');
fprintf('  - RCS from 4*pi*A^2/lambda^2; relative weighting Front/Back:Left/Right:\n');
fprintf('    Top/Bottom = 1:2.25:9. In a noise-free self-normalized image the\n');
fprintf('    ABSOLUTE scale cancels -- only this RELATIVE weighting changes things.\n');
fprintf('  - Step 1.5 occlusion: Back & Bottom faces culled (never face radar);\n');
fprintf('    Top now dominates alone -> a more honest, less symmetric image.\n');
fprintf('  - Resolution/PSLR depend on BW/aperture, NOT RCS (see calibration).\n');
fprintf('  - Absolute RCS level will matter in Step 2 once AWGN is added.\n');
fprintf('  Next: Step 2 -- add AWGN noise\n');

%% ===== FIGURE 11: Evaluation Text (Command Window output) =====

% Build text string from computed variables
w5vals_fig = [2.0, 0.5, 1.0, 1.0, 0.8, 0.3];

evalLines = {};
evalLines{end+1} = '====== EVALUATION — Week 6 Step 1: Physics-based RCS ======';
evalLines{end+1} = '';
evalLines{end+1} = '[Radar Parameters]';
evalLines{end+1} = sprintf('  Center freq : %.1f GHz   BW: %.1f MHz', fc/1e9, bw/1e6);
evalLines{end+1} = sprintf('  Wavelength  : %.4f m      PRF: %d Hz', lambda_e, prf);
evalLines{end+1} = sprintf('  Speed: %.0f m/s   Height: %.0f m   Aperture: %.0f m', speed, radarpos0(3), apertureLen);
evalLines{end+1} = '';
evalLines{end+1} = '[Box Target]';
evalLines{end+1} = sprintf('  Center: range=%.0f m, cross=%.0f m', boxCenter(1), boxCenter(2));
evalLines{end+1} = sprintf('  Size  : %.0f x %.0f x %.0f m   Scatterers: %d', boxL, boxW, boxH, numScatterers);
if useOcclusion
    evalLines{end+1} = sprintf('  Occlusion: visible %s', strjoin(faceLabels(faceVisible), ','));
    evalLines{end+1} = sprintf('             hidden  %s', strjoin(faceLabels(~faceVisible), ','));
else
    evalLines{end+1} = '  Occlusion: OFF (all 6 faces)';
end
evalLines{end+1} = '';
evalLines{end+1} = '[RCS Comparison — W5 vs W6]';
evalLines{end+1} = sprintf('  %-8s  %-10s  %-16s  %s', 'Face','W5 (m^2)','W6 plate(m^2)','W6 (dBsm)');
evalLines{end+1} = repmat('-', 1, 54);
for f = 1:6
    evalLines{end+1} = sprintf('  %-8s  %-10.1f  %-16.3e  %.1f', ...
        faceLabels{f}, w5vals_fig(f), sigmaPlate(f), 10*log10(sigmaPlate(f)));
end
evalLines{end+1} = '';
evalLines{end+1} = '[Theoretical Resolution]';
evalLines{end+1} = sprintf('  Range      : %.2f m  (c/2B)', rangeRes_th);
evalLines{end+1} = sprintf('  Cross-range: %.4f m = %.2f cm  (lambdaR/2L)', crRes_th, crRes_th*100);
evalLines{end+1} = '';
evalLines{end+1} = '[Calibration — TRUE resolution (single point, not box)]';
if useWindow
    evalLines{end+1} = sprintf('  Range window: %s ON (sidelobe control)', winType);
else
    evalLines{end+1} = '  Range window: none (rectangular)';
end
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
evalLines{end+1} = '[Step 1 Summary — W5 vs W6]';
evalLines{end+1} = sprintf('  TRUE resolution (calib): range %.2f m, cross %.3f m', calRangeRes, calCrossRes);
evalLines{end+1} = '    -> set by BW & aperture; RCS does NOT change it.';
if ~isempty(above3rE)
    mR = max(above3rE)-min(above3rE);
    evalLines{end+1} = sprintf('  Box range extent : %.2f -> %.2f m  (now ~ box L)', W5_rangeRes, mR);
end
if exist('snr_dBE','var')
    evalLines{end+1} = sprintf('  Peak/background  : %.1f -> %.1f dB  (%+.1f dB)', W5_snr, snr_dBE, snr_dBE-W5_snr);
end
if exist('yCentE','var')
    leR = abs(yCentE-boxCenter(1)); leC = abs(xCentE-boxCenter(2));
    evalLines{end+1} = sprintf('  Loc err          : R %.2f->%.2f m,  X %.2f->%.2f m', W5_locErrR, leR, W5_locErrC, leC);
end
evalLines{end+1} = '  Cause of any change: RELATIVE face weighting only; the';
evalLines{end+1} = '  noise-free image is scale-invariant in absolute RCS.';
evalLines{end+1} = '  Eval fix: resolution/PSLR now from point-target calibration,';
evalLines{end+1} = '  box used only for extent/localization/dynamic-range.';
evalLines{end+1} = '  Next: Step 2 — AWGN noise';
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


%% ================= SECTION 2: AWGN NOISE (try105)  (from try105_AWGN_Noise.m) =================
%%   (figure numbers offset by +100 so all sections' figures stay open)
%% SAR Backprojection — Week 6 | Step 2: AWGN Noise Study
%
%  ศึกษาผลของ noise: ใส่ controlled AWGN ที่ระดับ per-pulse raw echo
%  แล้ว sweep input SNR = [+20 +10 0 -10 -20] dB  (เทียบกับ noise-free)
%
%  วิธี: bypass receiver (ใช้ collector output = clean) แล้วบวก complex AWGN
%        ที่กำหนดกำลังเอง -> คุม SNR เป็น dB ตรงๆ และ reproducible (rng)
%
%  Key concept — coherent integration / processing gain:
%     image SNR  ~=  input SNR  +  range gain (10log10(T*B))
%                                +  azimuth gain (10log10(Npulses))
%  -> SAR ทน noise ได้ดีมาก เพราะรวมพลังงานหลายพัน pulse แบบ coherent
%
%  Forward sim รันครั้งเดียว; BP รันแค่ 2 ครั้ง (clean + unit-noise) เพราะ
%  BP เป็น linear -> ทุก SNR ได้จากการ scale + บวก field เดิม (เร็ว ~2-4 นาที)

clear;   % combined: keep prior figures open
rng(2024);   % reproducible noise realization

%% ===== RADAR CONFIG (เหมือน try104) =====
c  = physconst('LightSpeed'); fc = 4e9;
rangeResolution = 3; bw = c/(2*rangeResolution);
prf = 1000; aperture = 4; tpd = 3e-6; fs = 120e6;
waveform = phased.LinearFMWaveform('SampleRate',fs,'PulseWidth',tpd,'PRF',prf,'SweepBandwidth',bw);
speed = 100; flightDuration = 4;
slowTime = 1/prf; numpulses = flightDuration/slowTime + 1;
maxRange = 2500; truncrangesamples = ceil((2*maxRange/c)*fs);
antenna     = phased.CosineAntennaElement('FrequencyRange',[1e9 6e9]);
antennaGain = aperture2gain(aperture, c/fc);
transmitter = phased.Transmitter('PeakPower',50e3,'Gain',antennaGain);
radiator    = phased.Radiator('Sensor',antenna,'OperatingFrequency',fc,'PropagationSpeed',c);
collector   = phased.Collector('Sensor',antenna,'PropagationSpeed',c,'OperatingFrequency',fc);
channel     = phased.FreeSpace('PropagationSpeed',c,'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

%% ===== BOX + sub-patch RCS + occlusion (เหมือน try104) =====
boxCenter = [1000;0;0];
boxL=30; halfL=15; boxW=20; halfW=10; boxH=10; halfH=5;
lambda = c/fc; ptsPerFace = 36;
sigmaPlate = [4*pi*(boxW*boxH)^2/lambda^2*[1 1], ...
              4*pi*(boxL*boxH)^2/lambda^2*[1 1], ...
              4*pi*(boxL*boxW)^2/lambda^2*[1 1]];
faceRCS = sigmaPlate / ptsPerFace^2;

N_face = 6; u = linspace(-1,1,N_face);
scatPos = []; scatRCS = [];
faceLabels = {'Front','Back','Left','Right','Top','Bottom'};
[V,W] = meshgrid(u*halfW, u*halfH);
scatPos=[scatPos,[(boxCenter(1)-halfL)*ones(1,N_face^2);(boxCenter(2)+V(:))';(boxCenter(3)+halfH+W(:))']]; scatRCS=[scatRCS,faceRCS(1)*ones(1,N_face^2)];
scatPos=[scatPos,[(boxCenter(1)+halfL)*ones(1,N_face^2);(boxCenter(2)+V(:))';(boxCenter(3)+halfH+W(:))']]; scatRCS=[scatRCS,faceRCS(2)*ones(1,N_face^2)];
[U,W] = meshgrid(u*halfL, u*halfH);
scatPos=[scatPos,[(boxCenter(1)+U(:))';(boxCenter(2)-halfW)*ones(1,N_face^2);(boxCenter(3)+halfH+W(:))']]; scatRCS=[scatRCS,faceRCS(3)*ones(1,N_face^2)];
scatPos=[scatPos,[(boxCenter(1)+U(:))';(boxCenter(2)+halfW)*ones(1,N_face^2);(boxCenter(3)+halfH+W(:))']]; scatRCS=[scatRCS,faceRCS(4)*ones(1,N_face^2)];
[U,V] = meshgrid(u*halfL, u*halfW);
scatPos=[scatPos,[(boxCenter(1)+U(:))';(boxCenter(2)+V(:))';(boxCenter(3)+boxH)*ones(1,N_face^2)]]; scatRCS=[scatRCS,faceRCS(5)*ones(1,N_face^2)];
scatPos=[scatPos,[(boxCenter(1)+U(:))';(boxCenter(2)+V(:))';zeros(1,N_face^2)]]; scatRCS=[scatRCS,faceRCS(6)*ones(1,N_face^2)];
numScatterers = size(scatPos,2);

radarpos0=[0;-200;500]; radarvel0=[0;speed;0];
radarPosHistory = radarpos0 + radarvel0*((1:numpulses)-1)/prf;

% Step 1.5 occlusion: cull faces never facing the radar
faceNormal=[-1 0 0;1 0 0;0 -1 0;0 1 0;0 0 1;0 0 -1];
losAll = radarPosHistory - boxCenter; faceVisible = true(1,6);
for f=1:6, faceVisible(f) = any(faceNormal(f,:)*losAll > 0); end
for f=1:6, idx=(f-1)*ptsPerFace+(1:ptsPerFace); scatRCS(idx)=scatRCS(idx)*faceVisible(f); end
fprintf('Occlusion: visible %s | hidden %s\n', ...
    strjoin(faceLabels(faceVisible),','), strjoin(faceLabels(~faceVisible),','));

%% ===== FORWARD SIM — clean echo, NO receiver noise (run ONCE) =====
fprintf('Simulating clean raw echo (once, ~1 min)...\n');
boxPlatform   = phased.Platform('InitialPosition',scatPos,'Velocity',zeros(size(scatPos)));
boxTarget     = phased.RadarTarget('OperatingFrequency',fc,'MeanRCS',scatRCS);
radarPlatform = phased.Platform('InitialPosition',radarpos0,'Velocity',radarvel0);
rxsigClean = zeros(truncrangesamples,numpulses);
for ii=1:numpulses
    [radarpos,radarvel]=radarPlatform(slowTime);
    [tpos,tvel]=boxPlatform(slowTime);
    [~,tAng]=rangeangle(tpos,radarpos);
    sig=waveform(); sig=sig(1:truncrangesamples); sig=transmitter(sig);
    tAng(1,:)=zeros(1,numScatterers);
    sig=radiator(sig,tAng); sig=channel(sig,radarpos,tpos,radarvel,tvel); sig=boxTarget(sig);
    rxsigClean(:,ii)=collector(sig,tAng);   % <-- ไม่ผ่าน receiver = clean
end
fprintf('Clean echo done.\n');

% matched filter (no window)
refChirp=waveform(); chirpSamples=round(tpd*fs); refChirp=refChirp(1:chirpSamples);
matchedFilter=conj(flipud(refChirp)); mfLen=length(matchedFilter);

% reference signal power = mean per-pulse PEAK power (ใช้นิยาม input SNR)
pulsePeakPow = max(abs(rxsigClean).^2, [], 1);
sigRefPow    = mean(pulsePeakPow(pulsePeakPow>0));

% imaging grid
xScene=linspace(-100,100,401); yScene=linspace(850,1200,401);
[xGrid,yGrid]=meshgrid(xScene,yScene); Ny=numel(yScene); Nx=numel(xScene);
cornerMask=(abs(xGrid-boxCenter(2))>halfW*5)&(abs(yGrid-boxCenter(1))>halfL*5);
boxX_cr=boxCenter(2)+halfW*[-1 1 1 -1 -1]; boxY_r=boxCenter(1)+halfL*[-1 -1 1 1 -1];

procGain = 10*log10(numpulses) + 10*log10(tpd*bw);
fprintf('Processing gain ~ %.1f dB (azimuth %.1f + range %.1f)\n', ...
    procGain, 10*log10(numpulses), 10*log10(tpd*bw));

%% ===== BP DECOMPOSITION (linear) — run BP only TWICE =====
%  BP เป็น linear operator: BP(clean + noise) = BP(clean) + BP(noise)
%  ใช้ noise realization เดียว (กำลัง = sigRefPow คือ input SNR 0 dB) แล้ว
%  scale field ตาม SNR -> รัน BP แค่ 2 ครั้ง (clean + unit-noise) แม่นยำเป๊ะ
fprintf('Backprojection x2 (clean + unit-noise)...\n');
rc0 = rangeCompress(rxsigClean, matchedFilter);
bp0 = backproject(rc0, radarPosHistory, xGrid, yGrid, Ny, Nx, c, fs, fc, mfLen);
bp0mag = abs(bp0)/max(abs(bp0(:)));
[snr0, loc0, ~] = imgMetrics(bp0mag, xGrid, yGrid, cornerMask, boxCenter);
fprintf('  noise-free -> image SNR=%.1f dB, locErr=%.2f m\n', snr0, loc0);

unitNoise = sqrt(sigRefPow/2)*(randn(size(rxsigClean))+1i*randn(size(rxsigClean)));  % 0 dB power
rcN = rangeCompress(unitNoise, matchedFilter);
bpN = backproject(rcN, radarPosHistory, xGrid, yGrid, Ny, Nx, c, fs, fc, mfLen);

%% ===== SNR SWEEP =====
% example images use the first 5; curve extends lower to find the break point
snrList = [20 10 0 -10 -20 -30 -40 -50];   % input per-pulse SNR (dB)
nSNR    = numel(snrList);
imgSNR=zeros(1,nSNR); locErr=zeros(1,nSNR); bgFloor=zeros(1,nSNR); detected=false(1,nSNR);
bpStore=cell(1,nSNR);

for k=1:nSNR
    snr = snrList(k);
    bp  = bp0 + 10^(-snr/20) * bpN;     % scale unit-noise field to this input SNR
    bpMag = abs(bp)/max(abs(bp(:))); bpStore{k}=bpMag;
    [imgSNR(k), locErr(k), pkOK] = imgMetrics(bpMag, xGrid, yGrid, cornerMask, boxCenter);
    bgFloor(k) = mean(bpMag(cornerMask));
    detected(k) = (imgSNR(k) > 13) && pkOK;   % peak on box AND >13 dB above bg
    fprintf('  SNRin=%+3d dB -> image SNR=%5.1f dB | locErr=%5.2f m | detect=%d\n', ...
        snr, imgSNR(k), locErr(k), detected(k));
end

%% ===== FIGURE 1: Summary curves =====
figure(100+1); set(gcf,'Name','Fig1: AWGN Sweep Summary','Position',[80 80 980 420]);

subplot(1,2,1);
plot(snrList, imgSNR, 'o-','LineWidth',2,'MarkerFaceColor','b','DisplayName','Image SNR (measured)'); hold on;
yline(snr0,'k--','LineWidth',1.5,'DisplayName','noise-free baseline');
plot(snrList, snrList+procGain,'r:','LineWidth',1.5,'DisplayName','input + proc.gain (theory bound)');
yline(13,'g--','LineWidth',1,'DisplayName','detection threshold (13 dB)');
xlabel('Input SNR per pulse (dB)'); ylabel('Image SNR (dB)');
title('Image SNR vs Input SNR'); grid on; legend('Location','southeast');
set(gca,'XDir','reverse');   % noisy (low SNR) ไปทางขวา

subplot(1,2,2);
yyaxis left;  plot(snrList, locErr,'o-','LineWidth',2); ylabel('Localization error (m)');
yyaxis right; stem(snrList, double(detected),'filled','LineWidth',1.5); ylabel('Detected (1/0)'); ylim([-0.1 1.2]);
xlabel('Input SNR per pulse (dB)'); title('Localization & Detection vs Input SNR');
grid on; set(gca,'XDir','reverse');

sgtitle('Fig 1 | Step 2: AWGN Sweep — coherent integration makes SAR robust to noise');

%% ===== FIGURE 2: Example SAR images per SNR =====
figure(100+2); set(gcf,'Name','Fig2: SAR Images vs SNR','Position',[60 60 1200 720],'Color','k');
tl=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

panels = {bp0mag, bpStore{1}, bpStore{2}, bpStore{3}, bpStore{4}, bpStore{5}};
labels = {'noise-free', 'SNRin = +20 dB','SNRin = +10 dB','SNRin = 0 dB','SNRin = -10 dB','SNRin = -20 dB'};
for p=1:6
    ax=nexttile;
    imagesc(xScene,yScene,20*log10(panels{p}+eps)); set(ax,'YDir','normal');
    clim([-40 0]); colormap(ax,'jet'); hold on;
    plot(boxX_cr,boxY_r,'w--','LineWidth',1.2);
    title(labels{p},'Color','w');
    set(ax,'Color','k','XColor','w','YColor','w');
    if p==1||p==4, ylabel('Range (m)','Color','w'); end
    if p>=4, xlabel('Cross-Range (m)','Color','w'); end
end
cb=colorbar; cb.Color='w'; ylabel(cb,'dB','Color','w'); cb.Layout.Tile='east';
sgtitle('Fig 2 | SAR BP Image (dB) at decreasing SNR — box stays visible far below 0 dB','Color','w');

%% ===== FIGURE 3: Evaluation text =====
L={};
L{end+1}='====== STEP 2 — AWGN NOISE STUDY ======';
L{end+1}='';
L{end+1}='[Setup]';
L{end+1}=sprintf('  Controlled AWGN at raw-echo level; rng(2024) reproducible');
L{end+1}=sprintf('  Pulses: %d   Processing gain ~ %.1f dB', numpulses, procGain);
L{end+1}=sprintf('    (azimuth %.1f dB + range comp %.1f dB)', 10*log10(numpulses), 10*log10(tpd*bw));
L{end+1}='';
L{end+1}='[Image SNR vs Input SNR]';
L{end+1}=sprintf('  %-12s  %-12s  %-10s  %s','Input SNR','Image SNR','LocErr','Detect');
L{end+1}=repmat('-',1,48);
L{end+1}=sprintf('  %-12s  %-12.1f  %-10.2f  %s','noise-free',snr0,loc0,'yes');
for k=1:nSNR
    L{end+1}=sprintf('  %-+12d  %-12.1f  %-10.2f  %s', snrList(k), imgSNR(k), locErr(k), ternary(detected(k),'yes','NO'));
end
L{end+1}='';
L{end+1}='[Reading]';
L{end+1}='  - Image SNR >> input SNR thanks to coherent integration';
L{end+1}='    (processing gain). Box detectable even at strongly negative';
L{end+1}='    per-pulse SNR.';
L{end+1}='  - Resolution/PSLR (system) unchanged by noise -> see calibration';
L{end+1}='    in try104 (noise affects detectability, not resolution).';
L{end+1}='  - Localization error rises only when noise overwhelms the';
L{end+1}='    bright-region threshold.';
L{end+1}=repmat('=',1,48);

figure(100+3); set(gcf,'Name','Fig3: Step2 Evaluation Text','Color','k','Position',[100 80 720 560]);
ax=axes('Position',[0 0 1 1],'Visible','off','Color','k'); hold(ax,'on');
ny=numel(L); ys=1/(ny+2);
for li=1:ny
    s=L{li};
    if startsWith(s,'==')||startsWith(s,'--'), clr=[.6 .6 .6];
    elseif startsWith(s,'['), clr=[.3 .85 1];
    elseif contains(s,' NO'), clr=[1 .4 .4];
    else, clr=[.92 .92 .92]; end
    text(0.03,1-li*ys,s,'Units','normalized','Color',clr,'FontName','Courier','FontSize',9,'Interpreter','none','Parent',ax);
end
hold(ax,'off');

%% ===== SAVE FIGURES =====
figDir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'figure');
if ~exist(figDir,'dir'), mkdir(figDir); end
for fn=1:3
    exportgraphics(figure(100+fn), fullfile(figDir,sprintf('fig_step2_%02d.png',fn)),'Resolution',150);
end
fprintf('\nStep 2 done. Figures saved to %s\n', figDir);


%% ================= SECTION 3: CORNER LOCALIZATION (try109)  (from try109_CornerLocalization.m) =================
%%   (figure numbers offset by +200 so all sections' figures stay open)
%% Week 6 (addendum) | TRUE localization error via corner reflectors
%
%  ปัญหาเดิม: "location error" วัดจาก centroid ของกล่อง = ปนรูปทรง/RCS
%             ไม่ใช่ความแม่นในการระบุตำแหน่งจริง (เรียกผิด ควรเป็น centroid offset)
%
%  วิธีที่ถูก (มาตรฐาน SAR): วาง point/corner reflector ที่พิกัดที่รู้ค่า
%  -> image -> หา PEAK ของแต่ละตัว (sub-pixel) -> วัดว่าห่างจากพิกัดจริงเท่าไหร่
%  = localization error ที่แท้จริง (คาดว่าระดับซับเมตร)
%
%  ใช้ analytic range-compressed model (sinc @ slant range + azimuth phase) ตรวจแล้ว

clear;   % combined: keep prior figures open
c=physconst('LightSpeed'); fc=4e9; lambda=c/fc;
bw=c/(2*3); fs=120e6; prf=1000; v=100; h=500; flightDur=4;
Naz=flightDur*prf+1; t=linspace(0,flightDur,Naz); ry=-200+v*t;
rngres=c/(2*bw);

% ----- corner reflectors at the box's 4 GROUND corners (known coords) -----
%   box: range 1000±15 (985/1015), cross 0±10 (-10/+10), ground z=0
corners=[985 -10; 985 10; 1015 -10; 1015 10];   % [groundRange, cross]
Ncorner=size(corners,1);

% slant-range axis
rbins=(1090:c/2/fs:1170).'; Nr=numel(rbins);

% build range-compressed data (each reflector = sinc @ its slant range)
S=zeros(Nr,Naz);
for k=1:Ncorner
    xg=corners(k,1); yg=corners(k,2);
    R=sqrt(xg^2+(ry-yg).^2+h^2);
    S=S+sinc((rbins-R)/rngres).*exp(-1j*4*pi*R/lambda);
end

% ----- backprojection on a FINE ground grid -----
gr=linspace(975,1025,401);  cr=linspace(-20,20,401);   % ~0.125 / 0.10 m pixels
[CR,GR]=meshgrid(cr,gr); Ny=numel(gr); Nx=numel(cr);
bp=zeros(Ny,Nx); kc=4*pi*fc/c;
for a=1:Naz
    Rp=sqrt(GR.^2+(CR-ry(a)).^2+h^2);
    iv=interp1(rbins,S(:,a),Rp(:),'linear',0);
    bp=bp+reshape(iv,Ny,Nx).*exp(1j*kc*Rp);
end
bpMag=abs(bp)/max(abs(bp(:)));

% ----- locate each corner: peak + parabolic sub-pixel refine -----
fprintf('\n===== TRUE Localization Error (corner reflectors) =====\n');
fprintf('  %-22s %-12s %-12s %-10s\n','corner (true)','found','err(range)','err(cross)');
errs=zeros(Ncorner,1); found=zeros(Ncorner,2);
for k=1:Ncorner
    xg=corners(k,1); yg=corners(k,2);
    % search window +/-6 m around true corner
    win=(abs(GR-xg)<6)&(abs(CR-yg)<6);
    M=bpMag; M(~win)=0;
    [~,idx]=max(M(:)); [ri,ci]=ind2sub([Ny Nx],idx);
    % parabolic sub-pixel refinement (range = rows, cross = cols)
    rr=subpix(gr, M(:,ci), ri);
    cc=subpix(cr.', M(ri,:).', ci);
    found(k,:)=[rr cc];
    er=abs(rr-xg); ec=abs(cc-yg);
    errs(k)=hypot(er,ec);
    fprintf('  (%4.0f,%+4.0f) m         (%6.1f,%+5.1f) %-12.3f %-10.3f\n',xg,yg,rr,cc,er,ec);
end
fprintf('  ---------------------------------------------------------\n');
fprintf('  RMS localization error : %.3f m   (max %.3f m)\n',sqrt(mean(errs.^2)),max(errs));
fprintf('  -> compare with "centroid offset" of the box ~ 8 m (a SHAPE metric, not this)\n');

% ----- figure -----
figure(200+1); set(gcf,'Name','Corner localization','Color','k','Position',[80 80 760 700]);
imagesc(cr,gr,20*log10(bpMag+eps)); axis xy image; clim([-30 0]); colormap(jet);
cb=colorbar; cb.Color='w'; ylabel(cb,'dB','Color','w'); hold on;
plot(corners(:,2),corners(:,1),'wo','MarkerSize',14,'LineWidth',1.5);
plot(found(:,2),found(:,1),'w+','MarkerSize',12,'LineWidth',1.5);
for k=1:Ncorner
    text(corners(k,2),corners(k,1)+1.6,sprintf('%.2f m',errs(k)),'Color','w','FontSize',10,'HorizontalAlignment','center');
end
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('Cross-Range (m)','Color','w'); ylabel('Ground Range (m)','Color','w');
title({'Corner-reflector localization (o = true,  + = measured peak)',...
       sprintf('RMS error = %.3f m  <<  box centroid offset ~8 m',sqrt(mean(errs.^2)))},'Color','w');

figDir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'figure'); if ~exist(figDir,'dir'),mkdir(figDir);end
exportgraphics(figure(200+1),fullfile(figDir,'fig_corner_localization.png'),'Resolution',150);
fprintf('Figure saved.\n');

%% ---- local function: parabolic sub-pixel peak position ----

%%%% ===== LOCAL FUNCTIONS (all sections) =====

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

function rc = rangeCompress(rx, mf)
    rc = zeros(size(rx,1)+length(mf)-1, size(rx,2));
    for ii=1:size(rx,2)
        rc(:,ii) = conv(rx(:,ii), mf, 'full');
    end
end


function bp = backproject(rc, rph, xGrid, yGrid, Ny, Nx, c, fs, fc, mfLen)
    pixelPos=[yGrid(:)'; xGrid(:)'; zeros(1,numel(xGrid))];
    bp=zeros(Ny,Nx); np=size(rph,2);
    for ii=1:np
        rp=rph(:,ii); dv=pixelPos-rp; sr=sqrt(sum(dv.^2,1));
        si=2*sr/c*fs + mfLen;
        vm=(si>=1)&(si<=size(rc,1)); iv=zeros(1,numel(xGrid));
        if any(vm), iv(vm)=interp1(1:size(rc,1), rc(:,ii), si(vm), 'linear', 0); end
        bp=bp+reshape(iv.*exp(1j*4*pi*fc/c.*sr), Ny, Nx);
    end
end


function [snrdB, locErr, peakOK] = imgMetrics(bpMag, xGrid, yGrid, cornerMask, boxCenter)
    bg = mean(bpMag(cornerMask));
    snrdB = 20*log10(max(bpMag(:))/bg);
    % peak location (is the brightest pixel on the box, or a noise spike?)
    [~, pk] = max(bpMag(:));
    peakOK = (abs(yGrid(pk)-boxCenter(1)) < 40) & (abs(xGrid(pk)-boxCenter(2)) < 40);
    bm = bpMag > 0.5;
    if any(bm(:))
        yC=mean(yGrid(bm)); xC=mean(xGrid(bm));
        locErr=hypot(yC-boxCenter(1), xC-boxCenter(2));
    else
        locErr=NaN;
    end
end


function out = ternary(cond, a, b)
    if cond, out=a; else, out=b; end
end

function xp = subpix(ax, prof, ip)
    prof=abs(prof);
    if ip<=1 || ip>=numel(prof), xp=ax(ip); return; end
    L=prof(ip-1); Cc=prof(ip); Rr=prof(ip+1);
    d=(L-2*Cc+Rr);
    if d==0, xp=ax(ip); else
        off=0.5*(L-Rr)/d;                 % fractional offset in samples
        xp=ax(ip)+off*(ax(ip+1)-ax(ip));
    end
end
