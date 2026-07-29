%% W1-3 ADD-ON — เซฟรูป pipeline + วัด resolution แบบ "เป๊ะ" ไม่ใช่ประมาณ
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%
%  วิธีใช้:
%     1. รัน  W1-3_3point.m  ให้จบก่อน (อย่าปิดรูป อย่า clear)
%     2. รันไฟล์นี้ต่อ
%
%  ไฟล์นี้ไม่แก้อะไรใน W1-3_3point.m เลย แค่ใช้ตัวแปรที่ค้างอยู่ใน workspace
%
%  =====================================================================
%  ทำ 3 อย่าง
%
%  PART A — เซฟรูป 5 ขั้นของ pipeline ลง figure/  (สำหรับใส่สไลด์)
%           stage 1 สร้างใหม่ (chirp + geometry) เพราะเดิมไม่มีรูปขั้นนี้
%
%  PART B — วัด resolution ใหม่บน "กริดละเอียด 1 มิติ"
%           เหตุผล: กริดของภาพหลักหยาบเกินกว่าจะวัดได้แม่น
%             - แกน range  ห่าง 2.00 m  แต่ mainlobe กว้างแค่ ~3 m -> ได้แค่ 1-2 จุด
%             - แกน cross  ห่าง 1.00 m  แต่ mainlobe กว้าง ~0.09 m -> วัดไม่ได้เลย
%           วิธีแก้: ยิง BP ใหม่เฉพาะ "เส้นตัด 1 มิติ" ผ่านเป้าแต่ละตัว
%           ด้วย step 0.005-0.01 m -> ได้ IRF จริงไม่โดนกริดบัง (เร็ว ไม่กี่วินาที)
%
%  PART C — พิมพ์สรุปตัวเลขเป็นบล็อกเดียว ก๊อปไปใช้ได้เลย
%
%  =====================================================================
%  หมายเหตุสำคัญเรื่องนิยาม "ความกว้าง"  <-- อ่านก่อน
%
%  โค้ดเดิมใน W1-3_3point.m วัดที่ระดับ  max/2  ของ "แอมพลิจูด"
%  ซึ่ง = -6 dB ในเชิงกำลัง ไม่ใช่ -3 dB
%
%  นิยามมาตรฐานของ resolution คือความกว้างที่ระดับ "ครึ่งกำลัง" (-3 dB)
%  ซึ่งบนแกนแอมพลิจูดคือ  peak / sqrt(2) = 0.7071 * peak
%
%  สคริปต์นี้รายงาน "ทั้งสองค่า" เพื่อให้เทียบกับทฤษฎีได้ตรงนิยาม
%  และให้รู้ว่าตัวเลขเดิมมาจากนิยามไหน
%  =====================================================================

if ~exist('bpImageMag','var') || ~exist('rxsigRC','var')
    error(['ไม่พบตัวแปรจาก W1-3_3point.m ใน workspace — ' ...
           'กรุณารัน W1-3_3point.m ให้จบก่อน แล้วค่อยรันไฟล์นี้']);
end

figDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'figure');
if ~exist(figDir,'dir'), mkdir(figDir); end

fprintf('\n');
fprintf('==============================================================\n');
fprintf(' W1-3 ADD-ON : save pipeline figures + exact resolution\n');
fprintf('==============================================================\n');

%% ====================================================================
%  PART A — เซฟรูป 5 ขั้นของ pipeline
%  ====================================================================

% ---- stage 1: Radar Config (สร้างใหม่) --------------------------------
% ซ้าย  : LFM chirp — ความถี่ทันทีไล่ขึ้นเป็นเส้นตรงตลอดพัลส์
% ขวา   : เรขาคณิตการบิน 3 มิติ (เส้นทางบิน + เป้า 3 จุด)

fp1 = figure('Name','Pipeline 1 — Radar Config','Color','w', ...
             'Position',[80 80 1000 400]);

subplot(1,2,1);
tChirp = (0:1/fs:tpd-1/fs);
Kchirp = bw / tpd;                       % chirp rate (Hz/s)
fInst  = Kchirp * tChirp;                % instantaneous frequency (baseband)
plot(tChirp*1e6, fInst/1e6, 'LineWidth', 2, 'Color', [0.17 0.43 0.39]);
grid on; box on;
xlabel('Time within pulse (\mus)');
ylabel('Instantaneous frequency (MHz)');
title(sprintf('LFM chirp: BW = %.1f MHz over T_p = %.1f \\mus', bw/1e6, tpd*1e6));
xlim([0 tpd*1e6]);
text(0.05*tpd*1e6, 0.85*bw/1e6, ...
    sprintf('K = B/T_p = %.2f MHz/\\mus\n\\delta r = c/2B = %.3f m', ...
            Kchirp/1e12, c/(2*bw)), ...
    'FontSize', 10, 'BackgroundColor', [0.95 0.97 0.98], 'EdgeColor', [0.7 0.75 0.8]);

subplot(1,2,2);
plot3(radarPosHistory(2,:), radarPosHistory(1,:), radarPosHistory(3,:), ...
      'LineWidth', 2.5, 'Color', [0.13 0.19 0.35]); hold on;
plot3(targetpos0(2,:), targetpos0(1,:), targetpos0(3,:), ...
      'r*', 'MarkerSize', 11, 'LineWidth', 2);
% เส้นชี้จากกึ่งกลางเส้นทางบินไปเป้าแต่ละตัว
midIdx = round(numpulses/2);
for k = 1:size(targetpos0,2)
    plot3([radarPosHistory(2,midIdx) targetpos0(2,k)], ...
          [radarPosHistory(1,midIdx) targetpos0(1,k)], ...
          [radarPosHistory(3,midIdx) targetpos0(3,k)], ...
          ':', 'Color', [0.4 0.6 0.75], 'LineWidth', 1);
end
hold off; grid on; box on;
xlabel('Cross-range x (m)'); ylabel('Ground range y (m)'); zlabel('Height z (m)');
title(sprintf('Geometry: h = %d m, v = %d m/s, aperture = %d m', ...
              round(radarPosHistory(3,1)), speed, round(speed*flightDuration)));
legend('Flight path','Point targets','Location','northeast');
view(-35, 20);

exportgraphics(fp1, fullfile(figDir,'fig_pipe1_radar_config.png'), 'Resolution', 150);
fprintf('\n[A] saved  fig_pipe1_radar_config.png\n');

% ---- stage 2-5: เซฟรูปที่มีอยู่แล้ว -----------------------------------
%   fig 1  = ground truth        -> stage 2  Platform / Scene Setup
%   fig 5  = raw echo + เส้นทฤษฎี -> stage 3  Raw Echo Simulation
%   fig 8  = ก่อน/หลัง RC         -> stage 4  Range Compression
%   fig 10 = BP image (dB)       -> stage 5  Backprojection

pipeMap = { 1, 'fig_pipe2_scene_setup.png',      'stage 2 — scene / ground truth'
            5, 'fig_pipe3_raw_echo.png',          'stage 3 — raw echo + theory curves'
            8, 'fig_pipe4_range_compression.png', 'stage 4 — before/after range compression'
           10, 'fig_pipe5_backprojection.png',    'stage 5 — focused BP image (dB)' };

for k = 1:size(pipeMap,1)
    fh = findobj('Type','figure','Number',pipeMap{k,1});
    if isempty(fh)
        fprintf('[A] !! ไม่พบ figure %d — ข้าม %s\n', pipeMap{k,1}, pipeMap{k,2});
        continue;
    end
    set(fh, 'Color', 'w');
    exportgraphics(fh, fullfile(figDir, pipeMap{k,2}), 'Resolution', 150);
    fprintf('[A] saved  %-38s (%s)\n', pipeMap{k,2}, pipeMap{k,3});
end

%% ====================================================================
%  PART B — วัด resolution บนกริดละเอียด (1-D cuts)
%  ====================================================================
%
%  ทำ BP ซ้ำเฉพาะเส้นตัดผ่านเป้า ไม่ใช่ทั้งภาพ -> เร็วมาก
%  ใช้ forward model เดียวกับใน W1-3_3point.m ทุกประการ

fprintf('\n[B] วัด resolution บนกริดละเอียด ...\n');

targets_range = targetpos0(1,:);          % [800 1000 1300]
Ltrack  = speed * flightDuration;         % aperture length (m)
lambda  = c / fc;
nSampRC = size(rxsigRC,1);
mfLen   = length(matchedFilter);

% --- helper: BP บนชุดพิกเซลใด ๆ (สูตรเดียวกับไฟล์หลัก) ---
bpAt = @(P) localBP(P, radarPosHistory, rxsigRC, numpulses, c, fs, fc, mfLen, nSampRC);

% --- helper: วัดความกว้างที่ระดับ frac ของ peak ด้วยการหาจุดตัดแบบ linear ---
widthAt = @(ax, prof, frac) localWidth(ax, prof, frac);

resR3 = zeros(1,3); resR6 = zeros(1,3);
resC3 = zeros(1,3); resC6 = zeros(1,3);
pkR   = zeros(1,3); pkC   = zeros(1,3);
thR   = c/(2*bw);
thC   = zeros(1,3);

for k = 1:3
    Rk = targets_range(k);

    % ---------- เส้นตัดแนว range (x = 0) ----------
    yFine = (Rk-12 : 0.01 : Rk+12);
    Pr    = [yFine; zeros(1,numel(yFine)); zeros(1,numel(yFine))];
    profR = abs(bpAt(Pr));
    profR = profR / max(profR);

    resR3(k) = widthAt(yFine, profR, 1/sqrt(2));   % -3 dB (half power)
    resR6(k) = widthAt(yFine, profR, 0.5);         % -6 dB (half amplitude)
    [~,iR]   = max(profR);  pkR(k) = yFine(iR);

    % ---------- เส้นตัดแนว cross-range (y = Rk) ----------
    xFine = (-1.5 : 0.002 : 1.5);
    Pc    = [Rk*ones(1,numel(xFine)); xFine; zeros(1,numel(xFine))];
    profC = abs(bpAt(Pc));
    profC = profC / max(profC);

    resC3(k) = widthAt(xFine, profC, 1/sqrt(2));
    resC6(k) = widthAt(xFine, profC, 0.5);
    [~,iC]   = max(profC);  pkC(k) = xFine(iC);

    thC(k) = lambda * Rk / (2 * Ltrack);

    fprintf('    T%d (R = %4d m) done\n', k, Rk);
end

%% ====================================================================
%  PART C — สรุปตัวเลข
%  ====================================================================

fprintf('\n');
fprintf('==============================================================\n');
fprintf(' EXACT EVALUATION — วัดบนกริดละเอียด (range 0.01 m, cross 0.002 m)\n');
fprintf('==============================================================\n');

fprintf('\nพารามิเตอร์ที่ใช้คำนวณทฤษฎี\n');
fprintf('  c        = %.6e m/s\n', c);
fprintf('  fc       = %.4f GHz   -> lambda = %.6f m\n', fc/1e9, lambda);
fprintf('  BW       = %.6f MHz\n', bw/1e6);
fprintf('  aperture = %.1f m  (v = %.0f m/s x T = %.1f s)\n', Ltrack, speed, flightDuration);

fprintf('\n--- RANGE RESOLUTION ---\n');
fprintf('  ทฤษฎี   delta_r = c/(2B) = %.6f m\n', thR);
fprintf('  %-4s %10s %12s %12s\n','','-3dB (มาตรฐาน)','-6dB','error vs ทฤษฎี');
for k = 1:3
    fprintf('  T%d   %10.4f m %12.4f m %12.4f m\n', ...
        k, resR3(k), resR6(k), resR3(k)-thR);
end

fprintf('\n--- CROSS-RANGE RESOLUTION ---\n');
fprintf('  ทฤษฎี   delta_x = lambda*R/(2L)\n');
fprintf('  %-4s %8s %14s %12s %12s\n','','R (m)','ทฤษฎี','-3dB (มาตรฐาน)','-6dB');
for k = 1:3
    fprintf('  T%d   %8.0f %14.4f m %12.4f m %12.4f m\n', ...
        k, targets_range(k), thC(k), resC3(k), resC6(k));
end

fprintf('\n--- PEAK LOCATION (บนกริดละเอียด) ---\n');
fprintf('  %-4s %14s %14s %12s %12s\n','','range จริง','range วัดได้','error','cross error');
for k = 1:3
    fprintf('  T%d   %14.0f %14.3f m %12.4f m %12.4f m\n', ...
        k, targets_range(k), pkR(k), pkR(k)-targets_range(k), pkC(k));
end

fprintf('\n--- เทียบกับที่โค้ดเดิมรายงาน ---\n');
fprintf('  โค้ดเดิมวัดที่ max/2 ของแอมพลิจูด = คอลัมน์ "-6dB" ข้างบน\n');
fprintf('  และวัดบนกริดหยาบ (range step %.2f m, cross step %.2f m)\n', ...
        yScene(2)-yScene(1), xScene(2)-xScene(1));
fprintf('  ซึ่งปัดค่าให้เป็นขั้นละ %.2f m ในแนว range\n', (yScene(2)-yScene(1))/2);
fprintf('==============================================================\n\n');

save(fullfile(figDir,'w1-3_exact_eval.mat'), ...
     'resR3','resR6','resC3','resC6','pkR','pkC','thR','thC', ...
     'targets_range','lambda','Ltrack','bw','fc','c');
fprintf('saved: figure/w1-3_exact_eval.mat\n');
fprintf('>>> ก๊อปบล็อก EXACT EVALUATION ข้างบนส่งกลับมาได้เลย <<<\n\n');

%% ====================================================================
%  LOCAL FUNCTIONS
%  ====================================================================

function img = localBP(pixelPos, radarPosHistory, rxsigRC, numpulses, ...
                       c, fs, fc, mfLen, nSampRC)
% BP บนชุดพิกเซลใด ๆ — สูตรเดียวกับลูปหลักใน W1-3_3point.m
    Npix = size(pixelPos,2);
    img  = zeros(1,Npix);
    for ii = 1:numpulses
        rp = radarPosHistory(:,ii);
        d  = pixelPos - rp;
        R  = sqrt(sum(d.^2,1));
        idx = 2.*R./c.*fs + mfLen;
        m   = (idx >= 1) & (idx <= nSampRC);
        v   = zeros(1,Npix);
        if any(m)
            v(m) = interp1(1:nSampRC, rxsigRC(:,ii), idx(m), 'linear', 0);
        end
        img = img + v .* exp(1j*4*pi*fc/c .* R);
    end
end

function w = localWidth(ax, prof, frac)
% ความกว้างของ mainlobe ที่ระดับ frac*peak
% หาจุดตัดด้วย linear interpolation ระหว่างตัวอย่าง -> ไม่ถูกกริดปัด
    [pk, ip] = max(prof);
    thr = frac * pk;

    iL = find(prof(1:ip) < thr, 1, 'last');
    if isempty(iL)
        xL = ax(1);
    else
        xL = interp1(prof(iL:iL+1), ax(iL:iL+1), thr);
    end

    iR = ip - 1 + find(prof(ip:end) < thr, 1, 'first');
    if isempty(iR)
        xR = ax(end);
    else
        xR = interp1(prof(iR-1:iR), ax(iR-1:iR), thr);
    end

    w = xR - xL;
end
