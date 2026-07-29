%% W10 ข้อ 4 — PER-PASS PHASE REFINEMENT ที่ตัวรถ (แก้ ghost แนวสูง)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%
%  #########################################################################
%  ##  ผลการทดลอง: สมมติฐานผิด — วิธีนี้ไม่ได้ช่วย และทำให้แย่ลงด้วยซ้ำ    ##
%  #########################################################################
%  รันจริง 18 look, HH, รถ 3 คัน:
%      ghost%  A(11m) 5% -> 5%   |  B(28m) 37% -> 33%  |  J(55m) 27% -> 44%
%      z p90   A 1.25 -> 1.25    |  B 1.88 -> 2.00     |  J 1.50 -> 1.75
%      sharpness +14% / +26% / +41%   <-- คมขึ้นจริง แต่ ghost แย่ลง
%
%  ตีความ: sharpness "ถูกโกงได้" — มันวัดแค่ว่าพลังงานกระจุกตัวไหม ไม่สนว่ากระจุกตรงไหน
%  พอ kz มีแค่ 8 ตัวอย่าง autofocus จึงจับพลังงานไปกระจุกที่ sidelobe ได้คะแนนดีกว่า
%  (guard ที่ใส่ไว้ช่วยไม่ได้ เพราะ guard ใช้ sharpness เป็นเกณฑ์ตัวเดียวกัน — ทำงานแค่ 1/54)
%
%  สาเหตุจริงของ ghost (คำนวณจากมุม elevation จริงของทั้ง 8 pass):
%      elevation ของ 8 pass = 44.07 ถึง 45.74 องศา = ช่วงแค่ 1.68 องศา
%      -> ความละเอียดแนวสูง (Rayleigh) = 0.75 m  ** รถสูง 1.43 m = แค่ 1.9 resolution cell **
%      -> sidelobe แนวสูงอยู่ที่ z = +-0.91 m (-8.9 dB) และ +-1.68 m (-9.2 dB)
%      เกณฑ์วัดของเราคือ -20 dB ซึ่ง "ต่ำกว่า" sidelobe เยอะ
%      => จุดสะท้อนจริงทุกจุดสร้างสำเนาปลอมของตัวเองที่สูงขึ้นไป 0.91 และ 1.68 m
%         ที่ระดับเพียง -9 dB แล้วถูกนับเป็นรถทั้งหมด
%
%  นี่คือปัญหา "สุ่มตัวอย่างไม่พอ" ไม่ใช่ "เฟสผิด" — เฟสคงที่ราย pass แก้ไม่ได้ในเชิงหลักการ
%  และเป็นเหตุผลว่าทำไม z p90 ที่เคยดู "ตรงหลังคา" ถึงเป็นเรื่องบังเอิญ:
%  มันคือ sidelobe -9 dB ของการสะท้อนจากพื้น/ตัวถัง ไม่ใช่หลังคา
%
%  เก็บไฟล์นี้ไว้เป็นหลักฐานของการทดสอบสมมติฐาน (ดู Week10 handoff)
%  ถ้าจะลองต่อ: เปลี่ยนเกณฑ์วัดเป็น -8 dB (สูงกว่า sidelobe) แทนการแก้เฟส
%
%  ===================== ปัญหาที่จะแก้ =====================
%  autofocus ที่ AFRL ให้มา (data.af) อ้างอิงที่ "ศูนย์กลางฉาก"
%  รถที่อยู่ไกลจากศูนย์กลางจึงยังเหลือ phase error ตกค้าง
%  แกนที่เจ็บที่สุดคือแกนสูง เพราะ aperture แนวสูงมาจาก 8 pass เท่านั้น
%  (แกน L/W ได้ aperture เต็มวง 360 องศา จึงแทบไม่กระทบ)
%
%  หลักฐานจาก W10 ข้อ 1:
%    - Δ ความสูง สัมพันธ์กับระยะจากศูนย์ฉาก r = +0.85 (+0.11 m ต่อ 10 m)
%      ขณะที่ W ได้ r = -0.34 และ L ได้ -0.25 คือแทบไม่มีผล
%    - โปรไฟล์ความสูงของรถใกล้ (Malibu 11 m) มีช่วงที่สูงเกินหลังคาแค่ 5%
%      แต่รถไกล (Prizm 55 m) สูงเกิน 24% | Camry 26% | SantaFe 35%
%
%  ===================== ขีดจำกัดที่ต้องรู้ก่อน =====================
%  ** เฟสคงที่ราย pass ที่ "เชิงเส้นใน kz" แยกไม่ออกจากความสูงจริงของเป้า **
%  เพราะเป้าที่สูง z ทำให้เกิด phase ramp ข้าม pass พอดี ๆ กับ error แบบนั้น
%  => autofocus แบบนี้ทำให้ภาพ "คม" ได้ แต่ยก/กดความสูงสัมบูรณ์ไม่ได้
%     สคริปต์จึงถอดองค์ประกอบเชิงเส้นใน kz ออกจากค่าแก้ทุกครั้ง (detrend)
%     ไม่งั้นมันจะเลื่อนรถขึ้นลงตามใจชอบ (sharpness ไม่สนตำแหน่ง z)
%  สิ่งที่คาดว่าจะดีขึ้น = การเบลอ/ghost | สิ่งที่แก้ไม่ได้ = offset ความสูงคงที่
%
%  ===================== วิธี =====================
%  iterative dominant-scatterer autofocus (sharpness maximization)
%    1. binning k-space ด้วยเฟสปัจจุบัน -> ภาพ (adjoint / 1 FFT)
%    2. เก็บเฉพาะ voxel สว่างสุด afKeepN จุด -> "แบบจำลองเป้าที่เชื่อว่าจริง"
%    3. แปลงกลับเป็น k-space แล้วหาเฟสที่ทำให้ข้อมูลจริงของแต่ละ pass ตรงกับแบบจำลอง
%    4. detrend เชิงเส้นใน kz -> วนซ้ำ
%  ต้นทุน ~2 FFT ต่อรอบ (ไม่ใช่ต่อ pass) จึงถูกกว่าการ search เฟสแบบ brute force มาก
%
%  ===================== การทดลอง =====================
%  รถ 3 คัน: A ใกล้สุด 11 m (ตัวควบคุม — ควรเปลี่ยนน้อย)
%            B กลาง 28 m | J ไกลสุด 55 m (ควรดีขึ้นชัดสุด)
%  รันคู่ before/after ด้วย engine เดียวกันเป๊ะ ต่างแค่ใส่เฟสหรือไม่ใส่
%  ** ตัวควบคุมสำคัญ: L กับ W ต้องแทบไม่เปลี่ยน ถ้าเปลี่ยนเยอะ = วิธีนี้ไปยุ่งผิดแกน **

clear; clc; close all;
c = physconst('LightSpeed');

%% ===== (0) CONFIG =====
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % root ของ repo
discRoot = fullfile(repoRoot, 'data');  % ที่วางข้อมูล GOTCHA
assert(exist(fullfile(discRoot,'GOTCHA-CP_Disc1','DATA'),'dir') > 0, ...
    'ไม่พบโฟลเดอร์ GOTCHA-CP_Disc1 ใน %s', discRoot);
figDir = fullfile(repoRoot, 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end

% ดึงตารางรถจาก cache ของ W9 (คอลัมน์: name id cx cy heading L W roofSpec)
cacheFile = fullfile(figDir, 'w9cmp2_cache.mat');
assert(exist(cacheFile,'file') > 0, 'ต้องมี %s ก่อน', cacheFile);
SC   = load(cacheFile, 'veh');
vehAll = SC.veh;

pickIDs = {'A','B','J'};               % ใกล้ / กลาง / ไกล
sel = ismember(vehAll(:,2), pickIDs);
veh = vehAll(sel, :);
Nveh = size(veh,1);

polList    = {'HH'};                   % 1 pol พอสำหรับการทดลอง (เร็วขึ้น 2 เท่า)
subApWidth = 5;
subApStep  = 20;                       % 18 subaperture (การทดลอง ไม่ใช่ production)
NsigCFAR   = 30;   relFloor = 0.02;    % ** เหมือน W9 เป๊ะ ห้ามเปลี่ยน **
nIterFISTA = 200;
gateRng    = 3.5;  gateCrs  = 4.5;
topDB      = 40;   metDB    = 20;
N  = 80;  dr = 0.125;                  % grid เดียวกับ W9_CompareVehicles
dk = 2*pi/(N*dr);
axv = ((0:N-1) - N/2) * dr;

% --- พารามิเตอร์ autofocus ---
% ** afKeep ตั้งเป็น "จำนวน voxel" ไม่ใช่สัดส่วน และต้องน้อยกว่าจำนวนจุดสะท้อนจริง **
%    เหตุผลจากการทดสอบด้วยฉากสังเคราะห์ที่รู้คำตอบ (ดูหัวข้อ VALIDATION ท้ายไฟล์):
%    ถ้าแบบจำลองเป้าใหญ่เกินจำนวนจุดจริง มันจะ "ดูดเอา error ไปไว้ในตัวเอง"
%    แล้วยืนยันเฟสเดิม -> estimator ไม่ขยับเลย (phi = 0 กลายเป็นจุดคงที่)
%    รถที่ -20 dB มี voxel ~150-275 จุด -> 40 จึงปลอดภัย
afIter     = 20;       % รอบการวนประมาณเฟส (ทดสอบแล้วลู่เข้าภายใน ~20)
afKeepN    = 40;       % จำนวน voxel สว่างสุดที่ใช้เป็นแบบจำลองเป้า
afMinPass  = 4;        % ต้องมีอย่างน้อยกี่ pass ถึงจะประมาณเฟส
afMaxCorr  = pi;       % จำกัดขนาดค่าแก้ (กันหลุด)
afGuard    = true;     % ถ้าแก้แล้วภาพ "คมลดลง" ให้ถอยกลับไปใช้เฟสศูนย์ (ไม่ทำร้ายของเดิม)

%% ===== (1) probe =====
passList = [];
for pn = 1:8
    dsc = 1 + (pn == 8);
    if exist(fullfile(discRoot, sprintf('GOTCHA-CP_Disc%d',dsc), 'DATA', ...
             sprintf('pass%d',pn), 'HH'), 'dir')
        passList(end+1) = pn; %#ok<SAGROW>
    end
end
nSub = numel(0:subApStep:359);
fprintf(['W10 #4 | per-pass phase refinement\n' ...
    'รถ %d คัน | passes %s | %d subaperture x %d pol | grid %d^3 @ %.3f m\n'], ...
    Nveh, mat2str(passList), nSub, numel(polList), N, dr);
for v = 1:Nveh
    fprintf('  %-16s %-3s | ระยะจากศูนย์ฉาก %.1f m\n', veh{v,1}, veh{v,2}, ...
        hypot(veh{v,3}, veh{v,4}));
end

for v = 1:Nveh
    V(v).ctr3 = [veh{v,3}, veh{v,4}, 0];
    V(v).raw  = zeros(N,N,N);          % before (เหมือน W9)
    V(v).ref  = zeros(N,N,N);          % after  (ใส่เฟสแก้แล้ว)
end
phiLog = cell(Nveh,1);                 % เก็บค่าแก้เฟสไว้ดูภายหลัง
shpLog = nan(Nveh, nSub*numel(polList), 2);

%% ===== (2) MAIN LOOP =====
nLook = 0;  nGuard = 0;  tAll = tic;
for aDeg = 0 : subApStep : 359
  for ip = 1:numel(polList)
    pol = polList{ip};
    passes = loadSubap(discRoot, passList, pol, aDeg, subApWidth);
    if numel(passes) < afMinPass, continue; end
    nLook = nLook + 1;
    for v = 1:Nveh
        % ---- k-samples แยกราย pass (ยังไม่ binning) ----
        G = kSamples(passes, V(v).ctr3, N, dk, gateRng, gateCrs, c);
        if numel(G) < afMinPass, continue; end

        % ---- BEFORE: เฟสศูนย์ = เหมือน W9 ทุกประการ ----
        [yg0, mk0] = binG(G, zeros(1,numel(G)), N);
        x0 = fistaKspace(yg0, mk0, nIterFISTA, NsigCFAR, relFloor);
        V(v).raw = max(V(v).raw, abs(x0));

        % ---- ประมาณเฟสรายเที่ยวบิน ----
        [phi, s0, s1] = estPassPhase(G, N, afIter, afKeepN, afMaxCorr);
        if afGuard && s1 < s0
            phi = zeros(1, numel(G));  s1 = s0;   % ถอยกลับ: ห้ามทำให้แย่ลง
            nGuard = nGuard + 1;
        end
        phiLog{v} = [phiLog{v}; phi];
        shpLog(v, nLook, :) = [s0, s1];

        % ---- AFTER ----
        [yg1, mk1] = binG(G, phi, N);
        x1 = fistaKspace(yg1, mk1, nIterFISTA, NsigCFAR, relFloor);
        V(v).ref = max(V(v).ref, abs(x1));
    end
    fprintf('  az %3d %s | %d passes | %.0f s (ETA %.0f min)\n', aDeg, pol, ...
        numel(passes), toc(tAll), toc(tAll)/nLook*(nSub*numel(polList)-nLook)/60);
  end
end
fprintf('รวม %d looks | total %.1f min\n', nLook, toc(tAll)/60);
if afGuard
    fprintf('guard ทำงาน %d ครั้ง จาก %d (look ที่แก้แล้วภาพคมลดลง จึงถอยกลับ)\n', ...
        nGuard, nLook*Nveh);
end

%% ===== (3) วัดผล before vs after =====
[XX, YY, ZZ] = ndgrid(axv, axv, axv);
vb  = -3 : 0.25 : 3;   vbc = vb(1:end-1) + 0.125;
res = struct();
for v = 1:Nveh
    for w = 1:2
        Vol = V(v).raw;  if w == 2, Vol = V(v).ref; end
        M = measureVol(Vol, XX, YY, ZZ, veh(v,:), topDB, metDB, vb);
        if w == 1, res(v).before = M; else, res(v).after = M; end
    end
end

fprintf('\n================ ผลเทียบ before / after ================\n');
fprintf('%-16s %-3s %7s | %12s %12s | %11s %11s\n', 'Vehicle','ID','ระยะ', ...
    'ghost% ก่อน','ghost% หลัง','z p90 ก่อน','z p90 หลัง');
fprintf('%s\n', repmat('-', 1, 86));
for v = 1:Nveh
    a = res(v).before;  b = res(v).after;
    fprintf('%-16s %-3s %6.0fm | %11.0f%% %11.0f%% | %11.2f %11.2f\n', ...
        veh{v,1}, veh{v,2}, hypot(veh{v,3},veh{v,4}), ...
        a.ghostPct, b.ghostPct, a.zp90, b.zp90);
end

fprintf('\n--- ตัวควบคุม: L / W ต้องแทบไม่เปลี่ยน (วิธีนี้ควรแตะแค่แกนสูง) ---\n');
fprintf('%-16s %-3s | %14s %14s | %14s %14s\n', 'Vehicle','ID', ...
    'L ก่อน/หลัง','ΔL','W ก่อน/หลัง','ΔW');
fprintf('%s\n', repmat('-', 1, 82));
for v = 1:Nveh
    a = res(v).before;  b = res(v).after;
    fprintf('%-16s %-3s | %6.2f %6.2f %+14.2f | %6.2f %6.2f %+14.2f\n', ...
        veh{v,1}, veh{v,2}, a.L, b.L, b.L-a.L, a.W, b.W, b.W-a.W);
end
fprintf(['ถ้า |ΔW| > ~0.05 m (คือใหญ่กว่าความแม่นของ W ที่วัดได้) = ผิดคาด\n' ...
    'แปลว่า autofocus ไปแก้อย่างอื่นนอกจากแกนสูง ต้องตรวจก่อนเชื่อผล\n']);

fprintf('\n--- ขนาดค่าแก้เฟสที่ประมาณได้ (rad) ---\n');
fprintf('%-16s %-3s %10s %10s %10s\n', 'Vehicle','ID','median|phi|','p90|phi|','max|phi|');
for v = 1:Nveh
    if isempty(phiLog{v}), continue; end
    ap = abs(phiLog{v}(:));
    fprintf('%-16s %-3s %10.3f %10.3f %10.3f\n', veh{v,1}, veh{v,2}, ...
        median(ap), prctile2(ap,90), max(ap));
end
fprintf(['ค่าแก้ควรอยู่ระดับ "เศษหนึ่งของ pi" ถ้าใหญ่เกือบ pi ทุกตัว\n' ...
    'แปลว่า estimator หลุด (หรือ SNR ต่ำเกินจะประมาณได้) อย่าเชื่อผล\n']);

sh = shpLog(:,1:nLook,:);
fprintf('\n--- ความคมของภาพ (sharpness) เฉลี่ยข้าม look ---\n');
for v = 1:Nveh
    s0 = mean(sh(v,:,1),'omitnan');  s1 = mean(sh(v,:,2),'omitnan');
    fprintf('  %-16s %-3s  %.4g -> %.4g  (%+.1f%%)\n', veh{v,1}, veh{v,2}, ...
        s0, s1, 100*(s1-s0)/s0);
end

%% ===== (4) FIGURES =====
f1 = figure(1); clf; set(f1,'Color','w','Position',[50 50 1150 380*Nveh], ...
    'Name','W10 phase refine — โปรไฟล์ความสูง');
tiledlayout(Nveh, 2, 'TileSpacing','compact','Padding','compact');
for v = 1:Nveh
    for w = 1:2
        M = res(v).before;  ttl = 'ก่อน (เหมือน W9)';
        if w == 2, M = res(v).after; ttl = 'หลังแก้เฟสรายเที่ยวบิน'; end
        nexttile; hold on;
        scatter(M.pv, M.pz, 6 + 30*max(0,(M.pm+topDB)/topDB).^1.5, M.pz, 'filled', ...
            'MarkerFaceAlpha', 0.7);
        stairs(vbc, M.roof, 'k-', 'LineWidth', 1.4);
        yline(veh{v,8}, '--', sprintf('หลังคา ~%.2f', veh{v,8}), 'Color',[.3 .3 .3]);
        yline(0, ':', 'Color', [.5 .5 .5]);
        hold off; grid on; colormap(turbo); clim([-0.3 2.4]);
        xlim([-3 3]); ylim([-1 3]);
        xlabel('v — ตามยาวรถ (m)'); ylabel('z (m)');
        title(sprintf('%s (%s) %s | ghost %.0f%%', veh{v,1}, veh{v,2}, ttl, M.ghostPct));
    end
end
sgtitle(sprintf('W10 #4 | per-pass phase refinement — %d looks, %s', nLook, ...
    strjoin(polList,'+')));
exportgraphics(f1, fullfile(figDir,'fig_w10_phaserefine_profiles.png'), 'Resolution', 150);

f2 = figure(2); clf; set(f2,'Color','w','Position',[60 60 900 400], ...
    'Name','W10 phase refine — สรุป');
gp = [arrayfun(@(v) res(v).before.ghostPct, 1:Nveh).', ...
      arrayfun(@(v) res(v).after.ghostPct,  1:Nveh).'];
bar(gp); grid on;
set(gca, 'XTick', 1:Nveh, 'XTickLabel', ...
    arrayfun(@(v) sprintf('%s (%.0f m)', veh{v,2}, hypot(veh{v,3},veh{v,4})), ...
    1:Nveh, 'UniformOutput', false));
ylabel('% ช่วงที่โปรไฟล์สูงเกินหลังคา'); legend({'ก่อน','หลัง'}, 'Location','best');
title('W10 #4 | ghost แนวสูง — ยิ่งเตี้ยยิ่งดี (รถไกลควรดีขึ้นมากกว่ารถใกล้)');
exportgraphics(f2, fullfile(figDir,'fig_w10_phaserefine_ghost.png'), 'Resolution', 150);

%% ===== (5) SAVE =====
save(fullfile(figDir,'w10_phaserefine_cache.mat'), 'V','res','veh','axv','nLook', ...
    'phiLog','shpLog','polList','subApStep','-v7.3');
fprintf('\nเซฟแล้ว: figure/w10_phaserefine_cache.mat + fig_w10_phaserefine_*.png\n');
fprintf(['\nอ่านผลยังไง:\n' ...
    '  1. ghost%% ของ J (ไกลสุด) ต้องลดลงชัด และ A (ใกล้สุด) ควรเปลี่ยนน้อย\n' ...
    '     -> ถ้าเป็นแบบนั้น = ยืนยันว่าสาเหตุคือ phase error ที่โตตามระยะ\n' ...
    '  2. ถ้า A แย่ลง = วิธีนี้ทำร้ายภาพที่ดีอยู่แล้ว ต้องลด afKeep หรือเพิ่ม afMinPass\n' ...
    '  3. ถ้า ghost ไม่ลดเลยทุกคัน = สมมติฐานผิด สาเหตุอยู่ที่อื่น\n' ...
    '     (เช่น ระยะห่าง elevation ที่ไม่สม่ำเสมอ ซึ่งแก้ด้วยเฟสคงที่ไม่ได้)\n']);

%% ===== LOCAL FUNCTIONS =====

function G = kSamples(passes, ctr3, N, dk, gateRng, gateCrs, c)
% เหมือน gridSubap ของ W9 ทุกอย่าง ต่างแค่ "ไม่ binning" — คืน index+ค่า แยกราย pass
% เพื่อให้ใส่เฟสรายเที่ยวบินก่อนรวมได้
    Q = numel(passes);
    KXc = cell(Q,1); KYc = cell(Q,1); KZc = cell(Q,1); VLc = cell(Q,1);
    for q = 1:Q
        fp = passes(q).fp;  ant = passes(q).ant;  r0u = passes(q).r0used;
        freq = passes(q).freq;
        Nf = numel(freq);  df = freq(2) - freq(1);  P = size(fp,2);
        Rc  = sqrt(sum((ant - ctr3).^2, 2)).';
        fpc = fp .* exp(1j*(4*pi/c)*freq*(Rc - r0u));
        Nz  = Nf*4;
        rc  = fftshift(ifft(fpc, Nz, 1), 1);
        drax = ((0:Nz-1).' - Nz/2) * (c/(2*df)) / Nz;
        rc(abs(drax) > gateRng, :) = 0;
        sp  = fft(ifftshift(rc,1), [], 1);
        fpc = sp(1:Nf, :);
        if P >= 64
            Sp = fft(fpc, [], 2);
            fd = (0:P-1)/P;  fd(fd >= 0.5) = fd(fd >= 0.5) - 1;
            azp = unwrap(atan2(ant(:,2), ant(:,1)));
            dphi = abs(azp(end) - azp(1)) / max(P-1, 1);
            uax = fd * 2*pi / ((4*pi*mean(freq)/c) * max(dphi, 1e-9));
            Sp(:, abs(uax) > gateCrs) = 0;
            fpc = ifft(Sp, [], 2);
        end
        U  = ant - ctr3;  U = U ./ sqrt(sum(U.^2, 2));
        kf = (4*pi/c) * freq;
        KXc{q} = reshape(kf*U(:,1).', [], 1);
        KYc{q} = reshape(kf*U(:,2).', [], 1);
        KZc{q} = reshape(kf*U(:,3).', [], 1);
        VLc{q} = fpc(:);
    end
    % K0 ต้องคิดจากทุก pass รวมกัน (เหมือน W9) แล้วค่อยแยก index ราย pass
    K0 = [mean(cell2mat(KXc)), mean(cell2mat(KYc)), mean(cell2mat(KZc))];
    G = struct('lin',{},'val',{},'kzm',{});
    for q = 1:Q
        ix = round((KXc{q}-K0(1))/dk) + N/2 + 1;
        iy = round((KYc{q}-K0(2))/dk) + N/2 + 1;
        iz = round((KZc{q}-K0(3))/dk) + N/2 + 1;
        ok = ix>=1 & ix<=N & iy>=1 & iy<=N & iz>=1 & iz<=N;
        if nnz(ok) < 100, continue; end
        G(end+1) = struct('lin', sub2ind([N N N], ix(ok), iy(ok), iz(ok)), ...
                          'val', VLc{q}(ok), ...
                          'kzm', mean(KZc{q}(ok)));  %#ok<AGROW>
    end
end

function [yg, mk] = binG(G, phi, N)
% รวมทุก pass ลง k-grid หลังคูณเฟสแก้ (เฉลี่ยราย bin เหมือน W9)
    Q = numel(G);
    lin = cell(Q,1);  val = cell(Q,1);
    for q = 1:Q
        lin{q} = G(q).lin;
        val{q} = G(q).val * exp(-1j*phi(q));
    end
    lin = cell2mat(lin);  val = cell2mat(val);
    cnt = accumarray(lin, 1,          [N^3 1]);
    yre = accumarray(lin, real(val),  [N^3 1]);
    yim = accumarray(lin, imag(val),  [N^3 1]);
    mk  = cnt > 0;
    yg  = zeros(N^3, 1);
    yg(mk) = (yre(mk) + 1j*yim(mk)) ./ cnt(mk);
    yg = reshape(yg, [N N N]);  mk = reshape(mk, [N N N]);
end

function [phi, sh0, sh1] = estPassPhase(G, N, nIt, keepN, maxCorr)
% autofocus แบบ dominant-scatterer: วนสร้างแบบจำลองเป้าจาก voxel สว่างสุด keepN จุด
% แล้วหาเฟสรายเที่ยวบินที่ทำให้ข้อมูลจริงของแต่ละ pass ตรงกับแบบจำลองนั้น
%
% ** keepN ต้องน้อยกว่าจำนวนจุดสะท้อนจริง ไม่งั้นล้มเหลวแบบเงียบ ๆ **
%    ทดสอบกับฉากสังเคราะห์: ฉาก 8 จุด ใช้ keepN=8 -> คลาด 0.0002 rad
%    แต่ keepN=32 -> คลาด 1.78 rad (เท่ากับสุ่มล้วน) เพราะแบบจำลองใหญ่พอจะ
%    อธิบายภาพที่ยังเบลออยู่ได้ จึงไม่มีแรงผลักให้เฟสขยับ (phi=0 เป็นจุดคงที่)
%
% ** detrend เชิงเส้นใน kz ทุกรอบ: องค์ประกอบนั้นคือการเลื่อนความสูงทั้งก้อน
%    ซึ่งแยกจากความสูงจริงไม่ได้ และ sharpness ก็ไม่สนใจ (อยู่ใน null space)
%    ถ้าไม่ถอดออก ค่าจะลอยไปเรื่อย ๆ แล้วยกรถขึ้นลงตามใจชอบ **
    Q   = numel(G);
    phi = zeros(1, Q);
    Nc  = sqrt(N^3);
    Aj  = @(y) fftshift( fftn(ifftshift(y))) / Nc;   % k -> image
    Fw  = @(x) fftshift(ifftn(ifftshift(x))) * Nc;   % image -> k
    kzm = [G.kzm];  kzm = kzm(:);
    A   = [ones(Q,1), kzm];
    sh0 = NaN;
    for it = 1:nIt
        x  = Aj(binG(G, phi, N));
        ax = abs(x);
        s  = sum(ax(:).^4) / max(sum(ax(:).^2)^2, eps);   % sharpness
        if it == 1, sh0 = s; end
        % แบบจำลองเป้า = voxel สว่างสุด keepN จุดเท่านั้น
        as  = sort(ax(:), 'descend');
        thr = as(min(keepN, numel(as)));
        xs  = x;  xs(ax < thr) = 0;
        if ~any(xs(:)), break; end
        Fx  = Fw(xs);
        dphi = zeros(1, Q);
        for q = 1:Q
            r = sum( (G(q).val * exp(-1j*phi(q))) .* conj(Fx(G(q).lin)) );
            if abs(r) > 0, dphi(q) = angle(r); end
        end
        phi = phi + dphi;
        if Q >= 3                                  % detrend a + b*kz
            cf  = A \ unwrap(phi(:));
            phi = phi(:).' - (A*cf).';
        end
        phi = angle(exp(1j*phi));                  % wrap เข้า [-pi,pi]
        phi(abs(phi) > maxCorr) = 0;
        if max(abs(dphi)) < 1e-4, break; end
    end
    ax  = abs(Aj(binG(G, phi, N)));
    sh1 = sum(ax(:).^4) / max(sum(ax(:).^2)^2, eps);
end

function M = measureVol(Vol, XX, YY, ZZ, vrow, topDB, metDB, vb)
% วัดแบบเดียวกับ W9/W10 เพื่อให้เทียบกันได้: point cloud -> L, W, z p90, โปรไฟล์, ghost%
    Mn     = Vol / max(Vol(:));
    selTop = Mn >= 10^(-topDB/20);
    nb     = convn(double(selTop), ones(3,3,3), 'same') - selTop;
    sel    = selTop & nb >= 1;
    thH = deg2rad(vrow{5});
    Rr  = [cos(thH) -sin(thH); sin(thH) cos(thH)];
    uv  = Rr.' * [XX(sel).'; YY(sel).'];
    M.pu = uv(1,:).';  M.pv = uv(2,:).';  M.pz = ZZ(sel);  M.pm = 20*log10(Mn(sel));
    GL = vrow{6};  GW = vrow{7};  GH = vrow{8};
    inEx = abs(M.pu) <= 3 & abs(M.pv) <= 3 & M.pm >= -metDB;
    M.L = NaN; M.W = NaN;
    if nnz(inEx) >= 10
        M.L = (prctile2(M.pv(inEx),95) - prctile2(M.pv(inEx),5)) / 0.90;  % p5-95/.90
        M.W =  prctile2(M.pu(inEx),90) - prctile2(M.pu(inEx),10);         % p10-90
    end
    inb = abs(M.pu) <= GW/2+0.15 & abs(M.pv) <= GL/2+0.15 & M.pz >= -0.5 & M.pz <= 3;
    mz  = inb & M.pm >= -metDB;
    M.zp90 = NaN;
    if nnz(mz) >= 10, M.zp90 = prctile2(M.pz(mz), 90); end
    m2 = M.pm >= -metDB & abs(M.pu) <= GW/2 + 0.2;
    roof = nan(1, numel(vb)-1);
    for j = 1:numel(vb)-1
        zz = M.pz(m2 & M.pv >= vb(j) & M.pv < vb(j+1));
        if ~isempty(zz), roof(j) = prctile2(zz, 90); end
    end
    M.roof = roof;
    vbc = vb(1:end-1) + 0.125;
    okp = isfinite(roof) & abs(vbc) <= GL/2;
    M.ghostPct = NaN;
    if any(okp), M.ghostPct = 100 * mean(roof(okp) > GH + 0.15); end
end

function y = prctile2(x, p)
% percentile แบบ prctile ของ MATLAB (ไม่ต้องพึ่ง Statistics Toolbox)
x = sort(x(:));  n = numel(x);
if n == 0, y = nan(size(p)); return; end
if n == 1, y = repmat(x, size(p)); return; end
q  = 100 * ((1:n).' - 0.5) / n;
pq = p(:);
y  = interp1(q, x, pq, 'linear');
y(pq <= q(1))   = x(1);
y(pq >= q(end)) = x(end);
y  = reshape(y, size(p));
end

function passes = loadSubap(discRoot, passList, pol, aDeg, nFiles)
% เหมือน W9 ทุกประการ (คัดลอกมาเพื่อให้ไฟล์นี้รันเดี่ยวได้)
    passes = struct('fp',{},'ant',{},'r0used',{},'freq',{});
    for pn = passList
        dsc = 1 + (pn == 8);
        fp=[]; ant=[]; r0u=[]; fq=[];
        for k = 1:nFiles
            azi = mod(aDeg + k - 1, 360) + 1;
            fn = fullfile(discRoot, sprintf('GOTCHA-CP_Disc%d',dsc), 'DATA', ...
                 sprintf('pass%d',pn), pol, ...
                 sprintf('data_3dsar_pass%d_az%03d_%s.mat', pn, azi, pol));
            if ~exist(fn,'file'), continue; end
            S = load(fn);  d = S.data;
            if isempty(fp)
                fq = double(d.freq(:));
            elseif size(d.fp,1) ~= numel(fq)
                continue;
            end
            phc = reshape(double(d.af.ph_correct), 1, []);
            rcr = reshape(double(d.af.r_correct),  1, []);
            fp  = [fp, double(d.fp).*exp(1j*phc)]; %#ok<AGROW>
            ant = [ant; [double(d.x(:)), double(d.y(:)), double(d.z(:))]]; %#ok<AGROW>
            r0u = [r0u, reshape(double(d.r0),1,[]) + rcr]; %#ok<AGROW>
        end
        if isempty(fp), continue; end
        passes(end+1) = struct('fp',fp,'ant',ant,'r0used',r0u,'freq',fq); %#ok<AGROW>
    end
end

function x = fistaKspace(yg, mk, nIter, Nsig, relFloor)
% เหมือน W9 ทุกประการ
    d  = yg .* mk;  Nc = sqrt(numel(d));
    Fw = @(x) fftshift(ifftn(ifftshift(x))) * Nc;
    Aj = @(y) fftshift( fftn(ifftshift(y))) / Nc;
    x0 = Aj(d);
    sig = median(abs(x0(:))) / 1.1774;
    thr = max(Nsig*sig, relFloor*max(abs(x0(:))));
    x = zeros(size(d));  z = x;  t = 1;
    for it = 1:nIter
        w  = z - Aj(Fw(z).*mk - d);
        aw = abs(w);
        xn = max(aw - thr, 0) .* (w ./ max(aw, 1e-30));
        tn = (1 + sqrt(1 + 4*t^2))/2;
        z  = xn + ((t-1)/tn)*(xn - x);
        if mod(it,10)==0 && norm(xn(:)-x(:)) < 1e-4*max(norm(x(:)),1e-12)
            x = xn; break;
        end
        x = xn;  t = tn;
    end
end
