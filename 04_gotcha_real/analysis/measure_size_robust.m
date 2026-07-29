%% W10 ข้อ 1 — ROBUST METRICS: วัด L / W / ความสูงหลังคา ใหม่จาก cache ของ W9
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%
%  ปัญหาที่แก้ (ดู Week9_Handoff.md):
%    W9 วัด L/W ด้วย min/max extent ของ voxel >= -20 dB
%    -> voxel หลงจุดเดียว (ghost / sidelobe) กำหนดคำตอบทั้งหมด
%    วิธีแก้: ใช้ order statistic (percentile) แทนค่าสุดขั้ว
%
%  สคริปต์นี้ **ไม่รัน recon ใหม่** — อ่าน figure/w9cmp2_cache.mat ใช้เวลาไม่กี่วินาที
%
%  =======================================================================
%  บทเรียนจากการรันรอบแรก (v1 ของสคริปต์นี้) — เหตุผลที่ต้องแก้เกณฑ์
%  -----------------------------------------------------------------------
%  v1 เลือก estimator ด้วย MAE (ค่าคลาดเฉลี่ย) ซึ่ง **ผิดสำหรับงาน discrimination**
%  เพราะรถ 6 คันในฉากยาวใกล้กันมาก (spec 4.41-4.98 = ช่วงแค่ 0.57 m)
%  estimator ที่ "ทายค่ากลางทุกคัน" จึงได้ MAE ต่ำ ทั้งที่ไม่มีข้อมูลแยกแยะเลย
%
%  หลักฐาน: minmax@-20dB (ที่ W9 ใช้) ได้ MAE L = 0.24 m ดูดี
%           แต่ correlation กับความยาวจริง r = -0.04, slope = -0.05
%           => มันไม่ได้ "วัดความยาว" มันแค่คืนเลขใกล้ ๆ ค่าเฉลี่ยของทุกคัน
%
%  ดังนั้น v2 เลือกด้วย 3 ตัวพร้อมกัน:
%    r      = correlation ระหว่างค่าที่วัดกับ spec  -> "ติดตามของจริงไหม" (สำคัญสุด)
%    slope  = ความชัน (1.0 = สเกลถูก, 0 = ทายค่ากลาง)
%    sd(Δ)  = ความกระจายของ error -> ความละเอียดที่แยกคันได้จริง
%  MAE ยังรายงานอยู่ แต่ **ห้ามใช้ตัวเดียวตัดสิน**
%  =======================================================================
%
%  หมายเหตุ "percentile แล้วหาร gain" (#6, #7): percentile ตัดขอบรถจริงไปด้วย
%  ถ้าจุดกระจายสม่ำเสมอ p5-95 ครอบ 90% ของช่วงจริง -> หารด้วย 0.90 คืนได้
%  แต่ถ้าหางเป็น ghost (ไม่ใช่ตัวรถ) การหารจะขยาย ghost -> ดูตัวชี้วัด r_u / r_v ตาราง 4
%
%  ** ตัวเลขทุกตัวต้องมาจากการรันจริง — เลขใน Week9_Handoff.md มาจาก simulation **

clear; clc; close all;

%% ===== (0) LOAD CACHE =====
root      = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % root ของ repo
cacheFile = fullfile(root, 'figure', 'w9cmp2_cache.mat');
assert(exist(cacheFile, 'file') > 0, ...
    ['ไม่พบ %s\n' ...
     'ต้องรัน W9_CompareVehicles.m ให้จบก่อน (มันเซฟ cache ตอนท้าย)'], cacheFile);

S      = load(cacheFile);             % V, R, axv, veh, nLook
R      = S.R;   veh = S.veh;   nLook = S.nLook;
Nveh   = size(veh, 1);
dVox   = S.axv(2) - S.axv(1);         % ขนาด voxel (m) — เพดานความละเอียด
figDir = fullfile(root, 'figure');
if ~exist(figDir, 'dir'); mkdir(figDir); end

% ---- ที่มาของ spec (ตรวจแล้ว ตุลา W10) ----
% L, W : มาจาก AFRL โดยตรง — "Gotcha Spotlight Target Locations.xls" ชีต
%        "Vehicle Dimensions" | ตรวจไขว้กับพิกัดมุมรถ 4 จุด (ชีต Vehicles) แล้วตรงกัน
%        => เป็น ground truth จริง เอาไปคำนวณ correlation ได้
% H    : **ไม่มีในชุดข้อมูล** — ชีต Dimensions มีแค่ L กับ W ส่วน z ของมุมรถ ~ 0 m
%        (จุดที่ล้อแตะพื้น ไม่ใช่หลังคา) ค่า roof ในตารางนี้เราหามาจากภายนอก
%        => ใช้เป็น "การตรวจความสมเหตุสมผล" เท่านั้น ห้ามอ้างเป็น ground truth
GT_L = cell2mat(veh(:, 6));           % ยาว (m)  — AFRL xls
GT_W = cell2mat(veh(:, 7));           % กว้าง (m) — AFRL xls
GT_H = cell2mat(veh(:, 8));           % สูงหลังคา (m) — เราหาเอง ไม่ใช่ GT

% ---- ตรวจ roofSpec กับสเปกผู้ผลิตบนเว็บ (W10) ----
% ข้อมูลเก็บปี 2006 -> รถเป็นรุ่นปี ~1998-2006 | ค่าที่ค้นได้ (นิ้ว -> m):
%   A Malibu  57.5" = 1.46 (รุ่น 2004-07) | 56.4" = 1.43 (รุ่น 1997-2003)
%   B Camry   58.8" = 1.50 (XV30 2002-06) | 55.4" = 1.41 (XV20 1997-01)  << กำกวม
%   C TaurusW 58.0" = 1.47 (2000)
%   D Maxima  56.5" = 1.44 (2000)
%   E Sentra  55.5" = 1.41 (2002)
%   F SantaFe 66.0" = 1.68 (2003)
%   J Prizm   53.7" = 1.36 (2002)
% สรุป: ตาราง roofSpec เดิมตรงกับเว็บภายใน 0.02 m ทุกคัน ยกเว้น Camry ที่แยกรุ่นปีไม่ได้
% (AFRL ให้ L 4.75 ซึ่งสั้นกว่าทั้งสองรุ่น ~0.06 m จึงใช้ชี้รุ่นไม่ได้)
% -> ความไม่แน่นอนของ "ค่าอ้างอิง" เอง = 0.088 m สำหรับ Camry ซึ่งเท่ากับ 80%
%    ของช่วงความสูงจริงทั้งกลุ่ม (1.36-1.47 = 0.11 m) => ต่อให้เรดาร์วัดเป๊ะก็ตรวจสอบไม่ได้
GT_Hunc = zeros(Nveh,1);              % ความไม่แน่นอนของค่าอ้างอิงความสูง (m)
iCam = find(strcmp(veh(:,2),'B'), 1);
if ~isempty(iCam), GT_Hunc(iCam) = 0.088; end   % Camry: XV20 vs XV30
ids  = veh(:, 2).';

% G (SaturnIon) / H (VWJetta) = ไม่ใช่รถ (finding W9: ตำแหน่งนั้นว่างตอนเก็บข้อมูลจริง)
% ถ้า cache เป็นชุด 11 คันเดิม ต้องกันออกจากค่าเฉลี่ย ไม่งั้น clutter ปนสถิติ
isVeh = ~ismember(veh(:, 2), {'G','H'});
% "เก๋ง" = กลุ่มรูปทรงเดียวกัน ใช้สรุปหลัก (tractor/forklift/SUV รูปทรงต่างชั้น)
isCar = ismember(veh(:, 2), {'A','B','C','D','E','J'});

%% ===== (1) นิยาม ESTIMATOR =====
EST = struct( ...
    'name', {'minmax@-20dB', 'p5-95@-20dB', 'p10-90@-20dB', 'minmax@-10dB', ...
             'p5-95@-10dB',  'p5-95/.90',   'p10-90/.80'}, ...
    'dB',   {20,  20,  20,  10,  10,  20,  20}, ...
    'plo',  { 0,   5,  10,   0,   5,   5,  10}, ...
    'phi',  {100, 95,  90, 100,  95,  95,  90}, ...
    'gain', { 1,   1,   1,   1,   1, 0.90, 0.80});
nE = numel(EST);
iBase = 1;                            % #1 = ที่ W9 ใช้อยู่

ZEST = struct('name', {'z median', 'z p75', 'z p90', 'z p95', 'z max'}, ...
              'p',    {50, 75, 90, 95, 100});
nZ = numel(ZEST);

exR   = 3.0;  metDB = 20;  zLo = -0.5;  zHi = 3.0;  minPt = 10;

%% ===== (2) คำนวณทุก ESTIMATOR =====
Lm = nan(Nveh, nE);  Wm = nan(Nveh, nE);  nUse  = zeros(Nveh, nE);
Zm = nan(Nveh, nZ);                       nZuse = zeros(Nveh, 1);
outFracU = nan(Nveh, 1);  tailU = nan(Nveh, 1);
shapeRu  = nan(Nveh, 1);  shapeRv = nan(Nveh, 1);
prec     = nan(Nveh, 1);

vb    = -3 : 0.25 : 3;
roofP = nan(Nveh, numel(vb) - 1);
roofM = nan(Nveh, numel(vb) - 1);
roofBad = nan(Nveh, 1);   % % ช่องโปรไฟล์ที่สูงเกินหลังคา = ghost ปน
roofCab = nan(Nveh, 1);   % ความสูงห้องโดยสาร (จากโปรไฟล์ ไม่ใช่กองรวม)
vbc     = vb(1:end-1) + 0.125;

for v = 1:Nveh
    pu = R(v).pu;  pv = R(v).pv;  pz = R(v).pz;  pm = R(v).pm;
    if isempty(pu), continue; end
    if isfield(R(v), 'prec'), prec(v) = R(v).prec; end

    inEx = abs(pu) <= exR & abs(pv) <= exR;

    for e = 1:nE
        m = inEx & pm >= -EST(e).dB;
        nUse(v, e) = nnz(m);
        if nnz(m) >= minPt
            Lm(v,e) = (pctl(pv(m), EST(e).phi) - pctl(pv(m), EST(e).plo)) / EST(e).gain;
            Wm(v,e) = (pctl(pu(m), EST(e).phi) - pctl(pu(m), EST(e).plo)) / EST(e).gain;
        end
    end

    m20 = inEx & pm >= -metDB;
    if nnz(m20) >= minPt
        au = abs(pu(m20));
        outFracU(v) = 100 * mean(au > GT_W(v)/2);
        tailU(v)    = max(au) - pctl(au, 95);
        shapeRu(v)  = spanRatio(pu(m20));      % แกนกว้าง
        shapeRv(v)  = spanRatio(pv(m20));      % แกนยาว
    end

    inb = abs(pu) <= GT_W(v)/2 + 0.15 & abs(pv) <= GT_L(v)/2 + 0.15 & ...
          pz >= zLo & pz <= zHi;
    mz = inb & pm >= -metDB;
    nZuse(v) = nnz(mz);
    if nnz(mz) >= minPt
        for k = 1:nZ, Zm(v,k) = pctl(pz(mz), ZEST(k).p); end
    end

    m2 = pm >= -metDB & abs(pu) <= GT_W(v)/2 + 0.2;
    for j = 1:numel(vb) - 1
        zz = pz(m2 & pv >= vb(j) & pv < vb(j+1));
        if ~isempty(zz)
            roofM(v,j) = max(zz);
            roofP(v,j) = pctl(zz, 90);
        end
    end
    % ---- ความสูงที่ "รู้ว่ากำลังวัดส่วนไหนของรถ" ----
    % Zm ข้างบนกอง voxel ทั้งคันรวมกัน (ฝากระโปรง+หลังคา+ท้าย+ล้อ+ghost) แล้วหา percentile
    % => ไม่ใช่ความสูงหลังคา เป็นแค่สถิติของกองความสูงทั้งคัน
    % ส่วนโปรไฟล์ roofP รู้ตำแหน่งตามยาว จึงแยก "ห้องโดยสาร" ออกมาได้จริง
    onCar = abs(vbc) <= GT_L(v)/2;
    pr    = roofP(v, :);
    okp   = onCar & isfinite(pr);
    if any(okp)
        % contamination: ช่องที่สูงเกินหลังคาที่คาด = ghost (รถจริงสูงเกินตัวเองไม่ได้)
        roofBad(v) = 100 * mean(pr(okp) > GT_H(v) + 0.15);
        % roof จากห้องโดยสาร: กลางรถครึ่งหนึ่งของความยาว, p75 กันช่อง ghost เดี่ยว ๆ
        cab = okp & abs(vbc) <= GT_L(v)/4;
        if nnz(cab) >= 3, roofCab(v) = pctl(pr(cab), 75); end
    end
end

dL = Lm - GT_L;   dW = Wm - GT_W;   dH = Zm - GT_H;

%% ===== (3) SELF-CHECK =====
fprintf('===== SELF-CHECK: คำนวณ baseline ของ W9 ใหม่จาก cache =====\n');
okChk = true;  worst = 0;
for v = 1:Nveh
    if isfield(R(v), 'L') && ~isempty(R(v).L) && R(v).L > 0
        eL = abs(Lm(v,iBase) - R(v).L);  eW = abs(Wm(v,iBase) - R(v).W);
        worst = max([worst, eL, eW]);
        if eL > 1e-6 || eW > 1e-6, okChk = false; end
    end
end
if okChk
    fprintf('  OK — minmax@-20dB ตรงกับ cache ทุกคัน (ต่างสุด %.2e m)\n', worst);
else
    warning('baseline ไม่ตรงกับ cache — ตรวจ mask ก่อนใช้ผล');
end
fprintf('  voxel = %.3f m | cache มี %d แถว (นับเป็นรถจริง %d, เก๋ง %d)\n', ...
    dVox, Nveh, nnz(isVeh), nnz(isCar));

%% ===== (4)(5)(6) ตารางเทียบ =====
printTable('ตาราง 1 | ความยาว L', Lm, GT_L, ids, {EST.name}, isCar, isVeh, nLook, iBase);
printTable('ตาราง 2 | ความกว้าง W', Wm, GT_W, ids, {EST.name}, isCar, isVeh, nLook, iBase);
printTable('ตาราง 3 | ความสูงหลังคา', Zm, GT_H, ids, {ZEST.name}, isCar, isVeh, nLook, 1);
fprintf('หมายเหตุ: roof spec เป็นค่าอ้างอิงภายนอก (~) ไม่ใช่ ground truth ที่วัดหน้างาน\n');

%% ===== (7) เพดานความละเอียด — เช็คว่า "แยกได้ในหลักการ" หรือเปล่า =====
fprintf('\n===== ตาราง 4 | เพดานความละเอียด (voxel = %.3f m) =====\n', dVox);
fprintf('%-10s %10s %10s %12s\n', 'มิติ', 'ช่วง spec', 'เป็นกี่ voxel', 'แยกได้?');
fprintf('%s\n', repmat('-', 1, 46));
dimN = {'L', 'W', 'H'};  dimG = {GT_L, GT_W, GT_H};
for q = 1:3
    g  = dimG{q}(isCar);  sp = max(g) - min(g);
    fprintf('%-10s %9.2f m %10.1f %12s\n', dimN{q}, sp, sp/dVox, ...
        ternary(sp < 2*dVox, 'ไม่ได้', 'ได้'));
end
fprintf(['ถ้าช่วง spec ของกลุ่ม < ~2 voxel = แยกคันในกลุ่มไม่ได้ในเชิงหลักการ\n' ...
    'ไม่ว่าจะใช้ estimator อะไร (ข้อจำกัดของ grid ไม่ใช่ของอัลกอริทึม)\n']);

%% ===== (8) กลไก =====
fprintf('\n===== ตาราง 5 | หลักฐานเชิงกลไก (ชุด -20 dB) =====\n');
fprintf('%-16s %-3s %7s %9s %8s %9s %7s %7s\n', 'Vehicle','ID','nVox', ...
    'outside%','tail(m)','W minmax','r_u','r_v');
fprintf('%s\n', repmat('-', 1, 74));
for v = 1:Nveh
    fprintf('%-16s %-3s %7d %8.1f%% %8.2f %9.2f %7.3f %7.3f\n', veh{v,1}, veh{v,2}, ...
        nUse(v,iBase), outFracU(v), tailU(v), Wm(v,iBase), shapeRu(v), shapeRv(v));
end
fprintf([ ...
    'outside%% = สัดส่วน voxel นอกความกว้างจริง | tail = ระยะ p95->max ของ |u|\n' ...
    'r_u, r_v = span(p5-95)/span(p10-90) ต่อแกน | จุดสม่ำเสมอ -> 1.125\n' ...
    '  ต่ำกว่า 1.125 = มวลเกาะขอบ (หางบาง) | สูงกว่า = มวลเกาะกลาง (หางหนา = ghost เยอะ)\n' ...
    'แกนที่ r สูงกว่า = แกนที่ ghost แผ่มากกว่า -> ต้อง trim หนักกว่า\n']);
fprintf('-> มัธยฐานเก๋ง: r_u = %.3f, r_v = %.3f\n', ...
    median(shapeRu(isCar),'omitnan'), median(shapeRv(isCar),'omitnan'));

%% ===== (9) เลือก estimator แยกรายแกน =====
% เหตุผลที่ต้องแยก: ghost แผ่ไม่เท่ากันในแต่ละแกน (ดู bias ของ minmax ข้างล่าง)
% การบังคับใช้ estimator ตัวเดียวทั้ง L และ W = ยอมให้แกนหนึ่งพังเพื่ออีกแกน
fprintf('\n===== การเลือก estimator (แยกรายแกน) =====\n');
biasLmm = mean(dL(isCar,iBase),'omitnan');  biasWmm = mean(dW(isCar,iBase),'omitnan');
fprintf('หลักฐานว่า ghost แผ่ไม่สมมาตร: bias ของ minmax@-20dB\n');
fprintf('  แกนยาว (L) %+.2f m | แกนกว้าง (W) %+.2f m  -> ต่างกัน %.1f เท่า\n', ...
    biasLmm, biasWmm, abs(biasWmm/biasLmm));
fprintf('  => ghost แผ่ในแนว cross-range (กว้าง) เป็นหลัก -> W ต้อง trim หนักกว่า L\n\n');

[rL, slL, sdL] = deal(nan(1,nE), nan(1,nE), nan(1,nE));
[rW, slW, sdW] = deal(nan(1,nE), nan(1,nE), nan(1,nE));
for e = 1:nE
    [rL(e), slL(e)] = corrSlope(GT_L(isCar), Lm(isCar,e));  sdL(e) = std(dL(isCar,e),'omitnan');
    [rW(e), slW(e)] = corrSlope(GT_W(isCar), Wm(isCar,e));  sdW(e) = std(dW(isCar,e),'omitnan');
end
[rH, slH, sdH] = deal(nan(1,nZ), nan(1,nZ), nan(1,nZ));
for k = 1:nZ
    [rH(k), slH(k)] = corrSlope(GT_H(isCar), Zm(isCar,k));  sdH(k) = std(dH(isCar,k),'omitnan');
end

% เกณฑ์: ต้อง "ติดตามของจริง" ก่อน (r) แล้วค่อยดูความกระจาย (sd)
% ถ้า r ต่ำหมดทุกตัว = มิตินั้นวัดไม่ได้ ไม่ว่าเลือกอะไร -> ต้องบอกตรง ๆ
maeLv = mean(abs(dL(isCar,:)), 1, 'omitnan');
maeWv = mean(abs(dW(isCar,:)), 1, 'omitnan');
maeHv = mean(abs(dH(isCar,:)), 1, 'omitnan');
iBL = pickBest(rL, sdL, maeLv);
iBW = pickBest(rW, sdW, maeWv);
iBz = pickBest(rH, sdH, maeHv);

fprintf('%-14s %7s %7s %7s %7s  <- L\n', 'estimator', 'r', 'slope', 'sd', 'MAE');
for e = 1:nE
    fprintf('%s%-13s %7.3f %7.2f %7.3f %7.2f\n', markOf(e,iBL,iBase), EST(e).name, ...
        rL(e), slL(e), sdL(e), mean(abs(dL(isCar,e)),'omitnan'));
end
fprintf('\n%-14s %7s %7s %7s %7s  <- W\n', 'estimator', 'r', 'slope', 'sd', 'MAE');
for e = 1:nE
    fprintf('%s%-13s %7.3f %7.2f %7.3f %7.2f\n', markOf(e,iBW,iBase), EST(e).name, ...
        rW(e), slW(e), sdW(e), mean(abs(dW(isCar,e)),'omitnan'));
end
fprintf('\n%-14s %7s %7s %7s %7s  <- ความสูง\n', 'estimator', 'r', 'slope', 'sd', 'MAE');
for k = 1:nZ
    fprintf('%s%-13s %7.3f %7.2f %7.3f %7.2f\n', markOf(k,iBz,1), ZEST(k).name, ...
        rH(k), slH(k), sdH(k), mean(abs(dH(isCar,k)),'omitnan'));
end
fprintf('(-> = ชนะ, * = ที่ W9 ใช้อยู่ | r ใกล้ 1 + slope ใกล้ 1 = วัดของจริง)\n');

%% ===== (10) ข้อสรุปสำหรับรายงาน / สไลด์ 11 =====
fprintf('\n===== ข้อสรุปที่ต้องเอาไปแก้สไลด์ 11 (Feature Ranking) =====\n');
% เทียบกับ NULL MODEL = "ทายค่าเฉลี่ยของกลุ่มให้ทุกคันเท่ากัน"
% ถ้า estimator ชนะ null ไม่ได้ = ไม่มี skill เลย ต่อให้ MAE ดูต่ำแค่ไหน
% varExp = 1 - (sd ของ residual / sd ของ spec)^2  (หัก bias คงที่ออกก่อน เพราะ calibrate ได้)
dims = {'L', GT_L, Lm(:,iBL), dL(:,iBL), EST(iBL).name;
        'W', GT_W, Wm(:,iBW), dW(:,iBW), EST(iBW).name;
        'H', GT_H, Zm(:,iBz), dH(:,iBz), ZEST(iBz).name};
fprintf('%-3s %-14s %8s %9s %9s %9s  %s\n', ...
    'มิติ','estimator','sd(spec)','sd(resid)','อธิบาย%','MAEnull','คำตัดสิน');
fprintf('%s\n', repmat('-', 1, 96));
verdicts = cell(3,1);
for q = 1:3
    [nm, gt, ~, dd, en] = dims{q,:};
    g   = gt(isCar);  e = dd(isCar);
    sp  = max(g) - min(g);
    sds = std(g, 'omitnan');  sdr = std(e, 'omitnan');
    varExp  = 1 - (sdr/sds)^2;
    maeEst  = mean(abs(e), 'omitnan');
    maeNull = mean(abs(g - mean(g,'omitnan')), 'omitnan');
    if sp < 2*dVox
        verdict = 'แยกในกลุ่มไม่ได้ — เพดาน grid';
    elseif varExp >= 0.7 && maeEst < maeNull
        verdict = 'วัดได้จริง';
    elseif varExp > 0 && maeEst < maeNull
        verdict = 'พอมีข้อมูล แต่อ่อน';
    else
        verdict = 'ไม่มี skill — แย่กว่าเดาค่าเฉลี่ย';
    end
    verdicts{q} = verdict;
    fprintf('%-3s %-14s %8.3f %9.3f %8.1f%% %9.3f  %s\n', ...
        nm, en, sds, sdr, 100*varExp, maeNull, verdict);
end
fprintf(['\nอ่านตาราง: "อธิบาย%%" คือสัดส่วนความแปรปรวนของ spec ที่ estimator อธิบายได้\n' ...
    '  ติดลบ = แย่กว่าทายค่าเฉลี่ยกลุ่ม | MAEnull = MAE ของการทายค่าเฉลี่ย\n' ...
    '  bias คงที่ไม่นับเป็น error (แก้ด้วยการ calibrate ได้) จึงดู sd(resid) ไม่ใช่ MAE\n']);

% ---- leave-one-out: กันข้อสรุปถูกลากด้วยรถคันเดียว ----
fprintf('\n--- leave-one-out ของตัวชนะ (ตัดทีละคันแล้วดู r) ---\n');
for q = 1:3
    [nm, gt, ft, ~, en] = dims{q,:};
    g = gt(isCar);  f = ft(isCar);  idc = ids(isCar);
    rs = nan(1, numel(g));
    for i = 1:numel(g)
        k = true(size(g));  k(i) = false;
        rs(i) = corrSlope(g(k), f(k));
    end
    fprintf('  %s [%-13s] r เต็ม %+.3f | ตัดทีละคัน %+.3f ถึง %+.3f', ...
        nm, en, corrSlope(g, f), min(rs), max(rs));
    [~, iw] = min(rs);
    fprintf(' (อ่อนสุดตอนตัด %s)\n', idc{iw});
end
fprintf(['ถ้าช่วง r ตอนตัดทีละคันยังสูงใกล้เคียงกันหมด = ไม่ได้พึ่งรถคันใดคันหนึ่ง\n' ...
    'ถ้าร่วงแรงตอนตัดคันใดคันหนึ่ง = correlation มาจากจุดนั้นจุดเดียว อย่าเชื่อ\n']);
fprintf('\nตัวเลขทั้งหมดมาจากเก๋ง %d คัน — กลุ่มเล็ก ดูผลตรวจอิสระข้างล่างประกอบ\n', nnz(isCar));

%% ===== (10a) ความสูง: ตรวจว่าเรากำลังวัด "ส่วนไหนของรถ" =====
% รถไม่ได้สูงเท่ากันตลอดคัน (ฝากระโปรง ~1.0 / หลังคา ~1.45 / ท้าย ~1.2)
% Zm ที่ใช้ในตาราง 3 กอง voxel ทั้งคันรวมกันแล้วหา percentile -> ไม่รู้ว่าส่วนไหน
% ตารางนี้ใช้โปรไฟล์ตามยาว (roofP) ซึ่งรู้ตำแหน่ง จึงบอกได้ว่าโปรไฟล์เป็นรูปรถจริงไหม
fprintf('\n===== ตาราง 5a | ความสูง: โปรไฟล์เป็นรูปทรงรถหรือเป็น ghost =====\n');
fprintf('%-16s %-3s %8s %11s %10s %9s %8s\n', 'Vehicle','ID','ระยะ(m)', ...
    'สูงเกิน%','ห้องโดยสาร','Zm p90','spec~');
fprintf('%s\n', repmat('-', 1, 70));
for v = 1:Nveh
    if ~isVeh(v), continue; end
    fprintf('%-16s %-3s %8.0f %10.0f%% %10.2f %9.2f %8.2f\n', veh{v,1}, veh{v,2}, ...
        hypot(veh{v,3}, veh{v,4}), roofBad(v), roofCab(v), Zm(v,iBz), GT_H(v));
end
fprintf([ ...
    '"สูงเกิน%%" = สัดส่วนช่วงตามยาวที่โปรไฟล์สูงเกินหลังคาที่คาด >0.15 m\n' ...
    '  รถจริงสูงเกินตัวเองไม่ได้ -> ตัวเลขนี้คือ ghost ล้วน ๆ | ต่ำ = โปรไฟล์เชื่อได้\n' ...
    '"ห้องโดยสาร" = p75 ของโปรไฟล์เฉพาะกลางรถ (|v| <= L/4) = ความสูงหลังคาที่แท้จริง\n' ...
    '  ถ้าต่างจาก Zm p90 มาก แปลว่า Zm ไม่ได้กำลังวัดหลังคา แค่บังเอิญได้เลขใกล้กัน\n']);

%% ===== (10b) ผลของระยะห่างจากจุดโฟกัส (เจอตอนไล่หาสาเหตุที่ J สูงผิดปกติ) =====
% autofocus/เฟสอ้างอิงทำที่ tophat กลางฉาก -> phase error ที่เหลือโตตามระยะ
% แกนที่เปราะสุดคือ "ความสูง" เพราะ aperture แนวสูงมาจาก 8 pass เท่านั้น
% ส่วน L/W ได้ aperture เต็มวง 360 องศา จึงควรทนกว่ามาก -> ทดสอบสมมติฐานนี้
rngC = hypot(cell2mat(veh(:,3)), cell2mat(veh(:,4)));   % ระยะจากศูนย์ฉาก
fprintf('\n===== ตาราง 5b | ระยะจากจุดโฟกัส มีผลกับมิติไหน =====\n');
fprintf('%-4s %14s %14s %10s\n', 'มิติ', 'corr(ระยะ,Δ)', 'slope/10 m', 'sd(Δ)');
fprintf('%s\n', repmat('-', 1, 46));
rngDim = {'L', dL(:,iBL); 'W', dW(:,iBW); 'H', dH(:,iBz)};
slopeH = NaN;
for q = 1:3
    dd = rngDim{q,2};
    ok = isCar & isfinite(dd) & isfinite(rngC);
    if nnz(ok) >= 3
        cc = corrcoef(rngC(ok), dd(ok));
        pp = polyfit(rngC(ok), dd(ok), 1);
        if q == 3, slopeH = pp(1); end
        fprintf('%-4s %14.3f %14.3f %10.3f\n', rngDim{q,1}, cc(1,2), pp(1)*10, ...
            std(dd(ok),'omitnan'));
    end
end
fprintf(['ถ้าเฉพาะ H ที่ corr สูง = ยืนยันว่าเป็นปัญหาแกนสูงโดยเฉพาะ (aperture แคบสุด)\n' ...
    'ไม่ใช่ปัญหาของ estimator หรือของรถคันใดคันหนึ่ง\n']);
if ~isnan(slopeH)
    okH = isCar & isfinite(dH(:,iBz));
    pH  = polyfit(rngC(okH), dH(okH,iBz), 1);
    res = dH(okH,iBz) - polyval(pH, rngC(okH));
    fprintf('ถ้าแก้ bias เชิงระยะของ H: sd(Δ) %.3f -> %.3f m\n', ...
        std(dH(okH,iBz),'omitnan'), std(res,'omitnan'));
    fprintf(['** fit 2 พารามิเตอร์บนรถ %d คัน = สมมติฐาน ไม่ใช่ค่าแก้ที่ยืนยันแล้ว **\n'], nnz(okH));
    fprintf([ ...
        '** อัปเดต: ทดสอบแล้วด้วย W10_PhaseRefine.m — สาเหตุ "phase error" ผิด **\n' ...
        '   ทำ per-pass phase refinement แล้ว ghost แย่ลง (J: 27%% -> 44%%)\n' ...
        '   สาเหตุจริงคือ elevation aperture แคบ: 8 pass กินช่วงแค่ 1.68 องศา\n' ...
        '   -> ความละเอียดแนวสูง 0.75 m และ sidelobe อยู่ที่ +-0.91 m ระดับ -8.9 dB\n' ...
        '   ซึ่งสูงกว่าเกณฑ์วัด -20 dB มาก => ทุกจุดสะท้อนสร้างสำเนาปลอมเหนือตัวเอง\n' ...
        '   bias ที่โตตามระยะจึงน่าจะมาจาก SNR ที่ลดลง (FISTA กด sidelobe ได้น้อยลง)\n' ...
        '   ไม่ใช่ phase error สะสม | ดูหัวไฟล์ W10_PhaseRefine.m\n']);
end

%% ===== (11) VALIDATION อิสระ: Camry showcase (grid/looks คนละชุด) =====
% เลือก estimator จากเก๋ง 6 คัน @ 72 looks / 0.125 m -> เสี่ยง overfit
% ชุด showcase คือ Camry เดียวกัน แต่ 144 looks / 0.10 m = ชุดตรวจอิสระ
scFile = fullfile(figDir, 'w9v3_cache_B.mat');
if exist(scFile, 'file')
    fprintf('\n===== ตรวจอิสระ: Camry showcase (%s) =====\n', 'w9v3_cache_B.mat');
    SC = load(scFile);                 % Vcs, axv, ctr3, hd, target, nLook
    iB2 = find(strcmp(veh(:,2), 'B'), 1);
    [pu2, pv2, pz2, pm2] = cloudFromVol(SC.Vcs, SC.axv, deg2rad(SC.hd), 40);
    inEx2 = abs(pu2) <= exR & abs(pv2) <= exR;
    m2a = inEx2 & pm2 >= -EST(iBL).dB;
    m2b = inEx2 & pm2 >= -EST(iBW).dB;
    L2 = (pctl(pv2(m2a),EST(iBL).phi) - pctl(pv2(m2a),EST(iBL).plo)) / EST(iBL).gain;
    W2 = (pctl(pu2(m2b),EST(iBW).phi) - pctl(pu2(m2b),EST(iBW).plo)) / EST(iBW).gain;
    inb2 = abs(pu2) <= GT_W(iB2)/2+0.15 & abs(pv2) <= GT_L(iB2)/2+0.15 & pz2>=zLo & pz2<=zHi;
    H2 = pctl(pz2(inb2 & pm2 >= -metDB), ZEST(iBz).p);
    fprintf('  grid %.3f m, %d looks (เทียบกับ %.3f m, %d looks ที่ใช้เลือก)\n', ...
        SC.axv(2)-SC.axv(1), SC.nLook, dVox, nLook);
    fprintf('  %-6s %8s %8s %8s %9s\n', 'มิติ', 'spec', '72look', '144look', 'ต่างกัน');
    fprintf('  %-6s %8.2f %8.2f %8.2f %9.2f\n', 'L', GT_L(iB2), Lm(iB2,iBL), L2, L2-Lm(iB2,iBL));
    fprintf('  %-6s %8.2f %8.2f %8.2f %9.2f\n', 'W', GT_W(iB2), Wm(iB2,iBW), W2, W2-Wm(iB2,iBW));
    fprintf('  %-6s %8.2f %8.2f %8.2f %9.2f\n', 'H', GT_H(iB2), Zm(iB2,iBz), H2, H2-Zm(iB2,iBz));
    fprintf(['  ถ้า "ต่างกัน" เล็ก (< ~1 voxel) = estimator ทนต่อการเปลี่ยน grid/looks\n' ...
        '  ถ้าใหญ่ = ค่าที่วัดได้ขึ้นกับพารามิเตอร์ recon ไม่ใช่ตัวรถ -> อย่าเพิ่งเชื่อ\n']);
else
    fprintf('\n(ข้าม validation: ไม่พบ %s)\n', scFile);
end

%% ===== (12) ตารางสุดท้าย พร้อม paste =====
fprintf('\n===== ตาราง 6 | ผลชุดใหม่ (L:%s, W:%s, H:%s) =====\n', ...
    EST(iBL).name, EST(iBW).name, ZEST(iBz).name);
fprintf('%-16s %-3s | %6s %6s %6s | %6s %6s %6s | %6s %6s %6s | %6s\n', ...
    'Vehicle','ID','specL','fitL','dL','specW','fitW','dW','specH','fitH','dH','in-box');
fprintf('%s\n', repmat('-', 1, 97));
for v = 1:Nveh
    fprintf('%-16s %-3s | %6.2f %6.2f %+6.2f | %6.2f %6.2f %+6.2f | %6.2f %6.2f %+6.2f | %5.0f%%\n', ...
        veh{v,1}, veh{v,2}, GT_L(v), Lm(v,iBL), dL(v,iBL), ...
        GT_W(v), Wm(v,iBW), dW(v,iBW), GT_H(v), Zm(v,iBz), dH(v,iBz), prec(v));
end
if any(~isVeh)
    fprintf('(G/H = ไม่ใช่รถ ไม่ถูกนับในค่าเฉลี่ยใด ๆ — ดู finding W9)\n');
end

%% ===== (13) FIGURES =====
f1 = figure(1); clf; set(f1,'Color','w','Position',[60 60 1080 400],'Name','W10 estimator');
tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
nexttile; bar([rL; rW].'); set(gca,'XTick',1:nE,'XTickLabel',{EST.name},'XTickLabelRotation',30);
grid on; ylabel('r กับ spec'); ylim([-1 1]); yline(0.8,'--'); legend({'L','W'},'Location','best');
title('ติดตามของจริงไหม (สูง = ดี)');
nexttile; bar([sdL; sdW].'); set(gca,'XTick',1:nE,'XTickLabel',{EST.name},'XTickLabelRotation',30);
grid on; ylabel('sd ของ Δ (m)'); yline(dVox,'--','1 voxel'); legend({'L','W'},'Location','best');
title('ความกระจาย (ต่ำ = ดี)');
nexttile; bar(rH,'FaceColor',[0.85 0.55 0.15]);
set(gca,'XTick',1:nZ,'XTickLabel',{ZEST.name},'XTickLabelRotation',30); grid on;
ylabel('r กับ roof spec'); ylim([-1 1]); title('ความสูง');
sgtitle(sprintf('W10 #1 | เลือก estimator ด้วย r ไม่ใช่ MAE — เก๋ง %d คัน, %d looks', ...
    nnz(isCar), nLook));
exportgraphics(f1, fullfile(figDir,'fig_w10_est_quality.png'), 'Resolution', 150);

f2 = figure(2); clf; set(f2,'Color','w','Position',[60 60 1120 420],'Name','W10 fit vs spec');
tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
scatterFit(GT_L(isCar), [Lm(isCar,iBase) Lm(isCar,iBL)], 'L (m)', ...
    {EST(iBase).name, EST(iBL).name});
scatterFit(GT_W(isCar), [Wm(isCar,iBase) Wm(isCar,iBW)], 'W (m)', ...
    {EST(iBase).name, EST(iBW).name});
scatterFit(GT_H(isCar), [Zm(isCar,1) Zm(isCar,iBz)], 'ความสูง (m)', ...
    {ZEST(1).name, ZEST(iBz).name});
sgtitle('W10 #1 | วัดได้ vs spec — จุดควรเกาะเส้น 1:1 ถ้าวัดของจริง');
exportgraphics(f2, fullfile(figDir,'fig_w10_fit_vs_spec.png'), 'Resolution', 150);

f3 = figure(3); clf; set(f3,'Color','w','Position',[60 60 1150 460],'Name','W10 roofline');
cols = lines(Nveh); hold on;
for v = 1:Nveh
    if ~isVeh(v), continue; end
    stairs(vbc, roofM(v,:), ':', 'Color', min(cols(v,:)+0.35,1), 'LineWidth', 1, ...
        'HandleVisibility','off');
    stairs(vbc, roofP(v,:), 'Color', cols(v,:), 'LineWidth', 1.5 + (GT_H(v)>=2), ...
        'DisplayName', sprintf('%s (%s)', veh{v,1}, veh{v,2}));
end
hold off; grid on; xlim([-3 3]); ylim([-0.2 3.0]);
xlabel('v — ตามยาวรถ (m)'); ylabel('ความสูง (m)'); legend('Location','eastoutside');
title(sprintf('W10 #1 | roofline p90 (ทึบ) vs max เดิม (จุดจาง) — %d looks', nLook));
exportgraphics(f3, fullfile(figDir,'fig_w10_roofline_p90.png'), 'Resolution', 150);

%% ===== (14) EXPORT =====
mdFile = fullfile(figDir, 'w10_robust_table.md');
fid = fopen(mdFile, 'w', 'n', 'UTF-8');
fprintf(fid, '# W10 #1 — Robust metrics (%d looks, voxel %.3f m)\n\n', nLook, dVox);
fprintf(fid, 'estimator: L = **%s**, W = **%s**, roof = **%s**\n\n', ...
    EST(iBL).name, EST(iBW).name, ZEST(iBz).name);
fprintf(fid, '| ID | รถ | spec L | วัด L (Δ) | spec W | วัด W (Δ) | spec H | วัด H (Δ) | in-box |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|---|\n');
for v = 1:Nveh
    if ~isVeh(v), continue; end
    fprintf(fid, '| %s | %s | %.2f | %.2f (%+.2f) | %.2f | %.2f (%+.2f) | %.2f | %.2f (%+.2f) | %.0f%% |\n', ...
        veh{v,2}, veh{v,1}, GT_L(v), Lm(v,iBL), dL(v,iBL), GT_W(v), Wm(v,iBW), dW(v,iBW), ...
        GT_H(v), Zm(v,iBz), dH(v,iBz), prec(v));
end
fprintf(fid, '\n## คุณภาพ estimator (เก๋ง %d คัน) — r = ติดตาม spec, sd = ความกระจาย\n\n', nnz(isCar));
fprintf(fid, '### ความยาว L\n\n| estimator | r | slope | sd | MAE |\n|---|---|---|---|---|\n');
for e = 1:nE
    fprintf(fid, '| %s%s | %+.3f | %+.2f | %.3f | %.2f |\n', EST(e).name, ...
        ternary(e==iBase,' (W9 เดิม)',''), rL(e), slL(e), sdL(e), mean(abs(dL(isCar,e)),'omitnan'));
end
fprintf(fid, '\n### ความกว้าง W\n\n| estimator | r | slope | sd | MAE |\n|---|---|---|---|---|\n');
for e = 1:nE
    fprintf(fid, '| %s%s | %+.3f | %+.2f | %.3f | %.2f |\n', EST(e).name, ...
        ternary(e==iBase,' (W9 เดิม)',''), rW(e), slW(e), sdW(e), mean(abs(dW(isCar,e)),'omitnan'));
end
fprintf(fid, '\n### ความสูง\n\n| estimator | r | slope | sd | MAE |\n|---|---|---|---|---|\n');
for k = 1:nZ
    fprintf(fid, '| %s%s | %+.3f | %+.2f | %.3f | %.2f |\n', ZEST(k).name, ...
        ternary(k==1,' (W9 เดิม)',''), rH(k), slH(k), sdH(k), mean(abs(dH(isCar,k)),'omitnan'));
end
fprintf(fid, '\n> เพดานความละเอียด: voxel %.3f m | ช่วง spec ของเก๋ง L %.2f m, W %.2f m, H %.2f m\n', ...
    dVox, max(GT_L(isCar))-min(GT_L(isCar)), max(GT_W(isCar))-min(GT_W(isCar)), ...
    max(GT_H(isCar))-min(GT_H(isCar)));
fclose(fid);

save(fullfile(figDir, 'w10_robust_metrics.mat'), 'Lm','Wm','Zm','dL','dW','dH', ...
    'EST','ZEST','iBL','iBW','iBz','iBase','veh','GT_L','GT_W','GT_H','isCar','isVeh', ...
    'rL','slL','sdL','rW','slW','sdW','rH','slH','sdH','nUse','nZuse', ...
    'roofP','roofM','roofBad','roofCab','vbc','outFracU','tailU','shapeRu','shapeRv', ...
    'prec','nLook','dVox');

fprintf('\nเซฟแล้ว: %s\n  + w10_robust_metrics.mat + fig_w10_*.png (3 รูป)\n', mdFile);

%% ===== LOCAL FUNCTIONS =====
function y = pctl(x, p)
% percentile แบบเดียวกับ prctile ของ MATLAB | p=0 -> min, p=100 -> max เป๊ะ
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

function r = spanRatio(x)
% span(p5-95)/span(p10-90) | จุดสม่ำเสมอ -> 1.125 | สูงกว่า = หางหนา (ghost)
s = pctl(x, 90) - pctl(x, 10);
if s > 0, r = (pctl(x, 95) - pctl(x, 5)) / s; else, r = NaN; end
end

function [r, slope] = corrSlope(gt, fit)
% correlation + ความชัน ระหว่างค่าที่วัดกับ spec (ตัด NaN)
ok = isfinite(gt) & isfinite(fit);
if nnz(ok) < 3 || std(gt(ok)) == 0, r = NaN; slope = NaN; return; end
c = corrcoef(gt(ok), fit(ok));  r = c(1,2);
p = polyfit(gt(ok), fit(ok), 1); slope = p(1);
end

function i = pickBest(r, sd, mae)
% เลือกตัวที่ "ติดตามของจริง" ก่อน (r >= 0.8) แล้วเอาตัว sd ต่ำสุดในกลุ่มนั้น
% ถ้าไม่มีตัวไหนถึง 0.8 = มิตินั้นแยกคันไม่ได้อยู่แล้ว -> เปลี่ยนไปเอา "ค่าสัมบูรณ์แม่นสุด"
% (MAE ต่ำสุด) เพราะยังใช้บอกขนาดคร่าว ๆ / แยกข้ามชั้น (เก๋ง vs เครื่องจักร) ได้
good = find(r >= 0.8);
if isempty(good)
    [~, i] = min(mae);
else
    [~, j] = min(sd(good));  i = good(j);
end
end

function s = markOf(e, iBest, iBase)
% ใส่คอมมาให้ชัด: [s(1) '*'] กำกวม (MATLAB อาจอ่าน ' เป็น transpose)
s = '  ';
if e == iBest, s = '->'; end
if e == iBase, s = [s(1), '*']; end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function [pu, pv, pz, pm] = cloudFromVol(Vcs, axv, hdRad, topDB)
% สร้าง point cloud จาก voxel volume ด้วยกฎเดียวกับ W9 (top-N dB + ตัด voxel โดดเดี่ยว)
Mn     = Vcs / max(Vcs(:));
selTop = Mn >= 10^(-topDB/20);
nb     = convn(double(selTop), ones(3,3,3), 'same') - selTop;
sel    = selTop & nb >= 1;
[XX, YY, ZZ] = ndgrid(axv, axv, axv);
Rr = [cos(hdRad) -sin(hdRad); sin(hdRad) cos(hdRad)];
uv = Rr.' * [XX(sel).'; YY(sel).'];
pu = uv(1,:).';  pv = uv(2,:).';  pz = ZZ(sel);  pm = 20*log10(Mn(sel));
end

function scatterFit(gt, fits, ylab, names)
nexttile; hold on;
mk = {'o','s'};  cl = [0.75 0.35 0.15; 0.10 0.45 0.70];
lo = min([gt(:); fits(:)]) - 0.1;  hi = max([gt(:); fits(:)]) + 0.1;
plot([lo hi], [lo hi], 'k--', 'DisplayName', '1:1');
for j = 1:size(fits, 2)
    plot(gt, fits(:,j), mk{j}, 'MarkerFaceColor', cl(j,:), 'MarkerEdgeColor','none', ...
        'MarkerSize', 8, 'DisplayName', names{j});
end
hold off; grid on; axis equal; xlim([lo hi]); ylim([lo hi]);
xlabel('spec (m)'); ylabel(ylab); legend('Location','northwest'); title(ylab);
end

function printTable(ttl, fit, gt, ids, names, isCar, isVeh, nLook, iBase)
% แถว = estimator, คอลัมน์ = รถ | แสดง Δ + สรุปท้ายแถว
d = fit - gt;
nV = numel(ids);  nR = numel(names);
fprintf('\n===== %s | Δ = วัดได้ − spec (m) | %d looks =====\n', ttl, nLook);
fprintf('%-14s', 'estimator'); fprintf('%7s', ids{:});
fprintf(' |%7s%7s%7s%7s\n', 'MAEcar', 'bias', 'sd', 'r');
fprintf('%-14s', 'spec (m)');  fprintf('%7.2f', gt); fprintf(' |\n');
fprintf('%s\n', repmat('-', 1, 14 + 7*nV + 30));
for r = 1:nR
    rr = corrSlope(gt(isCar), fit(isCar, r));
    fprintf('%s%-12s', markOf(r, 0, iBase), names{r});
    fprintf('%+7.2f', d(:, r));
    fprintf(' |%7.2f%+7.2f%7.3f%+7.3f\n', mean(abs(d(isCar,r)),'omitnan'), ...
        mean(d(isCar,r),'omitnan'), std(d(isCar,r),'omitnan'), rr);
end
fprintf(['* = ที่ W9 ใช้ | MAEcar/bias/sd/r คิดจากเก๋งเท่านั้น (%d คัน)\n' ...
    'r = correlation กับ spec — ตัวชี้ว่า "วัดของจริง" หรือ "คืนค่ากลาง" (MAE ต่ำอย่างเดียวหลอกได้)\n'], ...
    nnz(isCar));
if any(~isVeh)
    fprintf('G/H ไม่ถูกนับในสถิติใด ๆ (ไม่ใช่รถ — finding W9)\n');
end
end
