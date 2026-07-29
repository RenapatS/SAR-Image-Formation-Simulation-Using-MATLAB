%% W9 v2 — SPARSE 3D บน GOTCHA ตามวิธีจริงของเปเปอร์ (k-space regularized L1-LS)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%
%  ทำไมต้อง v2: ผลของ v1 (สำรองไว้ที่ W9_GotchaSparse3D_v1backup.m) ไม่เป็นรูปรถ
%  เพราะ v1 ใช้วิธี "per-pixel elevation inversion" (ทำ 2D image ราย pass แล้วแก้
%  L1 ทีละพิกเซลในแนว elevation) — นั่นคือวิธีที่ 2 ของเปเปอร์ (multipass IFSAR,
%  SPIE2009 Sec.5) ซึ่งเปเปอร์เองสรุปว่าให้ "higher variability in 3D position"
%  ส่วนรูป Camry สวย ๆ (Fig 7 ที่เราอยากได้) มาจากวิธีที่ 1:
%
%    REGULARIZED Lp LS บน K-SPACE (SPIE2009 Sec.4 / JSTSP2011 Sec.III):
%    ต่อ 1 subaperture 5 องศา รวมข้อมูลทั้ง 8 pass เป็นก้อนเดียว แล้วแก้
%        min || y - A x ||^2 + lambda*||x||_1
%    ตรง ๆ บน voxel grid 3 มิติ (A = 3D Fourier จาก voxel -> k-space, ทำเร็ว
%    ด้วย nearest-neighbor gridding + FFT) แล้วรวม 72 subaperture x pol แบบ
%    noncoherent ด้วย MAX (eq.7) — แสดงผล "top 40 dB ของ voxel" แบบเปเปอร์
%
%  จุดที่ v1 พังเพิ่มเติม (เจอจากการวิเคราะห์ + จำลอง geometry จริง):
%   (1) elevation จริง 8 pass ไม่ uniform: 44.27,44.18,44.10,44.01,43.92,
%       43.53,43.01,43.06 องศา (เปเปอร์รายงานเอง) — ช่องว่าง 0.39/0.52 องศา
%       สร้าง near-ambiguity ที่ s~1.9-2.1 m -> "แฝดปลอม" z สูง 1.4 m
%       (คือหมอกที่ลอยเหนือหลังคาในรูป v1); การแก้รายพิกเซลด้วย 8 sample
%       เดี่ยว ๆ ต้านสิ่งนี้ไม่ไหว แต่การแก้ร่วมทั้ง volume บน k-space ไหว
%   (2) debias step (LS refit บน support, cond<50) ขยาย noise ได้ถึง 50 เท่า
%       ราย look -> โดน max-combining เก็บ -> จุดหลอนสว่าง
%   (3) lambda ตั้งจาก max ทั้ง patch -> พิกเซลอ่อนโดน threshold ตาย ->
%       point cloud โหรงเหรง
%   (4) steering ใช้ kc = 4*pi*freq(1)/c แต่เฟสข้าม pass ของภาพ BP จริง ๆ
%       ควบคุมโดยความถี่กลาง (9.6 GHz ไม่ใช่ 9.28) -> สเกลแกนสูงเพี้ยน ~3.3%
%   (5) แสดงผลที่ -10/-13 dB แต่เปเปอร์โชว์ top 40 dB
%
%  วิธีใหม่ตรวจแล้วบน synthetic replica ของ GOTCHA (elevation จริงทั้ง 8,
%  9.28-9.92 GHz, R~10.2 km, รถ+พื้น+รถข้าง ๆ นอกกรอบ):
%  median localization error = 0.09-0.10 m (~= voxel), ghost = 0
%
%  อ้างอิง:
%   [1] Austin, Ertin, Moses, "Sparse multipass 3D SAR imaging: applications
%       to the GOTCHA data set," Proc. SPIE 7337 (2009) — Sec.4, Fig.6-7
%   [2] Austin, Ertin, Moses, "Sparse Signal Methods for 3-D Radar Imaging,"
%       IEEE JSTSP 5(3) 2011 — Sec.V-B, Fig.12-13
%   พารามิเตอร์ตามเปเปอร์: subaperture 5 deg x 72 อัน (ไม่ overlap),
%   voxel 0.1 m, กล่อง [-5,5)^3 m, k-grid 100^3 (bandwidth 62.8318 rad/m),
%   p=1, รวมภาพด้วย max ข้าม subaperture+pol, แสดง top 40 dB
%
%  SPOTLIGHT: เปเปอร์ทำ k-space spotlight filtering กันสิ่งอื่นในฉาก alias
%  เข้ากรอบ 10 m (เพราะ NN downsample) — เราใช้ range gate + cross-range
%  (Doppler) gate รอบรถแทน (ผลเท่ากัน, ง่ายกว่า) | หมายเหตุ: ground offset
%  ถูกบีบเป็น slant ด้วย cos(elev~44) ~ 0.72 -> gate slant 4.0 m ~ ground 5.6 m
%
%  knobs หลัก: target ('B'=Camry ตรงกับเปเปอร์ | 'A'=Malibu | 'TOPHAT' ฯลฯ)
%              polList {'HH','VV'} (เปเปอร์บอกภาพรวม 2 pol สวยสุด)
%              subApStep (5 = ครบ 72 อันแบบเปเปอร์ | 10 = เร็วขึ้น 2 เท่า)
%              NsigCFAR (20 default; หมอกเยอะ -> 30, จุดน้อยไป -> 12)
%  RUNTIME: ~30-60 นาที (2 pol, step 5) | โหลดไฟล์คือคอขวดรอบแรก (cold cache)
%  ** ตัวเลขทุกตัวใน EVALUATION มาจากการรันจริง **

clear; clc; close all;
c = physconst('LightSpeed');

%% ===== (0) CONFIG =====
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % root ของ repo
discRoot = fullfile(repoRoot, 'data');  % ที่วางข้อมูล GOTCHA
assert(exist(fullfile(discRoot,'GOTCHA-CP_Disc1','DATA'),'dir') > 0, ...
    'ไม่พบโฟลเดอร์ GOTCHA-CP_Disc1 ใน %s', discRoot);

% เป้า: {id, name, cx, cy, heading(deg), L, W, roofSpec}
vehTab = {'A','ChevyMalibu',    9.9696,  -5.2239,   3.40747, 4.77, 1.74, 1.45;
          'B','ToyotaCamry',   20.6630, -18.7070, 182.80,    4.75, 1.74, 1.43;
          'C','FordTaurusWag', 12.4260, -18.2140, 185.04,    4.98, 1.86, 1.47;
          'F','HyundaiSantaFe',29.2360, -19.1770, 184.03,    4.45, 1.77, 1.68};
target = 'B';                    % 'B' = Toyota Camry (เทียบตรงกับเปเปอร์ Fig 7)
                                 % หรือ 'TOPHAT' = ตรวจ pipeline กับ calibration
                                 % target (ต้องได้วงแหวนรัศมี 1 m แบบ Fig 6)
polList    = {'HH','VV'};        % รวมแบบ max ข้าม pol (เปเปอร์ eq.7)
subApWidth = 5;                  % องศา/subaperture (เปเปอร์ = 5)
subApStep  = 5;                  % ระยะเลื่อน (5 = ครบ 72 อัน, 10 = เร็ว 2 เท่า)
% threshold ของ L1: อิง NOISE FLOOR ราย look (CFAR-style) ไม่อิง peak!
% บทเรียนจากรัน TOPHAT(ผ่าน) vs CAMRY(หมอก): tophat สว่างเท่ากันทุก aspect
% แต่รถมี broadside glint แรงมากเป็นบาง look -> ถ้า thr = %ของ peak:
% look ที่มี glint -> thr สูงเว่อร์ (เหลือแต่ glint), look มืด -> thr ต่ำ
% (เก็บ clutter เพียบ) แล้ว max-combine 144 look รวม junk เป็นหมอกเต็มกล่อง
% เปเปอร์ใช้ lambda=10 คงที่ทุก subaperture (absolute) = เจตนาเดียวกัน
NsigCFAR   = 30;                 % thr = NsigCFAR * sigma_noise(ราย look)
                                 % (เพิ่ม -> สะอาด/บางลง, ลด -> จุดแน่น/หมอกมากขึ้น)
relFloor   = 0.02;               % กันเศษ sidelobe ของ glint: thr >= 2% ของ peak
nIterFISTA = 200;
gateRng    = 3.5;                % slant range gate (m) ~ ground 4.9 m
gateCrs    = 4.5;                % cross-range gate (m)
% (แคบพอฆ่าเพื่อนบ้านห่าง ~5-6 m ก่อนมัน wrap เข้ากล่อง 10 m — บทเรียนจาก
%  รัน 11 คัน: Jetta โดน Malibu ที่ห่าง 5.5 m wrap มาทับจน in-box = 0%)
topDB      = 40;                 % แสดง top-40dB voxels แบบเปเปอร์

% voxel/k grid ตามเปเปอร์เป๊ะ: 0.1 m, กล่อง 10 m, 100^3
N  = 100;  dr = 0.10;
dk = 2*pi/(N*dr);                % = 0.62832 (N*dr*dk = 2*pi พอดี)
axv = ((0:N-1) - N/2) * dr;      % แกน voxel (ศูนย์กลางที่เป้า)

if strcmpi(target,'TOPHAT')
    ctr = [-15.91, 20.93];  hd = 0;  gtL = 2.0; gtW = 2.0; roofSpec = 0.75;
    tName = 'Tophat (calibration)';
else
    iv = find(strcmp(vehTab(:,1), target), 1);
    assert(~isempty(iv), 'target ต้องเป็น A/B/C/F หรือ TOPHAT');
    ctr = [vehTab{iv,3}, vehTab{iv,4}];  hd = vehTab{iv,5};
    gtL = vehTab{iv,6};  gtW = vehTab{iv,7};  roofSpec = vehTab{iv,8};
    tName = sprintf('%s (%s)', vehTab{iv,2}, target);
end
ctr3 = [ctr, 0];

%% ===== (1) probe passes + freq =====
passList = [];
for pn = 1:8
    dsc = 1 + (pn == 8);
    if exist(fullfile(discRoot, sprintf('GOTCHA-CP_Disc%d',dsc), 'DATA', ...
             sprintf('pass%d',pn), 'HH'), 'dir')
        passList(end+1) = pn; %#ok<SAGROW>
    end
end
assert(numel(passList) >= 4, 'พบแค่ %d pass (ต้อง >= 4)', numel(passList));
S0 = load(fullfile(discRoot,'GOTCHA-CP_Disc1','DATA','pass1','HH', ...
                   'data_3dsar_pass1_az001_HH.mat'));
freq = double(S0.data.freq(:));  Nf = numel(freq);
fprintf('W9v2 (paper method) | %s | passes %s | %d freq %.2f-%.2f GHz\n', ...
    tName, mat2str(passList), Nf, freq(1)/1e9, freq(end)/1e9);
nSub = numel(0:subApStep:359);
fprintf('subap %d deg x %d อัน x %d pol | grid %d^3 @ %.2f m | thr %.0f*sigma\n', ...
    subApWidth, nSub, numel(polList), N, dr, NsigCFAR);

%% ===== (2) MAIN LOOP: subaperture -> k-space grid -> FISTA -> max-combine =====
Vcs = zeros(N,N,N);              % sparse (L1)  — max ข้าม look (เปเปอร์ eq.7)
Vft = zeros(N,N,N);              % Fourier baseline (adjoint) — เทียบแบบ Fig 11
nLook = 0;  tAll = tic;  kSamp = [];
for aDeg = 0 : subApStep : 359
  for ip = 1:numel(polList)
    pol = polList{ip};
    passes = loadSubap(discRoot, passList, pol, aDeg, subApWidth);
    if numel(passes) < 4, continue; end
    [yg, mk, kS] = gridSubap(passes, ctr3, N, dk, gateRng, gateCrs, c);
    [xS, xF] = fistaKspace(yg, mk, nIterFISTA, NsigCFAR, relFloor);
    Vcs = max(Vcs, abs(xS));
    Vft = max(Vft, abs(xF));
    nLook = nLook + 1;
    if isempty(kSamp), kSamp = kS; end               % เก็บตัวอย่างไว้ทำ Fig 1
    fprintf('  az %3d %s | %d passes | nnz %.2f%% | %.0f s (ETA %.0f min)\n', ...
        aDeg, pol, numel(passes), 100*nnz(xS)/numel(xS), toc(tAll), ...
        toc(tAll)/nLook*(nSub*numel(polList)-nLook)/60);
  end
end
fprintf('รวม %d looks | total %.1f min\n', nLook, toc(tAll)/60);

%% ===== (3) METRICS (คิดใน car frame: u=กว้าง, v=ยาว) =====
Mn  = Vcs / max(Vcs(:));
Mf  = Vft / max(Vft(:));
thH = deg2rad(hd);  Rr = [cos(thH) -sin(thH); sin(thH) cos(thH)];
[XX, YY, ZZ] = ndgrid(axv, axv, axv);                % [ix,iy,iz] ตาม gridSubap

selTop = Mn >= 10^(-topDB/20);
% display denoise: ตัด voxel โดดเดี่ยว (ไม่มีเพื่อนบ้านใน 26-neighborhood)
nb  = convn(double(selTop), ones(3,3,3), 'same') - selTop;
selShow = selTop & nb >= 1;
uvw = Rr.' * [XX(selShow).'; YY(selShow).'];          % car frame
pu = uvw(1,:).';  pv = uvw(2,:).';  pz = ZZ(selShow);
pm = 20*log10(Mn(selShow));

inBox = abs(pu) <= gtW/2+0.15 & abs(pv) <= gtL/2+0.15 & pz >= -0.5 & pz <= 2.6;
sel20 = pm >= -20;                                    % ชุดเข้มไว้วัดขนาด
Lm = 0; Wm = 0;
if any(sel20)
    Lm = max(pv(sel20)) - min(pv(sel20));  Wm = max(pu(sel20)) - min(pu(sel20));
end
zin = pz(inBox & sel20);

fprintf('\n========== EVALUATION — W9v2 (%s) ==========\n', tName);
fprintf('[top %d dB] %d voxels (แสดง %d หลังตัดโดดเดี่ยว) | ในกรอบรถ %.0f%%\n', ...
    topDB, nnz(selTop), nnz(selShow), 100*mean(inBox));
fprintf('[extent @-20dB] L x W = %.2f x %.2f m (GT %.2f x %.2f)\n', Lm, Wm, gtL, gtW);
if ~isempty(zin)
    fprintf('[height ในกรอบ @-20dB] z median %.2f / max %.2f m (roof spec ~%.2f)\n', ...
        median(zin), max(zin), roofSpec);
end
fprintf('[หมายเหตุ] Fourier baseline เทียบแบบเปเปอร์ Fig 11 อยู่ใน Fig 2 ล่าง\n');

%% ===== (4) FIGURES =====
figDir = fullfile(repoRoot, 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end

% --- Fig 1: k-space sampling ของ 1 subaperture (แบบเปเปอร์ Fig 2) ---
figure(1); set(gcf,'Name','Fig1: k-space','Color','w','Position',[60 60 900 400]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
ds = 1:37:size(kSamp,1);                              % thin ให้พอ plot ไหว
plot(sqrt(kSamp(ds,1).^2+kSamp(ds,2).^2), kSamp(ds,3), '.', 'MarkerSize', 2);
grid on; xlabel('k_{xy} (rad/m)'); ylabel('k_z (rad/m)'); axis tight;
title('k-space: 8 pass = 8 เส้น elevation (ไม่ uniform!)');
nexttile;
plot(kSamp(ds,1), kSamp(ds,2), '.', 'MarkerSize', 2);
grid on; xlabel('k_x'); ylabel('k_y'); axis equal tight;
title(sprintf('มุมมองบน: arc %d deg', subApWidth));
sgtitle('Fig 1 | W9v2 — k-space ที่เก็บจริงต่อ 1 subaperture');

% --- Fig 2: SPARSE vs Fourier (side + top MIP) ---
figure(2); set(gcf,'Name','Fig2: sparse vs Fourier','Color','w','Position',[80 80 1100 700]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ttl = {sprintf('SPARSE L1 (วิธีเปเปอร์) — side MIP'), 'SPARSE — top MIP'; ...
       'Fourier baseline — side MIP (ฟุ้งแบบ Fig 11)', 'Fourier — top MIP'};
for r = 1:2
    if r==1, M3 = Mn; else, M3 = Mf; end
    sideM = squeeze(max(M3, [], 1)).';                % (z ต่อ y): max ข้าม x
    topM  = squeeze(max(M3, [], 3)).';                % (y ต่อ x)
    ax1 = nexttile; imagesc(axv, axv, 20*log10(sideM+1e-9));
    set(ax1,'YDir','normal'); clim([-topDB 0]); colormap(ax1,'jet');
    hold on; yline(0,'w:'); yline(roofSpec,'w--'); hold off;
    xlabel('y (m)'); ylabel('z (m)'); ylim([-2 3]); title(ttl{r,1});
    ax2 = nexttile; imagesc(axv, axv, 20*log10(topM+1e-9));
    set(ax2,'YDir','normal'); axis image; clim([-topDB 0]); colormap(ax2,'gray');
    xlabel('x (m)'); ylabel('y (m)'); title(ttl{r,2});
end
cb = colorbar; ylabel(cb,'dB'); cb.Layout.Tile = 'east';
sgtitle(sprintf('Fig 2 | W9v2 — %s: sparse vs Fourier (%d looks)', tName, nLook));

% --- Fig 3: paper-style point cloud 3 มุมมอง (top 40 dB) ---
figure(3); set(gcf,'Name','Fig3: paper-style views','Color','w','Position',[100 100 1400 430]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
sz = 4 + 46*max(0, (pm + topDB)/topDB).^1.5;
if strcmpi(target,'TOPHAT')
    gtx = cos(linspace(0,2*pi,60)); gty = sin(linspace(0,2*pi,60));
else
    gtx = [-gtW/2 gtW/2 gtW/2 -gtW/2 -gtW/2];
    gty = [-gtL/2 -gtL/2 gtL/2 gtL/2 -gtL/2];
end
for vw = 1:3
    ax3 = nexttile;
    scatter3(pu, pv, pz, sz, pz, 'filled', 'MarkerFaceAlpha', 0.75); hold on;
    plot3(gtx, gty, zeros(size(gtx)), '-', 'Color', [.55 .55 .55], 'LineWidth', 1.2);
    hold off; colormap(ax3, turbo); clim([-0.3 2.2]);
    xlabel('u กว้าง (m)'); ylabel('v ยาว (m)'); zlabel('z (m)');
    daspect([1 1 1]); grid on; box on;
    xlim([-3.5 3.5]); ylim([-3.5 3.5]); zlim([-1 2.8]);
    switch vw
        case 1, view(-40, 22); title('(a) 3D view');
        case 2, view(90, 0);   title('(b) Side view');
        case 3, view(0, 90);   title('(c) Top view');
    end
end
cb = colorbar; ylabel(cb,'z (m)'); cb.Layout.Tile = 'east';
sgtitle(sprintf('Fig 3 | W9v2 — %s: regularized L1-LS, top %d dB (แบบเปเปอร์ Fig 7)', ...
    tName, topDB));

exportgraphics(figure(1), fullfile(figDir,'fig_w9v3_kspace.png'),    'Resolution',150);
exportgraphics(figure(2), fullfile(figDir,'fig_w9v3_cs_vs_fourier.png'),'Resolution',150);
exportgraphics(figure(3), fullfile(figDir,sprintf('fig_w9v3_paperstyle_%s.png',target)),'Resolution',150);
save(fullfile(figDir, sprintf('w9v3_cache_%s.mat',target)), ...
     'Vcs','Vft','axv','ctr3','hd','target','nLook','-v7.3');
fprintf('\nFigures + cache saved to %s\n', figDir);


%%%% ===== LOCAL FUNCTIONS =====

function passes = loadSubap(discRoot, passList, pol, aDeg, nFiles)
% โหลด nFiles ไฟล์ (1 deg/ไฟล์) ต่อ pass + ใส่ autofocus (ph_correct, r_correct)
% หมายเหตุ: จำนวน freq sample (Nf) ไม่เท่ากันทุกไฟล์/ทุก pass -> เก็บ freq
% ราย pass แล้วให้ gridSubap ใช้ของใครของมัน (ไฟล์ที่ Nf ต่างจากไฟล์แรก
% ของ pass นั้นจะถูกข้าม — เกิดน้อยมาก)
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
                continue;                            % Nf ไม่ตรงกับไฟล์แรกของ pass
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

function [yg, mk, kSamp] = gridSubap(passes, ctr3, N, dk, gateRng, gateCrs, c)
% recenter เฟสไปที่เป้า (exact ราย pulse) -> range/cross-range gate (spotlight)
% -> NN-bin ลง k-grid กลางศูนย์ (เฉลี่ยต่อ cell) | คืน yg (N^3), mask, ตัวอย่าง k
    KX=[]; KY=[]; KZ=[]; VL=[];
    for q = 1:numel(passes)
        fp = passes(q).fp;  ant = passes(q).ant;  r0u = passes(q).r0used;
        freq = passes(q).freq;                              % ราย pass!
        Nf = numel(freq);  df = freq(2) - freq(1);
        P = size(fp, 2);
        Rc = sqrt(sum((ant - ctr3).^2, 2)).';               % 1xP exact
        fpc = fp .* exp(1j*(4*pi/c)*freq*(Rc - r0u));       % NfxP recenter
        % --- range gate (slant) รอบเป้า ---
        Nz = Nf*4;
        rc = fftshift(ifft(fpc, Nz, 1), 1);
        drax = ((0:Nz-1).' - Nz/2) * (c/(2*df)) / Nz;
        rc(abs(drax) > gateRng, :) = 0;
        sp = fft(ifftshift(rc, 1), [], 1);
        fpc = sp(1:Nf, :);
        % --- cross-range (Doppler) gate ---
        if P >= 64
            Sp = fft(fpc, [], 2);
            fd = (0:P-1)/P;  fd(fd >= 0.5) = fd(fd >= 0.5) - 1;
            azp = unwrap(atan2(ant(:,2), ant(:,1)));
            dphi = abs(azp(end) - azp(1)) / max(P-1, 1);
            uax = fd * 2*pi / ((4*pi*mean(freq)/c) * max(dphi, 1e-9));
            Sp(:, abs(uax) > gateCrs) = 0;
            fpc = ifft(Sp, [], 2);
        end
        % --- ตำแหน่ง k-space ราย sample: k = (4*pi*f/c)*u_hat(pulse) ---
        U = ant - ctr3;  U = U ./ sqrt(sum(U.^2, 2));
        kf = (4*pi/c) * freq;
        KX = [KX; reshape(kf*U(:,1).', [], 1)]; %#ok<AGROW>
        KY = [KY; reshape(kf*U(:,2).', [], 1)]; %#ok<AGROW>
        KZ = [KZ; reshape(kf*U(:,3).', [], 1)]; %#ok<AGROW>
        VL = [VL; fpc(:)];                      %#ok<AGROW>
    end
    kSamp = [KX, KY, KZ];
    K0 = mean(kSamp, 1);
    ix = round((KX - K0(1))/dk) + N/2 + 1;
    iy = round((KY - K0(2))/dk) + N/2 + 1;
    iz = round((KZ - K0(3))/dk) + N/2 + 1;
    ok = ix>=1 & ix<=N & iy>=1 & iy<=N & iz>=1 & iz<=N;
    lin = sub2ind([N N N], ix(ok), iy(ok), iz(ok));
    cnt = accumarray(lin, 1,            [N^3 1]);
    yre = accumarray(lin, real(VL(ok)), [N^3 1]);
    yim = accumarray(lin, imag(VL(ok)), [N^3 1]);
    mk = cnt > 0;
    yg = zeros(N^3, 1);
    yg(mk) = (yre(mk) + 1j*yim(mk)) ./ cnt(mk);
    yg = reshape(yg, [N N N]);  mk = reshape(mk, [N N N]);
end

function [x, xF] = fistaKspace(yg, mk, nIter, Nsig, relFloor)
% แก้ min ||mask.*(F x) - y||^2 + lam||x||_1 ด้วย FISTA
% F = centered 3D DFT (e^{+j k.r}) ผ่าน FFT (ortho) -> sigma_max(A) <= 1
% threshold อิง noise floor ราย look: sigma จาก median|A^H y| (Rayleigh:
% median = sigma*sqrt(ln 4) = 1.1774*sigma) -> thr = Nsig*sigma
% เสถียรข้าม look ไม่ว่าจะมี glint หรือไม่ (เทียบเท่า lambda คงที่ของเปเปอร์)
    d  = yg .* mk;  Nc = sqrt(numel(d));
    Fw = @(x) fftshift(ifftn(ifftshift(x))) * Nc;      % forward  e^{+jkr}
    Aj = @(y) fftshift( fftn(ifftshift(y))) / Nc;      % adjoint  e^{-jkr}
    xF = Aj(d);                                        % Fourier baseline
    sig = median(abs(xF(:))) / 1.1774;
    thr = max(Nsig*sig, relFloor*max(abs(xF(:))));     % = step*lam, step=0.5
    x = zeros(size(d));  z = x;  t = 1;
    for it = 1:nIter
        w  = z - Aj(Fw(z).*mk - d);                    % step*2*A^H(...) = 1*
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
