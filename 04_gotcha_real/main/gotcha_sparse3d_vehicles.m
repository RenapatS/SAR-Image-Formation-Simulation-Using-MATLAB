%% W9b v2 — VEHICLE-TYPE DISCRIMINATION ด้วยวิธีจริงของเปเปอร์ (k-space L1-LS)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%
%  v1 (สำรองที่ W9_CompareVehicles_v1backup.m) ใช้ per-pixel elevation L1 ->
%  ภาพไม่เป็นรูปรถ (เหตุผลเต็ม ๆ อยู่หัวไฟล์ W9_GotchaSparse3D.m)
%  v2 เปลี่ยน engine เป็น regularized L1-LS บน k-space ต่อ subaperture 5 องศา
%  แล้ว max-combine ข้าม subaperture + polarization (Austin/Ertin/Moses:
%  SPIE 7337 (2009) Sec.4 + IEEE JSTSP 5(3) 2011 Sec.V-B)
%
%  โครงรัน: ต่อ 1 subaperture โหลดไฟล์ครั้งเดียว (8 pass x 5 ไฟล์ x pol)
%  แล้ววน 11 คัน: recenter->gate->NN k-grid->FISTA->max-combine ราย voxel
%  metric วัดใน CAR FRAME (u=กว้าง, v=ยาว ตาม heading จริงจาก xls)
%
%  knobs: subApStep 10 = ใช้ 36 จาก 72 subaperture (default, ~ครึ่งเวลา)
%                    5 = ครบแบบเปเปอร์ (ภาพแน่นขึ้น, ช้า 2 เท่า)
%         polList {'HH','VV'} แบบเปเปอร์ | {'HH'} = เร็ว 2 เท่า
%  RUNTIME (cache อุ่น): step10+2pol ~ 1-1.5 ชม. | step10+HH ~ 30-45 นาที
%  ** ตัวเลขทุกตัวใน EVALUATION มาจากการรันจริง **

clear; clc; close all;
c = physconst('LightSpeed');

%% ===== (0) CONFIG =====
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % root ของ repo
discRoot = fullfile(repoRoot, 'data');  % ที่วางข้อมูล GOTCHA
assert(exist(fullfile(discRoot,'GOTCHA-CP_Disc1','DATA'),'dir') > 0, ...
    'ไม่พบโฟลเดอร์ GOTCHA-CP_Disc1 ใน %s', discRoot);

% {name, id, cx, cy, heading, L, W, roofSpec(อ้างอิงภายนอก~)}
veh = {'ChevyMalibu',   'A',   9.9696,  -5.2239,   3.40747, 4.77, 1.74, 1.45;
       'ToyotaCamry',   'B',  20.6630, -18.7070, 182.80,    4.75, 1.74, 1.43;
       'FordTaurusWag', 'C',  12.4260, -18.2140, 185.04,    4.98, 1.86, 1.47;
       'CASEtractor',   'C1', -0.9590, -17.4780,  97.83,    4.73, 3.07, 2.60;
       'HysterForkLift','C2', 24.9640,  -6.4460, 273.04,    4.31, 1.50, 2.10;
       'NissanMaxima',  'D',  31.4230, -28.8740,   3.68,    4.79, 1.76, 1.44;
       'NissanSentra',  'E',  22.6850, -28.3020,   3.79,    4.45, 1.71, 1.39;
       'HyundaiSantaFe','F',  29.2360, -19.1770, 184.03,    4.45, 1.77, 1.68;
       'ChevyPrizm',    'J',  35.4360, -41.7210, 183.93,    4.41, 1.54, 1.37};
% ** G (SaturnIon 14.497,-26.924) และ H (VWJetta 4.494,-4.517) ถูกถอดออก:
%    ตรวจกับภาพ 2D ทั้งของเรา (fig_w8_display_afrl) และ AFRL reference แล้ว
%    ตำแหน่งตาม xls ว่างเปล่า — รถ 2 คันนี้ไม่ได้อยู่ตรงนั้นตอนเก็บข้อมูลจริง
%    (3D recon เห็นแค่ clutter พื้น / แนวขอบถนน สอดคล้องกัน) **
Nveh = size(veh,1);
% เลขรูปจาก montage ของ W9_ExtractCarPhotos.m (ลำดับตามตาราง veh)
photoIdx = [7 2 11 10 9 4 3 12 8];

polList    = {'HH','VV'};
subApWidth = 5;
subApStep  = 10;                 % 10 = 36 subap (default) | 5 = ครบ 72
NsigCFAR   = 30;                 % thr = Nsig*sigma_noise ราย look (CFAR-style;
relFloor   = 0.02;               %  ดูเหตุผลใน W9_GotchaSparse3D.m) + floor 2% peak
nIterFISTA = 200;                % iterate นานพอให้ L1 แยก ambiguity twin ออก
                                 % (ยืนยันแล้วกับ Camry: 30/200 ให้ภาพระดับเปเปอร์)
% gate แคบกว่า v แรก: บทเรียนจากรัน 11 คัน — เพื่อนบ้านห่าง ~5.5 m (เช่น
% Malibu ข้าง Jetta) รอด gate เดิม (slant 4.0 ~ ground 5.6) แล้ว WRAP เข้า
% กล่อง 10 m โผล่ฝั่งตรงข้าม สว่างกว่ารถเป้า -> G/H ล่ม (in-box 0%)
% slant 3.5 ~ ground 4.9 m + cross 4.5 -> ฆ่าทุกอย่างนอก ~5 m ก่อน wrap
gateRng    = 3.5;  gateCrs = 4.5;
topDB      = 40;                 % แสดง top-40dB แบบเปเปอร์
metDB      = 20;                 % วัดขนาด/roofline ที่ -20 dB (ชุดเข้ม)

N  = 80;  dr = 0.125;            % กล่อง 10 m เท่าเปเปอร์ (voxel หยาบกว่านิด
dk = 2*pi/(N*dr);                %  เพื่อความเร็ว x11 คัน; = 0.62832 เท่าเดิม)
axv = ((0:N-1) - N/2) * dr;

%% ===== (1) probe =====
passList = [];
for pn = 1:8
    dsc = 1 + (pn == 8);
    if exist(fullfile(discRoot, sprintf('GOTCHA-CP_Disc%d',dsc), 'DATA', ...
             sprintf('pass%d',pn), 'HH'), 'dir')
        passList(end+1) = pn; %#ok<SAGROW>
    end
end
S0 = load(fullfile(discRoot,'GOTCHA-CP_Disc1','DATA','pass1','HH', ...
                   'data_3dsar_pass1_az001_HH.mat'));
freq = double(S0.data.freq(:));
for v = 1:Nveh
    V(v).ctr3 = [veh{v,3}, veh{v,4}, 0];
    V(v).Vcs  = zeros(N,N,N);
end
nSub = numel(0:subApStep:359);
fprintf('W9b v2 | %d คัน | passes %s | subap %d x %d pol | grid %d^3 @ %.3f m\n', ...
    Nveh, mat2str(passList), nSub, numel(polList), N, dr);

%% ===== (2) MAIN LOOP =====
nLook = 0;  tAll = tic;
for aDeg = 0 : subApStep : 359
  for ip = 1:numel(polList)
    pol = polList{ip};
    passes = loadSubap(discRoot, passList, pol, aDeg, subApWidth);
    if numel(passes) < 4, continue; end
    for v = 1:Nveh
        [yg, mk] = gridSubap(passes, V(v).ctr3, N, dk, gateRng, gateCrs, c);
        xS = fistaKspace(yg, mk, nIterFISTA, NsigCFAR, relFloor);
        V(v).Vcs = max(V(v).Vcs, abs(xS));
    end
    nLook = nLook + 1;
    fprintf('  az %3d %s | %d passes | %.0f s (ETA %.0f min)\n', aDeg, pol, ...
        numel(passes), toc(tAll), ...
        toc(tAll)/nLook*(nSub*numel(polList)-nLook)/60);
  end
end
fprintf('รวม %d looks | total %.1f min\n', nLook, toc(tAll)/60);

%% ===== (3) METRICS ใน CAR FRAME =====
[XX, YY, ZZ] = ndgrid(axv, axv, axv);
for v = 1:Nveh
    Mn = V(v).Vcs / max(V(v).Vcs(:));
    selTop = Mn >= 10^(-topDB/20);
    nb = convn(double(selTop), ones(3,3,3), 'same') - selTop;
    sel = selTop & nb >= 1;                          % ตัด voxel โดดเดี่ยว
    thH = deg2rad(veh{v,5});
    Rr  = [cos(thH) -sin(thH); sin(thH) cos(thH)];
    uv  = Rr.' * [XX(sel).'; YY(sel).'];
    R(v).pu = uv(1,:).';  R(v).pv = uv(2,:).';  R(v).pz = ZZ(sel);
    R(v).pm = 20*log10(Mn(sel));
    s20 = R(v).pm >= -metDB;
    % วัดขนาดเฉพาะในเขตวิเคราะห์ +-3 m รอบรถ (กันเศษขอบกล่อง/เพื่อนบ้านปน)
    sEx = s20 & abs(R(v).pu) <= 3 & abs(R(v).pv) <= 3;
    R(v).L = 0; R(v).W = 0;
    if any(sEx)
        R(v).L = max(R(v).pv(sEx)) - min(R(v).pv(sEx));
        R(v).W = max(R(v).pu(sEx)) - min(R(v).pu(sEx));
    end
    inb = abs(R(v).pu) <= veh{v,7}/2+0.15 & abs(R(v).pv) <= veh{v,6}/2+0.15 ...
          & R(v).pz >= -0.5 & R(v).pz <= 3.0;
    R(v).nTop = nnz(sel);  R(v).prec = 100*mean(inb);
    zin = R(v).pz(inb & s20);
    R(v).zMed = NaN; R(v).zMax = NaN;
    if ~isempty(zin), R(v).zMed = median(zin); R(v).zMax = max(zin); end
    % roofline: max z ต่อ v-bin (เฉพาะ |u| ในกรอบ, จุด >= -metDB)
    vb = -3:0.25:3;  roof = nan(1, numel(vb)-1);
    m2 = s20 & abs(R(v).pu) <= veh{v,7}/2 + 0.2;
    for j = 1:numel(vb)-1
        zz = R(v).pz(m2 & R(v).pv >= vb(j) & R(v).pv < vb(j+1));
        if ~isempty(zz), roof(j) = max(zz); end
    end
    R(v).vb = vb(1:end-1) + 0.125;  R(v).roof = roof;
end

fprintf('\n========== EVALUATION — W9b v2 (k-space L1, %d looks) ==========\n', nLook);
for v = 1:Nveh
    fprintf('[%-14s %-2s] GT %.2fx%.2f roof~%.2f | @-%ddB L %.2f W %.2f | ', ...
        veh{v,1}, veh{v,2}, veh{v,6}, veh{v,7}, veh{v,8}, metDB, R(v).L, R(v).W);
    fprintf('vox@%ddB %4d (ในกรอบ %.0f%%) | z med/max %.2f / %.2f m\n', ...
        topDB, R(v).nTop, R(v).prec, R(v).zMed, R(v).zMax);
end
fprintf('\n[Discrimination] เทียบ (L, roofline, z max) ข้ามคัน: เก๋ง ~1.4 m,\n');
fprintf('SUV ~1.7 m, เครื่องจักรสูง >2 m ควรแยกได้จาก data จริง\n');

%% ===== (4) FIGURES =====
figDir = fullfile(repoRoot, 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end
hasPhoto = ~isempty(photoIdx);

for v = 1:Nveh
    fh = figure(10+v); clf;
    set(fh,'Name',sprintf('%s (%s)',veh{v,1},veh{v,2}),'Color','w', ...
        'Position',[50 50 (1150+280*hasPhoto) 380]);
    tiledlayout(1, 3+hasPhoto, 'TileSpacing','compact','Padding','compact');
    Lh = veh{v,6}/2;  Wh = veh{v,7}/2;
    box = [-Wh Wh Wh -Wh -Wh; -Lh -Lh Lh Lh -Lh];
    sz = 4 + 40*max(0, (R(v).pm + topDB)/topDB).^1.5;
    if hasPhoto
        nexttile;
        pf = fullfile(figDir, 'car_photos', sprintf('ppt_img_%02d.jpg', photoIdx(v)));
        if exist(pf,'file'), imshow(imread(pf)); else, axis off; end
        title('ภาพถ่าย AFRL');
    end
    axp = nexttile;                                  % 3D
    scatter3(R(v).pu, R(v).pv, R(v).pz, sz, R(v).pz, 'filled', 'MarkerFaceAlpha',0.75);
    hold on; plot3(box(1,:), box(2,:), zeros(1,5), '-', 'Color',[.6 .6 .6]); hold off;
    colormap(axp, turbo); clim([-0.3 2.4]);
    xlabel('u (m)'); ylabel('v (m)'); zlabel('z (m)');
    daspect([1 1 1]); view(-40, 22); grid on;
    xlim([-3 3]); ylim([-3 3]); zlim([-1 3]);
    title(sprintf('3D sparse (top %d dB)', topDB));
    axp = nexttile;                                  % side
    scatter(R(v).pv, R(v).pz, sz, R(v).pz, 'filled', 'MarkerFaceAlpha',0.75);
    hold on; yline(0,':', 'Color',[.5 .5 .5]);
    yline(veh{v,8},'--', sprintf('~%.2f spec',veh{v,8}), 'Color',[.4 .4 .4], ...
          'LabelHorizontalAlignment','left','FontSize',8);
    xline(-Lh,':','Color',[.5 .5 .5]); xline(Lh,':','Color',[.5 .5 .5]); hold off;
    colormap(axp, turbo); clim([-0.3 2.4]); grid on;
    xlabel('v (m)'); ylabel('z (m)'); xlim([-3 3]); ylim([-1 3]);
    title(sprintf('side: z med %.2f / max %.2f', R(v).zMed, R(v).zMax));
    axp = nexttile;                                  % top
    scatter(R(v).pu, R(v).pv, sz, R(v).pz, 'filled', 'MarkerFaceAlpha',0.75);
    hold on; plot(box(1,:), box(2,:), '--', 'Color',[.3 .8 .8]); hold off;
    colormap(axp, turbo); clim([-0.3 2.4]); grid on; axis equal;
    xlabel('u (m)'); ylabel('v (m)'); xlim([-3 3]); ylim([-3 3]);
    title(sprintf('top: %.1fx%.1f @-%ddB (จริง %.1fx%.1f)', ...
        R(v).L, R(v).W, metDB, veh{v,6}, veh{v,7}));
    sgtitle(sprintf('%s (%s) — k-space sparse L1, %d looks | GT %.2f x %.2f m', ...
        veh{v,1}, veh{v,2}, nLook, veh{v,6}, veh{v,7}));
    exportgraphics(fh, fullfile(figDir, sprintf('fig_w9cmp2_veh_%s.png', veh{v,2})), ...
        'Resolution', 150);
end

figure(3); clf; set(gcf,'Name','rooflines','Color','w','Position',[100 100 1150 460]);
cols = lines(Nveh);  hold on;
for v = 1:Nveh
    tall = veh{v,8} >= 2.0;
    stairs(R(v).vb, R(v).roof, 'Color', cols(v,:), 'LineWidth', 1.5+tall, ...
        'DisplayName', sprintf('%s (%s)', veh{v,1}, veh{v,2}));
end
hold off; grid on; xlim([-3 3]); ylim([-0.2 3.0]);
xlabel('v — ตามยาวรถ (m)'); ylabel(sprintf('top-z ที่ \\geq -%d dB (m)', metDB));
legend('Location','eastoutside');
title(sprintf('Fig 3 | W9b v2 — roofline ทุกคัน (k-space sparse, %d looks)', nLook));
exportgraphics(figure(3), fullfile(figDir,'fig_w9cmp2_rooflines.png'), 'Resolution',150);

save(fullfile(figDir,'w9cmp2_cache.mat'), 'V','R','axv','veh','nLook','-v7.3');
fprintf('\nFigures + cache saved to %s\n', figDir);


%%%% ===== LOCAL FUNCTIONS (engine เดียวกับ W9_GotchaSparse3D.m) =====

function passes = loadSubap(discRoot, passList, pol, aDeg, nFiles)
% Nf ไม่เท่ากันทุกไฟล์/pass -> เก็บ freq ราย pass (ดูหมายเหตุใน W9_GotchaSparse3D)
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

function [yg, mk] = gridSubap(passes, ctr3, N, dk, gateRng, gateCrs, c)
    KX=[]; KY=[]; KZ=[]; VL=[];
    for q = 1:numel(passes)
        fp = passes(q).fp;  ant = passes(q).ant;  r0u = passes(q).r0used;
        freq = passes(q).freq;
        Nf = numel(freq);  df = freq(2) - freq(1);
        P = size(fp, 2);
        Rc = sqrt(sum((ant - ctr3).^2, 2)).';
        fpc = fp .* exp(1j*(4*pi/c)*freq*(Rc - r0u));
        Nz = Nf*4;
        rc = fftshift(ifft(fpc, Nz, 1), 1);
        drax = ((0:Nz-1).' - Nz/2) * (c/(2*df)) / Nz;
        rc(abs(drax) > gateRng, :) = 0;
        sp = fft(ifftshift(rc, 1), [], 1);
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
        U = ant - ctr3;  U = U ./ sqrt(sum(U.^2, 2));
        kf = (4*pi/c) * freq;
        KX = [KX; reshape(kf*U(:,1).', [], 1)]; %#ok<AGROW>
        KY = [KY; reshape(kf*U(:,2).', [], 1)]; %#ok<AGROW>
        KZ = [KZ; reshape(kf*U(:,3).', [], 1)]; %#ok<AGROW>
        VL = [VL; fpc(:)];                      %#ok<AGROW>
    end
    K0 = [mean(KX), mean(KY), mean(KZ)];
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

function x = fistaKspace(yg, mk, nIter, Nsig, relFloor)
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
