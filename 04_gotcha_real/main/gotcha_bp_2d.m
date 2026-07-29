%% W8 — REAL DATA: GOTCHA parking lot (X-band) ด้วย circular BP pipeline จาก W7
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%
%  เป้าหมาย: เอา pipeline ที่พิสูจน์บนกล่องจำลอง (W7_CircularSAR) มารันบน
%  phase history จริง (AFRL GOTCHA, pass 1, az 0-4 deg, HH) → ภาพลานจอดรถจริง
%
%  ทำ 4 อย่าง:
%   (1) BP บนข้อมูลจริง (โครงเดียวกับ gotcha_BP_reference แต่ vectorized)
%   (2) WINDOWING: Taylor -35 dB บนแกนความถี่ (กด range sidelobe เหมือน W7 v2)
%       + เดโม azimuth window (กด cross-range sidelobe แลกกับ resolution)
%   (3) AUTOFOCUS: ใช้ data.af ที่มากับข้อมูล — จาก readme AFRL:
%       r_correct = แก้ r0 (m), ph_correct = แก้เฟสรายพัลส์ (rad)
%       ** ตรวจสอบแล้ว (prototype): ต้องใช้ "คู่กัน" คือ
%          r0 <- r0 + r_correct  และ  fp <- fp .* exp(+1j*ph_correct)
%          ใช้ผิดเครื่องหมาย/ใช้เฟสอย่างเดียว: ภาพพังทันที (entropy ~9.8)
%          บน subset 4 deg นี้ ผลของ af = "registration" เป็นหลัก:
%          เลื่อน CR ~0.7 m เข้าที่ (r_correct มี bulk offset ~0.28 m) แต่
%          focus แทบไม่เปลี่ยน (peak +0.08 dB, IRF กว้างเท่าเดิม) —
%          สมเหตุผล: aperture สั้น (~4 deg) motion error สะสมไม่ทัน
%          af จะสำคัญจริงเมื่อ mosaic หลาย subaperture รอบ 360 deg **
%   (4) SIM vs REAL: จำลอง point target ที่ตำแหน่ง corner reflector ด้วย
%       geometry + freq grid จริง (forward model เดียวกับ simPassFreq ของ W7)
%       → เทียบ IRF จริง/จำลอง/ทฤษฎี
%
%  โบนัส (ท้ายไฟล์): 3D tomography จาก 8 elevation passes จริง — data ครบแล้ว
%  (อ่าน az001-004 HH ตรงจาก data/GOTCHA-CP_Disc1 และ Disc2; pass8 มี 3 ไฟล์
%   เพราะ az004 ต้นฉบับอ่านไม่ได้ — ไม่กระทบ) elev 44.06-45.75 deg -> baseline
%   1.68 deg (~300 m) -> height res ~0.3-0.5 m (พิสูจน์บน data จริงแล้ว)
%
%  ระบบจริง: fc 9.6 GHz | BW 624 MHz | elev ~45.7 deg | aperture 4 deg (469 พัลส์)
%  ทฤษฎี: range res = c/2BW = 0.24 m | cross res ~ lambda/(2*dtheta) = 0.22 m
%  RUNTIME: main ~1-2 นาที + bonus tomography ~5-10 นาที (8 passes x 100k voxels)

clear; clc; close all;
c = physconst('LightSpeed');

%% ===== (0) LOAD: ต่อทุกไฟล์ az เป็น aperture เดียว =====
% ภาพ 2 มิติส่วนนี้ใช้ pass 1 · azimuth 4 องศา · HH เท่านั้น
dataDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'GOTCHA-CP_Disc1', 'DATA', 'pass1', 'HH');
files = dir(fullfile(dataDir, 'data_3dsar_pass1_az00[1-4]_HH.mat'));
assert(~isempty(files), 'ไม่พบข้อมูลใน %s', dataDir);

fp=[]; x=[]; y=[]; z=[]; r0=[]; rCor=[]; phCor=[]; freq=[];
for k = 1:numel(files)
    S = load(fullfile(dataDir, files(k).name));  d = S.data;
    fp  = [fp, double(d.fp)];   freq = double(d.freq(:));
    x   = [x, double(d.x)];  y = [y, double(d.y)];  z = [z, double(d.z)];
    r0  = [r0, double(d.r0)];
    rCor  = [rCor,  double(d.af.r_correct)];
    phCor = [phCor, double(d.af.ph_correct)];
end
[Nf, Np] = size(fp);  df = freq(2)-freq(1);  BW = Nf*df;  fcen = mean(freq);
lambda = c/fcen;
azSpan = atan2d(y(end),x(end)) - atan2d(y(1),x(1));       % ~4 deg
resRgTheory = c/(2*BW);
resCrTheory = lambda/(2*deg2rad(abs(azSpan)));
fprintf('GOTCHA pass1: %d freq x %d pulses | fc %.2f GHz | BW %.0f MHz | az span %.2f deg\n', ...
    Nf, Np, fcen/1e9, BW/1e6, azSpan);
fprintf('Theory: range res %.3f m | cross res %.3f m | elev %.2f deg\n', ...
    resRgTheory, resCrTheory, atan2d(z(1), hypot(x(1),y(1))));

%% ===== (1)-(3) เตรียม 3 เคสเทียบ: raw / +Taylor / +Taylor+autofocus =====
wFreq = taylorwin(Nf, 4, -35);                 % (2) window แบบเดียวกับ W7 v2
fpAF  = fp .* exp(1j*phCor);                   % (3) เฟสรายพัลส์ (broadcast ทุก freq)
r0AF  = r0 + rCor;                             % (3) แก้ระยะอ้างอิง

Zp = 8;  Nzp = Nf*Zp;
drAxis = ((0:Nzp-1).' - Nzp/2) * (c/(2*df)) / Nzp;   % differential range axis
% ** สำคัญ: phase-alignment ของ BP ต้องใช้ freq(1) (min freq) ไม่ใช่ fcen **
% เพราะเฟสตกค้างหลัง IFFT คือ exp(-j*4*pi*freq(1)*dr/c) — convention เดียวกับ
% bpBasic ของ AFRL (minF) | ใช้ fcen (แบบ gotcha_BP_reference เดิม) ภาพยังโฟกัส
% แต่ตำแหน่งเพี้ยนตาม geometry ~0.7 m (พิสูจน์: sim point เลื่อน 0.68 m,
% พอเปลี่ยนเป็น freq(1) แล้วลงตำแหน่งเป๊ะ 0.00 m)
kc = 4*pi*freq(1)/c;

% range-compressed profiles (Nzp x Np) ต่อเคส
rcRaw = fftshift(ifft(fp,             Nzp, 1), 1);
rcWin = fftshift(ifft(fp  .* wFreq,   Nzp, 1), 1);
rcAF  = fftshift(ifft(fpAF.* wFreq,   Nzp, 1), 1);

%% ===== BP ทั้งฉาก 3 เคส (grid พื้น z=0, 100x100 m @ 0.25 m) =====
N = 400;  Wsc = 50;  g = linspace(-Wsc, Wsc, N);
[GX, GY] = meshgrid(g, g);
fprintf('\nBackprojection %dx%d px x %d pulses x 3 cases ...\n', N, N, Np);
tBP = tic;
imgRaw = bpGrid(rcRaw, r0,   x,y,z, GX,GY, drAxis, kc);
imgWin = bpGrid(rcWin, r0,   x,y,z, GX,GY, drAxis, kc);
imgAF  = bpGrid(rcAF,  r0AF, x,y,z, GX,GY, drAxis, kc);
fprintf('done (%.1f s)\n', toc(tBP));

[entRaw, pmRaw] = imgMetrics(imgRaw);
[entWin, pmWin] = imgMetrics(imgWin);
[entAF,  pmAF ] = imgMetrics(imgAF);

%% ===== FIGURE 1: raw -> windowed -> autofocused =====
figure(1); set(gcf,'Name','Fig1: window+autofocus','Color','k','Position',[60 60 1380 480]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
panels = {imgRaw, sprintf('raw BP\\newlineentropy %.2f | peak/mean %.1f dB', entRaw, pmRaw);
          imgWin, sprintf('+ Taylor -35 dB\\newlineentropy %.2f | peak/mean %.1f dB', entWin, pmWin);
          imgAF,  sprintf('+ autofocus (af)\\newlineentropy %.2f | peak/mean %.1f dB', entAF, pmAF)};
for k = 1:3
    ax = nexttile;
    M = abs(panels{k,1});  M = M/max(M(:));
    imagesc(g, g, 20*log10(M+1e-6)); set(ax,'YDir','normal'); clim([-35 0]); colormap(ax,'gray');
    axis image; title(panels{k,2},'Color','w');
    xlabel('x (m)','Color','w'); ylabel('y (m)','Color','w');
    set(ax,'Color','k','XColor','w','YColor','w');
end
cb = colorbar; cb.Color='w'; ylabel(cb,'dB','Color','w'); cb.Layout.Tile='east';
sgtitle('Fig 1 | W8 — REAL GOTCHA parking lot: ผลของ windowing + autofocus','Color','w');

%% ===== หา corner reflector (จุดสว่างสุดในภาพ af) =====
Mbest = abs(imgAF);
[~, pk] = max(Mbest(:));  [iy, ix] = ind2sub(size(Mbest), pk);
px = g(ix);  py = g(iy);
fprintf('\nCorner reflector (brightest point): x=%.2f m, y=%.2f m\n', px, py);

%% ===== (4) SIM vs REAL: point target ที่ตำแหน่ง CR ด้วย geometry จริง =====
% forward model เดียวกับ simPassFreq (W7): fp(f,p) = exp(-j*4*pi*f*(R-r0)/c)
Rsim  = sqrt((px-x).^2 + (py-y).^2 + z.^2);            % 1 x Np (จุดบนพื้น z=0)
fpSim = exp(-1j*4*pi/c * freq * (Rsim - r0AF));        % Nf x Np (outer product)
rcSim = fftshift(ifft(fpSim .* wFreq, Nzp, 1), 1);

% zoom patch +-2 m @ 2 cm รอบ CR — real / sim / (เดโม azimuth window)
gzx = px + linspace(-2, 2, 201);  gzy = py + linspace(-2, 2, 201);
[ZX, ZY] = meshgrid(gzx, gzy);
zmReal = bpGrid(rcAF,  r0AF, x,y,z, ZX,ZY, drAxis, kc);
zmSim  = bpGrid(rcSim, r0AF, x,y,z, ZX,ZY, drAxis, kc);
wAz    = hamming(Np).';                                % (2) เดโม azimuth taper
zmAzW  = bpGrid(rcAF .* wAz, r0AF, x,y,z, ZX,ZY, drAxis, kc);

zmReal = abs(zmReal)/max(abs(zmReal(:)));
zmSim  = abs(zmSim) /max(abs(zmSim(:)));
zmAzW  = abs(zmAzW) /max(abs(zmAzW(:)));

% -3 dB widths (แกน x ~ range เพราะ radar อยู่ทาง az~0-4 deg; แกน y ~ cross)
[wRgReal, wCrReal] = irfWidths(zmReal, gzx, gzy);
[wRgSim,  wCrSim ] = irfWidths(zmSim,  gzx, gzy);
[wRgAzW,  wCrAzW ] = irfWidths(zmAzW,  gzx, gzy);

%% ===== FIGURE 2: zoom CR — real vs sim vs azimuth-windowed =====
figure(2); set(gcf,'Name','Fig2: sim vs real IRF','Color','w','Position',[80 80 1380 460]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
zp = {zmReal, sprintf('REAL corner reflector\\newlinerange %.2f | cross %.2f m', wRgReal, wCrReal);
      zmSim,  sprintf('SIM point (geometry จริง)\\newlinerange %.2f | cross %.2f m', wRgSim, wCrSim);
      zmAzW,  sprintf('REAL + azimuth Hamming\\newlinerange %.2f | cross %.2f m', wRgAzW, wCrAzW)};
for k = 1:3
    ax = nexttile;
    imagesc(gzx, gzy, 20*log10(zp{k,1}+1e-9)); set(ax,'YDir','normal');
    clim([-35 0]); colormap(ax,'jet'); axis image;
    xlabel('x (m)'); ylabel('y (m)'); title(zp{k,2});
end
cb = colorbar; ylabel(cb,'dB'); cb.Layout.Tile='east';
sgtitle('Fig 2 | W8 — IRF จริง vs จำลอง (pipeline เดียวกัน 100%)');

%% ===== FIGURE 3: range/cross cuts ทับกัน =====
[~, jy] = max(max(zmReal,[],2));  [~, jx] = max(max(zmReal,[],1));
figure(3); set(gcf,'Name','Fig3: IRF cuts','Color','w','Position',[100 100 980 420]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(gzx-px, 20*log10(zmReal(jy,:)+1e-9),'b','LineWidth',1.4); hold on;
plot(gzx-px, 20*log10(zmSim(jy,:) +1e-9),'r--','LineWidth',1.4); hold off;
grid on; ylim([-40 0]); xlim([-2 2]);
xlabel('\Delta x — range (m)'); ylabel('dB'); legend('real','sim');
title(sprintf('range cut: real %.2f / sim %.2f / theory %.2f m', ...
    wRgReal, wRgSim, resRgTheory*1.3));                 % x1.3 = Taylor broadening
nexttile;
plot(gzy-py, 20*log10(zmReal(:,jx)+1e-9),'b','LineWidth',1.4); hold on;
plot(gzy-py, 20*log10(zmSim(:,jx) +1e-9),'r--','LineWidth',1.4);
plot(gzy-py, 20*log10(zmAzW(:,jx) +1e-9),'Color',[0 .6 0],'LineWidth',1.2); hold off;
grid on; ylim([-40 0]); xlim([-2 2]);
xlabel('\Delta y — cross-range (m)'); ylabel('dB');
legend('real (uniform az)','sim','real + az Hamming');
title(sprintf('cross cut: real %.2f / sim %.2f / theory %.2f m', ...
    wCrReal, wCrSim, resCrTheory));
sgtitle('Fig 3 | W8 — IRF cuts ผ่านจุด corner reflector');

%% ===== EVALUATION =====
fprintf('\n========== EVALUATION — W8: real GOTCHA vs sim ==========\n');
fprintf('\n[Focus quality]           entropy   peak/mean\n');
fprintf('  raw BP                 : %6.2f   %6.1f dB\n', entRaw, pmRaw);
fprintf('  + Taylor window        : %6.2f   %6.1f dB\n', entWin, pmWin);
fprintf('  + autofocus (af)       : %6.2f   %6.1f dB\n', entAF, pmAF);

% af บน subset นี้ = registration เป็นหลัก -> วัด "การเลื่อนตำแหน่ง CR"
MnoAF = abs(imgWin);
[~, pk2] = max(MnoAF(:));  [iy2, ix2] = ind2sub(size(MnoAF), pk2);
fprintf('\n[Autofocus = registration] CR เลื่อน (%.2f, %.2f) -> (%.2f, %.2f) m = %.2f m\n', ...
    g(ix2), g(iy2), px, py, hypot(g(ix2)-px, g(iy2)-py));
fprintf('  (entropy ทั้งฉากแทบไม่เปลี่ยน — aperture 4 deg สั้นเกินกว่า motion error จะสะสม)\n');
fprintf('\n[-3 dB IRF @ corner reflector]   range      cross\n');
fprintf('  theory (Taylor x1.3 / uniform): %5.2f m   %5.2f m\n', resRgTheory*1.3, resCrTheory);
fprintf('  sim point (geometry จริง)     : %5.2f m   %5.2f m\n', wRgSim, wCrSim);
fprintf('  REAL corner reflector         : %5.2f m   %5.2f m\n', wRgReal, wCrReal);
fprintf('  REAL + azimuth Hamming        : %5.2f m   %5.2f m  (sidelobe ลด, กว้างขึ้น)\n', wRgAzW, wCrAzW);
fprintf('\n[Note] sim กับ real ใช้ pipeline เดียวกันทุกขั้น — ต่างกันแค่ fp\n');
fprintf('       kc=4*pi*freq(1)/c ทำให้ sim ลงตำแหน่งที่วางเป๊ะ -> เทียบ real ได้ตรงจุด\n');
fprintf('       ถ้า real กว้าง/มี shoulder กว่า sim = clutter/multipath จริงรอบ CR\n');

%% ===== FIGURE 5: display stretch + เทียบภาพอ้างอิง AFRL =====
% "ภาพมืด" ไม่ใช่ปัญหาของ data — เป็นเรื่อง DISPLAY: เรา normalize 0 dB ที่
% CR peak (trihedral สว่างมาก) พอ floor -35 dB รถทั้งหมด (~-30..-45 dB) เลยจมดำ
% SAR ต้อง map dynamic range ก่อนดูเสมอ — เลื่อน clim ลงไปครอบช่วง clutter พอ
% ภาพอ้างอิง: AFRL ทำภาพฉากเดียวกันไว้ใน DOCUMENTATION/Challenge_Pictures_Images.ppt
% (ถ้ามี ให้วางเป็น data/AFRL_reference_2Dimage.png — แกน +-50 m เท่ากัน
%  จุด CR ของ AFRL อยู่ (-16, 21) ตรงกับที่เราวัด (-15.91, 20.93) ภายใน ~0.1 m)
figure(5); set(gcf,'Name','Fig5: display + AFRL reference','Color','w','Position',[70 70 1150 470]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
MA = abs(imgAF)/max(abs(imgAF(:)));
imagesc(g, g, 20*log10(MA+1e-6)); set(gca,'YDir','normal'); axis image;
clim([-50 -25]); colormap(gca,'jet');           % << เลื่อนช่วงลงมาที่ clutter
cb5 = colorbar; ylabel(cb5,'dB (0 = CR peak)');
xlabel('x (m)'); ylabel('y (m)');
title('data เดิม เปลี่ยนแค่ display: clim [-50,-25] dB');
nexttile;
% รูปอ้างอิงของ AFRL แกะจากเอกสารของ AFRL เอง ไม่ได้แนบใน repo — ข้ามได้ถ้าไม่มี
refFile = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'AFRL_reference_2Dimage.png');
if exist(refFile,'file')
    imshow(imread(refFile)); title('AFRL reference (Challenge\_Pictures\_Images.ppt)');
else
    axis off; title('ไม่พบ AFRL\_reference\_2Dimage.png');
end
sgtitle('Fig 5 | W8 — display stretch: รถโผล่ครบเหมือนภาพอ้างอิง AFRL');

%% ===== FIGURE 9: GROUND TRUTH OVERLAY + DETECTION CHECK =====
% รถ staged 11 คัน + เสาไฟ 7 ต้น จาก DOCUMENTATION/"Gotcha Spotlight Target
% Locations.xls" (พิกัด center + heading + ขนาดจริง) — ตอบคำถาม "เจอครบไหม/เกินไหม"
% หมายเหตุ: แถวรถฝั่งซ้าย (x<0) เป็นรถทั่วไปในลาน ไม่อยู่ในรายการ ground truth
% vehGT: {name, id, cx, cy, heading(deg), L(m), W(m)}
vehGT = {'ChevyMalibu','A',9.970,-5.224,3.41,4.77,1.74; 'ToyotaCamry','B',20.663,-18.707,182.80,4.75,1.74;
         'FordTaurusWag','C',12.426,-18.214,185.04,4.98,1.86; 'CASEtractor','C1',-0.959,-17.478,97.83,4.73,3.07;
         'HysterForkLift','C2',24.964,-6.446,273.04,4.31,1.50; 'NissanMaxima','D',31.423,-28.874,3.68,4.79,1.76;
         'NissanSentra','E',22.685,-28.302,3.79,4.45,1.71; 'HyundaiSantaFe','F',29.236,-19.177,184.03,4.45,1.77;
         'SaturnIon','G',14.497,-26.924,5.92,4.63,1.73; 'VWJetta','H',4.494,-4.517,4.24,4.36,1.66;
         'ChevyPrizm','J',35.436,-41.721,183.93,4.41,1.54};
poleGT = [22.496 -68.803; -20.669 -65.721; -40.940 -28.698; 10.607 -44.584;
          -4.368 -21.937; 24.160 -2.202; -24.150 1.190];
MAdB = 20*log10(abs(imgAF)/max(abs(imgAF(:))) + eps);
detThr = -42;                                   % dB rel CR peak (clutter ~ -45..-50)
fprintf('\n[Ground truth detection — staged vehicles]\n');
fprintf('  %-16s %-3s  max dB   status\n', 'vehicle', 'ID');
figure(9); set(gcf,'Name','Fig9: ground truth overlay','Color','w','Position',[100 60 760 700]);
imagesc(g, g, MAdB); set(gca,'YDir','normal'); axis image;
clim([-50 -25]); colormap(gray); hold on;
nDet = 0;
for v = 1:size(vehGT,1)
    [cx,cy,hd,L,W] = deal(vehGT{v,3},vehGT{v,4},vehGT{v,5},vehGT{v,6},vehGT{v,7});
    thv = deg2rad(hd);  Rv = [cos(thv) -sin(thv); sin(thv) cos(thv)];
    % max dB ในกรอบรถ (+0.8 m margin) — หมุนพิกัด pixel เข้า frame รถ
    [GXm, GYm] = meshgrid(g, g);
    q = Rv.' * [GXm(:).'-cx; GYm(:).'-cy];
    inBox = abs(q(1,:)) <= W/2+0.8 & abs(q(2,:)) <= L/2+0.8;
    mV = max(MAdB(inBox));
    ok = mV >= detThr;  nDet = nDet + ok;
    stTxt = {'MISS', 'DETECTED'};
    fprintf('  %-16s %-3s  %6.1f   %s\n', vehGT{v,1}, vehGT{v,2}, mV, stTxt{ok+1});
    crnV = Rv*[-W/2 W/2 W/2 -W/2 -W/2; -L/2 -L/2 L/2 L/2 -L/2] + [cx; cy];
    plot(crnV(1,:), crnV(2,:), 'c-', 'LineWidth', 1.5);
    text(cx+1.2, cy+1.2, vehGT{v,2}, 'Color','y', 'FontSize',11, 'FontWeight','bold');
end
plot(poleGT(:,1), poleGT(:,2), 'rx', 'MarkerSize', 9, 'LineWidth', 2);
plot(px, py, 'c+', 'MarkerSize', 14, 'LineWidth', 2); text(px+1.2, py+1.2, 'CAL', 'Color','c');
hold off; xlabel('x (m)'); ylabel('y (m)');
xlim([-50 50]); ylim([-50 50]);
title(sprintf('Fig 9 | ground truth overlay — detected %d/%d staged vehicles', nDet, size(vehGT,1)));
fprintf('  -> detected %d/%d | จุดสว่างอื่น = รถทั่วไปในลาน (ไม่มี GT), เสาไฟ, glint ราย aspect\n', ...
    nDet, size(vehGT,1));

%% ===== SAVE FIGURES =====
figDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end
exportgraphics(figure(1), fullfile(figDir,'fig_w8_window_autofocus.png'),'Resolution',150);
exportgraphics(figure(2), fullfile(figDir,'fig_w8_sim_vs_real.png'),     'Resolution',150);
exportgraphics(figure(3), fullfile(figDir,'fig_w8_irf_cuts.png'),        'Resolution',150);
exportgraphics(figure(5), fullfile(figDir,'fig_w8_display_afrl.png'),    'Resolution',150);
exportgraphics(figure(9), fullfile(figDir,'fig_w8_gt_overlay.png'),      'Resolution',150);
fprintf('\nFigures saved to %s\n', figDir);

%% ===== BONUS: 3D TOMOGRAPHY ถ้ามีหลาย pass =====
% GOTCHA เต็มมี 8 elevation passes (ต้องลงทะเบียน AFRL SDMS: sdms.afrl.af.mil)
% วางไฟล์ data_3dsar_passN_az*_HH.mat เพิ่มในโฟลเดอร์เดียวกันแล้วรันใหม่ —
% ส่วนนี้จะทำงานเองเมื่อพบ >= 2 passes
passList = [];
for pn = 1:8
    if ~isempty(dir(fullfile(dataDir, sprintf('data_3dsar_pass%d_az*_HH.mat', pn))))
        passList(end+1) = pn; %#ok<SAGROW>
    end
end
if numel(passList) < 2
    fprintf('\n[BONUS] พบ %d pass (ต้องมี >= 2 จึงทำ 3D tomography ได้) — ข้าม\n', numel(passList));
    fprintf('        ชุดเต็ม 8 passes: ลงทะเบียนที่ sdms.afrl.af.mil (collection GOTCHA)\n');
else
    % ** ตรวจสอบบน data จริง 8 passes แล้ว (prototype):
    %    - single pass = layover ridge เอียงตาม LOS (ไม่มี height resolution)
    %    - เฟสข้าม pass ที่ CR เป็น ramp เรียบ (af ทำให้แต่ละ pass self-coherent)
    %      แต่ต้อง CALIBRATE เฟสข้าม pass ก่อน coherent sum: เอาเฟสของ voxel
    %      ที่ CR (peak ของ pass แรก) ลบออกราย pass -> z-profile สะอาด
    %    - ได้ height res วัดจริง ~0.3 m (ทฤษฎี ~0.5 m จาก baseline 1.68 deg) **
    fprintf('\n[BONUS] พบ %d passes: %s — 3D tomography รอบ CR\n', ...
        numel(passList), mat2str(passList));
    vgx = px + (-2 : 0.1 : 2);  vgy = py + (-2 : 0.1 : 2);  vgz = -3 : 0.15 : 6;
    [VX, VY, VZ] = meshgrid(vgx, vgy, vgz);      % (row=y, col=x, page=z)
    NpassB = numel(passList);
    volP  = cell(1, NpassB);  elevs = zeros(1, NpassB);
    tTomo = tic;
    for q = 1:NpassB
        fl = dir(fullfile(dataDir, sprintf('data_3dsar_pass%d_az*_HH.mat', passList(q))));
        fpP=[]; xP=[]; yP=[]; zP=[]; r0P=[]; freqP=[];
        for k = 1:numel(fl)
            S = load(fullfile(dataDir, fl(k).name));  d = S.data;
            fpk = double(d.fp) .* exp(1j*double(d.af.ph_correct));
            fpP = [fpP, fpk]; xP=[xP,double(d.x)]; yP=[yP,double(d.y)]; zP=[zP,double(d.z)];
            r0P = [r0P, double(d.r0)+double(d.af.r_correct)];
            freqP = double(d.freq(:));
        end
        % ** แต่ละ pass มี Nf/df ไม่เท่ากัน (424-434 / 1.44-1.47 MHz) —
        %    ต้องสร้าง window + range axis ใหม่ราย pass (f0 เท่ากันทุก pass
        %    -> kc = 4*pi*freq(1)/c ใช้ตัวเดิมได้) **
        NfP = size(fpP,1);  dfP = freqP(2)-freqP(1);  NzpP = NfP*Zp;
        wFreqP  = taylorwin(NfP, 4, -35);
        drAxisP = ((0:NzpP-1).' - NzpP/2) * (c/(2*dfP)) / NzpP;
        rcP = fftshift(ifft(fpP .* wFreqP, NzpP, 1), 1);
        volP{q} = bpVol(rcP, r0P, xP,yP,zP, VX,VY,VZ, drAxisP, kc);
        elevs(q) = atan2d(mean(zP), mean(hypot(xP,yP)));
        fprintf('  pass %d done (Nf=%d, elev %.3f deg | %.0f s)\n', ...
            passList(q), NfP, elevs(q), toc(tTomo));
    end
    % inter-pass phase calibration ที่ CR voxel (peak ของ pass แรก)
    [~, iCR] = max(abs(volP{1}(:)));
    volCal = zeros(size(VX));  volInc = zeros(size(VX));
    for q = 1:NpassB
        phq   = angle(volP{q}(iCR));
        volCal = volCal + volP{q} * exp(-1j*phq);       % calibrated coherent
        volInc = volInc + abs(volP{q}).^2;              % incoherent (control)
    end
    Bel = deg2rad(max(elevs)-min(elevs)) * mean(r0AF);
    fprintf('  elev baseline %.2f deg (~%.0f m) -> height res theory ~%.2f m\n', ...
        max(elevs)-min(elevs), Bel, lambda*mean(r0AF)/(2*Bel));

    % วัด -3dB height width ผ่าน CR
    [iyC, ixC, izC] = ind2sub(size(VX), iCR);
    zProf = squeeze(abs(volCal(iyC, ixC, :)));
    wZ = w3(vgz, zProf);
    fprintf('  measured -3dB height width @ CR = %.2f m (single pass: resolve ไม่ได้เลย)\n', wZ);

    % Fig 4: x-z slice single pass (layover) vs 8-pass calibrated + z-profile
    figure(4); set(gcf,'Name','Fig4: BONUS 3D tomography','Color','w','Position',[80 80 1380 430]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    sl1 = squeeze(abs(volP{1}(iyC, :, :))).';   % z x x
    slC = squeeze(abs(volCal(iyC, :, :))).';
    nexttile;
    imagesc(vgx-px, vgz, 20*log10(sl1/max(sl1(:))+1e-9)); set(gca,'YDir','normal');
    clim([-25 0]); colormap(gca,'jet'); xlabel('\Delta x (m)'); ylabel('z (m)');
    title(sprintf('single pass %d — layover ridge', passList(1)));
    nexttile;
    imagesc(vgx-px, vgz, 20*log10(slC/max(slC(:))+1e-9)); set(gca,'YDir','normal');
    clim([-25 0]); colormap(gca,'jet'); xlabel('\Delta x (m)'); ylabel('z (m)');
    title(sprintf('%d passes calibrated — height resolved', NpassB));
    nexttile;
    plot(vgz, 20*log10(zProf/max(zProf)+1e-9), 'b', 'LineWidth', 1.5); hold on;
    z1 = squeeze(abs(volP{1}(iyC, ixC, :)));
    plot(vgz, 20*log10(z1/max(z1)+1e-9), 'r--', 'LineWidth', 1.2); hold off;
    grid on; ylim([-30 1]); xlabel('z (m)'); ylabel('dB');
    legend(sprintf('%d passes (%.2f m)', NpassB, wZ), 'single pass');
    title('vertical profile @ CR');
    sgtitle(sprintf('Fig 4 | BONUS — 3D tomography จาก data จริง %d passes (kc=freq(1), CR-calibrated)', NpassB));
    exportgraphics(figure(4), fullfile(figDir,'fig_w8_bonus_tomo3d.png'),'Resolution',150);
    fprintf('  Fig 4 saved. Tomography total %.0f s\n', toc(tTomo));
end

%% ===== BONUS 4: FULL-CIRCLE MULTI-ASPECT 3D (วิธีที่ทำให้ "เห็นเป็นรถ") =====
% ทำไม BP 3D แถวรถ (4 deg) ยังไม่เป็นรูปรถทั้งที่ pipeline ถูก: เราใช้ azimuth 4 deg จาก 360
% (~1% ของข้อมูล) — รถเป็น ANISOTROPIC scatterer แต่ละมุมเห็นคนละ glint
% ภาพสวยของ AFRL ใช้ complete circular aperture -> ต้องทำแบบ W7_CircularSAR:
%   subaperture 4 deg ทุก ๆ aspectStep รอบวง x 8 passes
%   coherent ภายใน aspect (calibrate เฟสข้าม pass ที่ tophat ซึ่งเห็นได้ 360 deg)
%   แล้ว MULTILOOK ข้าม aspect (เฉลี่ย intensity — ห้าม normalize ราย aspect)
% ** radar config แก้ไม่ได้ (data บันทึกมาแล้ว) — ตัวแปรที่คุมได้คือ aspect
%    coverage, ตำแหน่ง patch (ใช้ ground truth จาก xls), และ recon algorithm **
% หมายเหตุ: อ่านไฟล์ตรงจาก GOTCHA-CP_Disc1/2 | aspectStep=15 -> ~24 aspects
% (ลดเป็น 4 เพื่อใช้ครบทุก subaperture: ช้าลง ~4 เท่า แต่ดีขึ้นอีก)
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % root ของ repo
discRoot = fullfile(repoRoot, 'data');  % ที่วางข้อมูล GOTCHA
if ~exist(fullfile(discRoot,'GOTCHA-CP_Disc1','DATA'), 'dir')
    fprintf('\n[BONUS4] ไม่พบโฟลเดอร์ GOTCHA-CP_Disc1 — ข้าม full-circle 3D\n');
else
    fprintf('\n[BONUS4] full-circle multi-aspect 3D รอบรถ\n');
    % ** เป้า: Chevy Malibu (ID "A") — ground truth จาก DOCUMENTATION/
    %    "Gotcha Spotlight Target Locations.xls": center (9.9696, -5.2239),
    %    heading 3.41 deg, ขนาดจริง L 4.77 x W 1.74 m
    %    ตรวจแล้ว (full circle 90 aspects บนไฟล์ ~53%): พลังงานอยู่ในกรอบจริง
    %    เห็นเส้น glint ข้างรถ + footprint -10 dB ~3.4 x 3.2 m
    %    (L สั้นกว่าจริง: หัว-ท้ายสะท้อนอ่อน | W เกิน: sidelobe) **
    carCtr = [9.9696, -5.2239];   carHd = 3.40747;   % deg
    carL   = 4.77;                carW  = 1.74;      % m (ground truth)
    aspectStep = 4;   subApFiles = 4;   % 4 = เต็มวง 90 aspects (~2-3 นาที)
                                        % เปลี่ยนเป็น 15 ถ้าอยากได้เร็ว ๆ
    vgx4 = carCtr(1) + (-2 : 0.15 : 2);
    vgy4 = carCtr(2) + (-3 : 0.15 : 3);
    vgz4 = -1 : 0.15 : 2.6;
    [VX4, VY4, VZ4] = meshgrid(vgx4, vgy4, vgz4);
    [CX, CY, CZ]    = meshgrid(px, py, 0.30);          % tophat voxel (calibration)
    Icar = zeros(size(VX4));  nAsp = 0;  tB4 = tic;
    for aDeg = 0 : aspectStep : 359
        volA = zeros(size(VX4));  usedP = 0;
        for pn = passList
            dsc = 1 + (pn == 8);                        % pass8 อยู่ Disc2
            fpP=[]; xP=[]; yP=[]; zP=[]; r0P=[]; freqP=[];
            for k = 1:subApFiles
                azi = mod(aDeg + k - 1, 360) + 1;
                fn  = fullfile(discRoot, sprintf('GOTCHA-CP_Disc%d',dsc), 'DATA', ...
                      sprintf('pass%d',pn), 'HH', ...
                      sprintf('data_3dsar_pass%d_az%03d_HH.mat', pn, azi));
                if ~exist(fn,'file'), continue; end
                S = load(fn);  d = S.data;
                fpP = [fpP, double(d.fp).*exp(1j*double(d.af.ph_correct))];
                xP=[xP,double(d.x)]; yP=[yP,double(d.y)]; zP=[zP,double(d.z)];
                r0P = [r0P, double(d.r0)+double(d.af.r_correct)];
                freqP = double(d.freq(:));
            end
            if isempty(fpP), continue; end
            NfP = size(fpP,1);  dfP = freqP(2)-freqP(1);  NzpP = NfP*Zp;
            wFq = taylorwin(NfP, 4, -35);
            drA = ((0:NzpP-1).' - NzpP/2) * (c/(2*dfP)) / NzpP;
            rcP = fftshift(ifft(fpP .* wFq, NzpP, 1), 1);
            vP  = bpVol(rcP, r0P, xP,yP,zP, VX4,VY4,VZ4, drA, kc);
            cP  = bpVol(rcP, r0P, xP,yP,zP, CX,CY,CZ,    drA, kc);
            volA = volA + vP * exp(-1j*angle(cP));      % coherent ข้าม pass
            usedP = usedP + 1;
        end
        if usedP == 0, continue; end
        Icar = Icar + abs(volA).^2;  nAsp = nAsp + 1;   % multilook ข้าม aspect
        fprintf('  aspect az %3d deg | %d passes | %.0f s\n', aDeg, usedP, toc(tB4));
    end
    Mcar = sqrt(Icar / max(nAsp,1));  Mcar = Mcar / max(Mcar(:));
    fprintf('  รวม %d aspects (%.0f s)\n', nAsp, toc(tB4));

    % กรอบรถจริง (หมุนตาม heading)
    thH  = deg2rad(carHd);
    Rrot = [cos(thH) -sin(thH); sin(thH) cos(thH)];
    crn  = Rrot * [ -carW/2 carW/2 carW/2 -carW/2 -carW/2;
                    -carL/2 -carL/2 carL/2 carL/2 -carL/2 ] + carCtr(:);

    % วัด footprint -10 dB เทียบ ground truth
    mip4 = max(Mcar, [], 3);
    [iyF, ixF] = find(20*log10(mip4+eps) >= -10);
    fprintf('  [Ground truth check] footprint -10dB: L(y) %.2f m | W(x) %.2f m', ...
        vgy4(max(iyF))-vgy4(min(iyF)), vgx4(max(ixF))-vgx4(min(ixF)));
    fprintf('  (จริง %.2f x %.2f)\n', carL, carW);

    figure(8); set(gcf,'Name','Fig8: BONUS4 full-circle car 3D','Color','w','Position',[130 90 1250 520]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile;
    imagesc(vgx4, vgy4, 20*log10(mip4+1e-9)); set(gca,'YDir','normal'); axis image;
    clim([-14 0]); colormap(gca,'gray'); hold on;
    plot(crn(1,:), crn(2,:), 'c--', 'LineWidth', 2); hold off;
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('Chevy Malibu (A): %d aspects multilook + กรอบจริง', nAsp));
    nexttile;
    thr4 = -9;                                         % multilook floor สูง -> threshold ตื้น
    idx4 = find(20*log10(Mcar+eps) >= thr4 & VZ4 >= -0.5);
    scatter3(VX4(idx4), VY4(idx4), VZ4(idx4), 14, VZ4(idx4), 'filled', 'MarkerFaceAlpha', 0.55);
    hold on; plot3(crn(1,:), crn(2,:), zeros(1,5), 'c--', 'LineWidth', 2); hold off;
    colormap(gca, turbo); clim([-0.5 2.0]);
    cb8 = colorbar; ylabel(cb8,'z (m)');
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
    daspect([1 1 1]); grid on; box on; view(-50, 30); rotate3d on;
    title(sprintf('point cloud \\geq %d dB + กรอบจริง', thr4));
    sgtitle('Fig 8 | BONUS 4 — full-circle 3D บนรถที่มี ground truth (Malibu 4.77 x 1.74 m)');
    exportgraphics(figure(8), fullfile(figDir,'fig_w8_bonus_fullcircle.png'),'Resolution',150);
end


%%%% ===== LOCAL FUNCTIONS =====

function img = bpGrid(rc, r0v, x, y, z, GX, GY, drAxis, kc)
% 2D BP บนพื้น z=0, แกน differential range (vectorized ต่อ pulse)
    Nzp = numel(drAxis); ddr = drAxis(2)-drAxis(1); d1 = drAxis(1);
    img = zeros(size(GX));
    for p = 1:numel(r0v)
        dr = sqrt((GX-x(p)).^2 + (GY-y(p)).^2 + z(p)^2) - r0v(p);
        fi = (dr - d1)/ddr + 1;
        i0 = floor(fi);  w = fi - i0;
        v  = (i0 >= 1) & (i0 < Nzp);
        Sa = rc(:,p);  iv = zeros(size(GX));
        iv(v) = (1-w(v)).*Sa(i0(v)) + w(v).*Sa(i0(v)+1);
        img = img + iv .* exp(1j*kc*dr);
    end
end

function vol = bpVol(rc, r0v, x, y, z, VX, VY, VZ, drAxis, kc)
% 3D BP บน voxel grid (สำหรับ tomography)
    Nzp = numel(drAxis); ddr = drAxis(2)-drAxis(1); d1 = drAxis(1);
    vol = zeros(size(VX));
    for p = 1:numel(r0v)
        dr = sqrt((VX-x(p)).^2 + (VY-y(p)).^2 + (VZ-z(p)).^2) - r0v(p);
        fi = (dr - d1)/ddr + 1;
        i0 = floor(fi);  w = fi - i0;
        v  = (i0 >= 1) & (i0 < Nzp);
        Sa = rc(:,p);  iv = zeros(size(VX));
        iv(v) = (1-w(v)).*Sa(i0(v)) + w(v).*Sa(i0(v)+1);
        vol = vol + iv .* exp(1j*kc*dr);
    end
end

function [ent, pm] = imgMetrics(img)
% entropy (ต่ำ = โฟกัสดี) + peak/mean (สูง = คอนทราสต์ดี)
    m  = abs(img);  m = m/max(m(:));
    P  = m.^2 / sum(m(:).^2);
    ent = -sum(P(:).*log(P(:)+1e-12));
    pm  = 20*log10(max(m(:))/mean(m(:)));
end

function [wx, wy] = irfWidths(M, gx, gy)
% -3 dB width ผ่านจุด peak: แกน x (~range) และแกน y (~cross)
    [~, pk] = max(M(:));  [iy, ix] = ind2sub(size(M), pk);
    wx = w3(gx, M(iy,:));
    wy = w3(gy, M(:,ix));
end

function w = w3(ax, cut)
    ax = ax(:); m = abs(cut(:)); d = 20*log10(m/max(m)+eps);
    a = ax(d >= -3);
    if isempty(a), w = 0; else, w = max(a)-min(a); end
end
