%% GOTCHA real-SAR Backprojection — ตัวอย่างสั้นที่สุด (~50 บรรทัด)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%
%  พิสูจน์แล้ว (Week 6): BP เดียวกับ try104/try108 ใช้กับข้อมูลจริงได้
%  ข้อมูล: AFRL GOTCHA (parking lot, X-band) — subset จาก RITSAR
%          อ่านตรงจาก data/GOTCHA-CP_Disc1/ — โหลด az001-004 ของ pass 1 (HH)
%
%  โครงสร้างไฟล์ .mat: struct 'data' มี field:
%     fp   (424 x 117) complex  = phase history (freq x pulse)
%     freq (424 x 1)            = ความถี่ (Hz)  -> X-band ~9.6 GHz, BW 624 MHz
%     x,y,z (1 x 117)           = ตำแหน่ง antenna 3D ต่อ pulse (m)
%     r0   (1 x 117)            = ระยะ antenna -> scene center
%     af                        = autofocus corrections (r_correct, ph_correct)
%
%  ขั้นถัดไป: gotcha_bp_2d.m ใส่ windowing + autofocus + เทียบ sim vs real
%             แล้วต่อด้วย gotcha_sparse3d_showcase.m สำหรับภาพ 3 มิติ

clear; clc; close all;
c = physconst('LightSpeed');

dataDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'GOTCHA-CP_Disc1', 'DATA', 'pass1', 'HH');
files = dir(fullfile(dataDir, 'data_3dsar_pass1_az00[1-4]_HH.mat'));

fp=[]; x=[]; y=[]; z=[]; r0=[]; freq=[];
for k = 1:numel(files)
    S = load(fullfile(dataDir, files(k).name));  d = S.data;
    fp = [fp, d.fp];  freq = d.freq;
    x = [x, d.x];  y = [y, d.y];  z = [z, d.z];  r0 = [r0, d.r0];
end
[Nf, Np] = size(fp);  df = freq(2)-freq(1);  fcen = mean(freq);
fprintf('GOTCHA: %d freq x %d pulses | fc=%.2f GHz | BW=%.0f MHz | range res ~%.2f m\n', ...
    Nf, Np, fcen/1e9, Nf*df/1e6, c/(2*Nf*df));

% --- range compression (zero-padded IFFT over frequency) ---
Zp = 8;
rc = fftshift(ifft(fp, Nf*Zp, 1), 1);
dr_axis = ((0:Nf*Zp-1).' - Nf*Zp/2) * (c/(2*df)) / (Nf*Zp);   % differential range (m)

% --- backprojection onto ground grid (z = 0) ---
N = 300;  W = 50;  g = linspace(-W, W, N);  [GX, GY] = meshgrid(g, g);
img = zeros(N, N);  kc = 4*pi*fcen/c;
for p = 1:Np
    R  = sqrt((GX-x(p)).^2 + (GY-y(p)).^2 + z(p)^2);
    dr = R(:) - r0(p);
    iv = interp1(dr_axis, rc(:,p), dr, 'linear', 0);
    img = img + reshape(iv .* exp(1j*kc*dr), N, N);
end

M = 20*log10(abs(img)/max(abs(img(:))) + 1e-6);
figure; imagesc(g, g, M); axis xy image; colormap(gray); clim([-30 0]); colorbar;
title('REAL GOTCHA SAR — parking lot (BP from .mat phase history)');
xlabel('cross-range (m)'); ylabel('range (m)');
fprintf('Done. Real SAR image formed from GOTCHA phase history.\n');
