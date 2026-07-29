%% W10 #2b — PRIZM CHECK: รถคัน J กว้างเท่าไหร่กันแน่
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%
%  ###################################################################
%  ##  ผลรัน: วิธี edge-density ในไฟล์นี้ "ไม่ผ่านตัวคุม" — เป็นโมฆะ  ##
%  ##  B วัดได้ 1.219 (ควร ~0.87) | E ได้ 0.750 (ควร 0.855)          ##
%  ##  ห้ามใช้คำตัดสินจากสคริปต์นี้                                  ##
%  ##                                                                 ##
%  ##  หลักฐานที่ดีกว่า (จาก estimator W p10-90 ที่ validate แล้ว):     ##
%  ##  เก๋งทุกคันอ่าน "เกิน" ความกว้างจริง +0.03..+0.07 เป็นระบบ       ##
%  ##  J อ่าน 1.59: ถ้าจริง 1.54 -> +0.05 เข้าแพทเทิร์น               ##
%  ##             ถ้าจริง 1.69 -> -0.10 ผิดแพทเทิร์นคันเดียว          ##
%  ##  => น้ำหนักเทไปทาง xls (1.54) แต่ไม่เด็ดขาด | ดูรูปถ่ายประกอบ    ##
%  ###################################################################
%
%  ปม (เจอตอนตรวจ spec กับเว็บ, W10):
%    - AFRL xls ชีต "Vehicle Dimensions": J = Chevy Prizm กว้าง 1.54 m
%    - AFRL สำรวจพิกัดมุมรถ 4 จุด (ชีต Vehicles): คำนวณกลับได้ 1.54 m เช่นกัน
%    - แต่สเปกโรงงาน Chevy Prizm 1998-2002 = 66.7 นิ้ว = 1.69 m
%    - รถอีก 5 คัน xls ต่ำกว่าสเปกโรงงานแค่ 0.00-0.07 m | J ต่ำกว่าถึง 0.154 m = โดดคันเดียว
%    - เรดาร์เราวัดได้ 1.59 (bias กลุ่ม +0.05 -> ประมาณของจริง 1.54) = เข้าข้าง xls
%
%  คำถามที่เรดาร์ตอบได้: "ขอบข้างรถอยู่ตรงไหน"
%  เส้นสว่างที่ขอบข้างรถคือ double-bounce ระหว่างตัวถังกับพื้น ซึ่งอยู่ตำแหน่งจริง
%  (ไม่โดน layover) -> หา peak ของความหนาแน่น voxel ตามแกน |u| แล้วเทียบ
%  สมมติฐาน A: ขอบที่ |u| = 1.54/2 = 0.77   สมมติฐาน B: |u| = 1.69/2 = 0.845
%  ต่างกัน 7.5 cm ต่อข้าง = 0.6 voxel -> ขอบเขตความสามารถพอดี ๆ ต้องอ่านอย่างระวัง
%  ใช้ B (Camry) กับ E (Sentra) เป็นตัวคุม: รู้แล้วว่า xls ตรงสเปกโรงงาน
%
%  ใช้ cache เดิม ไม่รัน recon

clear; clc; close all;

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % root ของ repo
S = load(fullfile(root, 'figure', 'w9cmp2_cache.mat'));
R = S.R;  veh = S.veh;
figDir = fullfile(root, 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end

% {id, W ตาม xls, W ตามสเปกโรงงาน(เว็บ, ตรวจ W10)}
tests = {'B', 1.74, 1.80;      % ตัวคุม 1 (Camry — xls ใกล้สเปก)
         'E', 1.71, 1.71;      % ตัวคุม 2 (Sentra — ตรงเป๊ะ)
         'J', 1.54, 1.69};     % ตัวปัญหา
nT = size(tests,1);
metDB = 20;
uax = 0 : 0.03125 : 1.4;       % แกน |u| ละเอียด 1/4 voxel (จุดจริงกระจายใน bin)

fprintf('===== W10 #2b | ตำแหน่งขอบข้างรถจากความหนาแน่น voxel =====\n\n');
fprintf('%-4s %10s %10s | %12s %14s %12s\n', 'ID', 'xls W/2', 'spec W/2', ...
    'ขอบที่วัดได้', 'เข้าข้าง xls?', 'nVox');
fprintf('%s\n', repmat('-', 1, 72));

f1 = figure(1); clf; set(f1,'Color','w','Position',[60 60 1100 360],'Name','Prizm width check');
tiledlayout(1, nT, 'TileSpacing','compact','Padding','compact');
edgeHat = nan(nT,1);
for t = 1:nT
    v  = find(strcmp(veh(:,2), tests{t,1}), 1);
    GL = veh{v,6};
    pu=R(v).pu; pv=R(v).pv; pz=R(v).pz; pm=R(v).pm;
    % เฉพาะจุดเข้ม ในช่วงความยาวรถ ระดับต่ำ (z < 0.8 = แถวเส้นขอบล่าง ตัด ghost ลอย)
    m = pm >= -metDB & abs(pv) <= GL/2 & pz >= -0.4 & pz <= 0.8;
    au = abs(pu(m));
    % ความหนาแน่นถ่วงน้ำหนักด้วยพลังงาน (จุดสว่าง = เส้น double-bounce เด่นกว่า)
    wgt = 10.^(pm(m)/20);
    dens = zeros(size(uax));
    for i = 1:numel(uax)
        sel = abs(au - uax(i)) <= 0.0625;          % หน้าต่าง +-ครึ่ง voxel
        dens(i) = sum(wgt(sel));
    end
    dens = dens / max(dens);
    % ขอบ = จุดที่ density ร่วงต่ำกว่าครึ่งหนึ่งของ peak ด้านนอกสุด
    [~, ipk] = max(dens);
    iEdge = find(dens(ipk:end) < 0.5, 1) + ipk - 1;
    if isempty(iEdge), iEdge = numel(uax); end
    edgeHat(t) = uax(iEdge);

    hxls = tests{t,2}/2;  hspec = tests{t,3}/2;
    dx = abs(edgeHat(t)-hxls);  ds = abs(edgeHat(t)-hspec);
    if abs(hxls - hspec) < 0.02
        verdict = '(แยกไม่ได้-ค่าใกล้กัน)';
    elseif dx < ds, verdict = 'xls';
    else,           verdict = 'spec โรงงาน';
    end
    fprintf('%-4s %10.3f %10.3f | %12.3f %14s %12d\n', tests{t,1}, hxls, hspec, ...
        edgeHat(t), verdict, nnz(m));

    nexttile; hold on;
    plot(uax, dens, 'LineWidth', 1.3);
    xline(hxls, '--', sprintf('xls %.2f', hxls), 'Color', [0.15 0.45 0.7], ...
        'LabelOrientation','horizontal');
    if abs(hxls-hspec) >= 0.02
        xline(hspec, '--', sprintf('spec %.2f', hspec), 'Color', [0.75 0.35 0.15], ...
            'LabelOrientation','horizontal');
    end
    xline(edgeHat(t), '-', 'วัดได้', 'Color', [0.2 0.2 0.2]);
    hold off; grid on; xlim([0 1.3]); ylim([0 1.05]);
    xlabel('|u| จากแกนกลางรถ (m)'); ylabel('ความหนาแน่นพลังงาน (norm)');
    title(sprintf('%s (%s)', veh{v,1}, veh{v,2}));
end
sgtitle('W10 #2b | ขอบข้างรถ: เส้นทึบ=วัดได้ | ระวัง: ช่องต่าง 7.5 cm = 0.6 voxel เท่านั้น');
exportgraphics(f1, fullfile(figDir,'fig_w10_prizm_width.png'), 'Resolution', 150);

fprintf(['\nวิธีอ่าน:\n' ...
  '1. ตัวคุม B/E ต้องให้ขอบใกล้ค่า xls ของมัน (คลาด <= ~ครึ่ง voxel) ไม่งั้นวิธีวัดนี้เชื่อไม่ได้\n' ...
  '2. ถ้า J เข้าข้าง 0.77 (xls) = ของจริงในลานแคบกว่า Prizm มาตรฐาน 15 cm\n' ...
  '   -> สามแหล่งอิสระ (ชีต xls, สำรวจมุมรถ, เรดาร์) ตรงกัน: ป้ายชื่อรุ่นอาจไม่ตรงรถจริง\n' ...
  '   (แบบเดียวกับ finding G/H ที่เอกสารกับข้อมูลขัดกัน) — หรือรถถูกถอดกระจก/ชิ้นข้าง\n' ...
  '3. ถ้า J เข้าข้าง 0.845 (spec) = xls พิมพ์ผิดทั้งสองชีต และ bias เรดาร์ของ J ต่างจากกลุ่ม\n' ...
  '4. J อยู่ไกลสุด (55 m) SNR ต่ำสุด — ถ้า density เตี้ยแบน ไม่มี peak ชัด ให้ถือว่าสรุปไม่ได้\n']);
fprintf('\nรูปถ่าย AFRL ของ J: figure/car_photos/ppt_img_08.jpg — เปิดดูประกอบ (เทียบสัดส่วนกับ Prizm จริง)\n');
fprintf('เซฟรูป: fig_w10_prizm_width.png\n');
