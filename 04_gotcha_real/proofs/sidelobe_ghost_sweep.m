%% W10 #1b — SIDELOBE CHECK: กวาดเกณฑ์ dB ดูว่า ghost แนวสูงยุบตรงไหน
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%
%  ที่มา: W10_PhaseRefine พิสูจน์แล้วว่า ghost แนวสูงไม่ใช่ phase error
%  แต่เป็น sidelobe ของ elevation aperture ที่แคบ (8 pass / ~1.7-2.7 องศา):
%      sidelobe แรก  z +-0.91 m  ระดับ -8.9 dB
%      sidelobe สอง  z +-1.68 m  ระดับ -9.2 dB
%  เกณฑ์วัดเดิม -20 dB อยู่ "ใต้" sidelobe -> สำเนาปลอมถูกนับเป็นรถทั้งหมด
%
%  การทดสอบ: ถ้าคำอธิบายนี้ถูก ghost%% ต้องยุบฮวบเมื่อยกเกณฑ์ขึ้นเหนือ ~-9 dB
%  และต้องยุบ "พร้อมกันทุกคัน" (เพราะเป็นสมบัติของ aperture ไม่ใช่ของรถ)
%  ถ้ายกเกณฑ์แล้ว ghost ไม่ยุบ = มีอย่างอื่นอีก (เช่น multipath จริง)
%
%  ใช้ cache เดิม ไม่รัน recon — จบในไม่กี่วินาที
%  ข้อจำกัดที่รู้ล่วงหน้า: เกณฑ์สูง = voxel เหลือน้อย โดยเฉพาะรถไกล (SNR ต่ำ)
%  ตารางจึงพิมพ์ nVox กำกับทุกช่อง — ช่องที่จุดน้อยอ่านค่าไม่ได้ อย่าตีความ

clear; clc; close all;

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % root ของ repo
S = load(fullfile(root, 'figure', 'w9cmp2_cache.mat'));   % R, veh, nLook
R = S.R;  veh = S.veh;  Nveh = size(veh,1);
figDir = fullfile(root, 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end

GT_L = cell2mat(veh(:,6));  GT_W = cell2mat(veh(:,7));  GT_H = cell2mat(veh(:,8));
isVeh = ~ismember(veh(:,2), {'G','H'});

dbList = [8 10 12 14 16 20 25 30];       % เกณฑ์ที่กวาด (cache มีถึง -40)
vb  = -3:0.25:3;  vbc = vb(1:end-1)+0.125;

ghost = nan(Nveh, numel(dbList));
zp90  = nan(Nveh, numel(dbList));
roofC = nan(Nveh, numel(dbList));        % ความสูงห้องโดยสาร (กลางรถ p75 ของโปรไฟล์)
nvox  = zeros(Nveh, numel(dbList));

for v = 1:Nveh
    if ~isVeh(v), continue; end
    pu=R(v).pu; pv=R(v).pv; pz=R(v).pz; pm=R(v).pm;
    if isempty(pu), continue; end
    for e = 1:numel(dbList)
        thr = -dbList(e);
        m2 = pm >= thr & abs(pu) <= GT_W(v)/2 + 0.2;
        nvox(v,e) = nnz(m2);
        if nnz(m2) < 8, continue; end
        roof = nan(1, numel(vb)-1);
        for j = 1:numel(vb)-1
            zz = pz(m2 & pv >= vb(j) & pv < vb(j+1));
            if ~isempty(zz), roof(j) = pctl(zz, 90); end
        end
        onCar = abs(vbc) <= GT_L(v)/2 & isfinite(roof);
        if nnz(onCar) >= 4
            ghost(v,e) = 100 * mean(roof(onCar) > GT_H(v) + 0.15);
            cab = onCar & abs(vbc) <= GT_L(v)/4;
            if nnz(cab) >= 3, roofC(v,e) = pctl(roof(cab), 75); end
        end
        inb = abs(pu) <= GT_W(v)/2+0.15 & abs(pv) <= GT_L(v)/2+0.15 & ...
              pz >= -0.5 & pz <= 3 & pm >= thr;
        if nnz(inb) >= 8, zp90(v,e) = pctl(pz(inb), 90); end
    end
end

%% ===== ตาราง =====
fprintf('===== ghost%% (ช่วงโปรไฟล์ที่สูงเกินหลังคา) vs เกณฑ์ dB =====\n');
fprintf('sidelobe แนวสูงอยู่ที่ -8.9 dB -> ทำนาย: ghost ยุบเมื่อเกณฑ์ >= -8..-10 dB\n\n');
fprintf('%-4s', 'ID');
fprintf('  @-%-2d dB', dbList); fprintf('\n');
fprintf('%s\n', repmat('-', 1, 4 + 9*numel(dbList)));
for v = 1:Nveh
    if ~isVeh(v), continue; end
    fprintf('%-4s', veh{v,2});
    for e = 1:numel(dbList)
        if nvox(v,e) < 8
            fprintf('  %6s ', 'n/a');
        else
            fprintf('  %5.0f%% ', ghost(v,e));
        end
    end
    fprintf('\n');
end
fprintf('\nnVox ประกอบ (อ่านค่าไม่ได้ถ้าน้อย):\n');
fprintf('%-4s', 'ID'); fprintf('  @-%-2d dB', dbList); fprintf('\n');
for v = 1:Nveh
    if ~isVeh(v), continue; end
    fprintf('%-4s', veh{v,2}); fprintf('  %6d ', nvox(v,:)); fprintf('\n');
end

fprintf('\n===== ความสูงห้องโดยสาร (โปรไฟล์กลางรถ p75) vs เกณฑ์ =====\n');
fprintf('%-4s %6s |', 'ID', 'spec~'); fprintf('  @-%-2d dB', dbList); fprintf('\n');
fprintf('%s\n', repmat('-', 1, 13 + 9*numel(dbList)));
for v = 1:Nveh
    if ~isVeh(v), continue; end
    fprintf('%-4s %6.2f |', veh{v,2}, GT_H(v));
    for e = 1:numel(dbList)
        if isnan(roofC(v,e)), fprintf('  %6s ', '-');
        else, fprintf('  %6.2f ', roofC(v,e)); end
    end
    fprintf('\n');
end
fprintf(['\nวิธีอ่าน:\n' ...
  '1. ghost%% ยุบพร้อมกันทุกคันแถวเกณฑ์ -8..-12 dB = ยืนยัน sidelobe (คุณสมบัติ aperture)\n' ...
  '2. ถ้าห้องโดยสารที่เกณฑ์สูงเข้าใกล้ spec = จุดสว่างสุดคือตัวรถจริง sidelobe แค่จางกว่า\n' ...
  '3. รถไกล (D 43m, J 55m) n/a เร็วกว่ารถใกล้ = SNR ต่ำ -> ยืนยันกลไก bias ตามระยะ\n' ...
  '** เกณฑ์สูง = จุดน้อย ตัวเลขผันผวน อย่าใช้เป็นตัววัดหลัก — นี่คือการทดสอบกลไก ไม่ใช่ estimator ใหม่ **\n']);

%% ===== รูป =====
f1 = figure(1); clf; set(f1,'Color','w','Position',[60 60 1000 420],'Name','sidelobe check');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
cols = lines(Nveh);
for v = 1:Nveh
    if ~isVeh(v), continue; end
    ok = nvox(v,:) >= 8;
    plot(dbList(ok), ghost(v,ok), '-o', 'Color', cols(v,:), ...
        'DisplayName', sprintf('%s', veh{v,2}));
end
xline(8.9, '--', 'sidelobe -8.9', 'LabelOrientation','horizontal');
hold off; grid on; xlabel('เกณฑ์ (-dB)'); ylabel('ghost %');
legend('Location','eastoutside'); title('ghost ยุบตรงไหน');
nexttile; hold on;
for v = 1:Nveh
    if ~isVeh(v), continue; end
    ok = isfinite(roofC(v,:));
    plot(dbList(ok), roofC(v,ok), '-o', 'Color', cols(v,:), 'HandleVisibility','off');
    yline(GT_H(v), ':', 'Color', cols(v,:), 'HandleVisibility','off');
end
hold off; grid on; xlabel('เกณฑ์ (-dB)'); ylabel('ห้องโดยสาร (m)');
title('เส้นประ = spec ของแต่ละคัน');
sgtitle('W10 #1b | sidelobe check — จาก cache, ไม่รัน recon');
exportgraphics(f1, fullfile(figDir,'fig_w10_sidelobe_sweep.png'), 'Resolution', 150);
fprintf('\nเซฟรูป: fig_w10_sidelobe_sweep.png\n');

%% ===== LOCAL =====
function y = pctl(x, p)
x = sort(x(:)); n = numel(x);
if n == 0, y = nan(size(p)); return; end
if n == 1, y = repmat(x, size(p)); return; end
q = 100*((1:n).' - 0.5)/n;  pq = p(:);
y = interp1(q, x, pq, 'linear');
y(pq <= q(1)) = x(1);  y(pq >= q(end)) = x(end);
y = reshape(y, size(p));
end
