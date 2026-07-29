%% W9 helper — วัด "ตำแหน่งรถจากภาพ 3D" เทียบพิกัด xls (ไม่ต้องรัน recon ใหม่)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  โหลด figure/w9cmp2_cache.mat (เซฟไว้ตอนจบ W9_CompareVehicles) แล้วคำนวณ
%  centroid (median — ทนต่อ junk) ของ voxel เข้ม >= -20 dB ในเขต +-3 m
%  หมุนจาก car frame กลับเป็น global -> ได้พิกัดวัด (x,y) + ระยะคลาดจาก xls
%  หมายเหตุ: recon ยึดพิกัด xls เป็นศูนย์ ดังนั้นค่านี้คือ "offset ที่ข้อมูลจริง
%  เห็น" แบบเดียวกับที่ tophat โชว์ว่า calP อยู่บนขอบวงแหวน (เลื่อน ~1 m)

clear; clc;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % root ของ repo
S = load(fullfile(root, 'figure', 'w9cmp2_cache.mat'));   % V, R, axv, veh, nLook
veh = S.veh;  R = S.R;  Nveh = size(veh,1);

fprintf('%-16s %-3s | xls (x, y)          | วัดได้ (x, y)        | offset (m)\n', 'Vehicle', 'ID');
fprintf(repmat('-', 1, 78)); fprintf('\n');
for v = 1:Nveh
    cx = veh{v,3};  cy = veh{v,4};  hd = deg2rad(veh{v,5});
    m  = R(v).pm >= -20 & abs(R(v).pu) <= 3 & abs(R(v).pv) <= 3 & R(v).pz >= -0.5;
    if nnz(m) < 10
        fprintf('%-16s %-3s | (%8.3f, %8.3f) | %-20s | รถไม่อยู่ในฉาก/จุดน้อย (%d)\n', ...
            veh{v,1}, veh{v,2}, cx, cy, 'n/a', nnz(m));
        continue;
    end
    % median ใน car frame -> หมุนกลับ global (Rr = [cos -sin; sin cos])
    uc = median(R(v).pu(m));  vc = median(R(v).pv(m));
    Rr = [cos(hd) -sin(hd); sin(hd) cos(hd)];
    d  = Rr * [uc; vc];                       % offset ใน global
    mx = cx + d(1);  my = cy + d(2);
    fprintf('%-16s %-3s | (%8.3f, %8.3f) | (%8.3f, %8.3f) | %.2f\n', ...
        veh{v,1}, veh{v,2}, cx, cy, mx, my, hypot(d(1), d(2)));
end
fprintf(['\nตีความ: offset = ระยะระหว่างศูนย์กลางจุดสะท้อนจริงกับพิกัด xls\n' ...
    '(<0.3 m = ตรงระดับ voxel | ถ้าราว 1 m+ แปลว่าจุดอ้างอิง xls ไม่ใช่กลางรถ)\n']);
