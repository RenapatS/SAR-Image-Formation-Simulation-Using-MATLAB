%% W9 helper — แกะรูปถ่ายจริงที่ฝังใน Challenge_Pictures_Images.ppt (AFRL)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  .ppt ยุค 97-2003 เป็น OLE binary — รูป JPEG ฝังอยู่ในไฟล์ตรง ๆ
%  แค่ scan หา JPEG marker (FFD8FF ... FFD9) แล้ว dump ออกมา
%  ผลลัพธ์: figure/car_photos/ppt_img_NN.jpg + montage ให้ดูว่ารูปไหนคือรถอะไร
%  จากนั้นไปตั้ง photoFiles ใน W9_CompareVehicles.m ตาม index ที่เห็น

clear; clc; close all;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % root ของ repo
ppt  = fullfile(root, 'data', 'GOTCHA-CP_Disc1', 'DOCUMENTATION', 'Challenge_Pictures_Images.ppt');
assert(exist(ppt,'file') > 0, 'ไม่พบ %s', ppt);

outDir = fullfile(root, 'figure', 'car_photos');
if ~exist(outDir,'dir'); mkdir(outDir); end

fid = fopen(ppt, 'rb');  data = fread(fid, inf, '*uint8').';  fclose(fid);
soi = strfind(data, uint8([255 216 255]));      % JPEG start
eoi = strfind(data, uint8([255 217]));          % JPEG end
n = 0;  files = {};
for s = soi
    e = eoi(find(eoi > s+3, 1));                % EOI ตัวแรกหลัง SOI
    if isempty(e), continue; end
    blob = data(s : e+1);
    if numel(blob) < 15000, continue; end       % ข้าม thumbnail เล็ก ๆ
    if n > 0 && s < lastEnd, continue; end      % กัน JPEG ซ้อนใน JPEG เดิม
    n = n + 1;  lastEnd = e+1;
    files{n} = fullfile(outDir, sprintf('ppt_img_%02d.jpg', n)); %#ok<SAGROW>
    f2 = fopen(files{n}, 'wb');  fwrite(f2, blob);  fclose(f2);
end
fprintf('แกะได้ %d รูป -> %s\n', n, outDir);

% montage พร้อมเลข index
nc = ceil(sqrt(n));  nr = ceil(n/nc);
figure('Name','photos in ppt','Color','w','Position',[50 50 1400 800]);
tiledlayout(nr, nc, 'TileSpacing','tight','Padding','tight');
for k = 1:n
    nexttile;
    try, imshow(imread(files{k})); catch, axis off; end
    title(sprintf('%d', k), 'FontSize', 12, 'FontWeight','bold');
end
sgtitle('เลือกรูปของ Malibu (A) / Santa Fe (F) / Forklift (C2) แล้วจดเลข index');
fprintf('ดูรูปแล้วตั้งใน W9_CompareVehicles.m เช่น photoIdx = [5 8 11];\n');
