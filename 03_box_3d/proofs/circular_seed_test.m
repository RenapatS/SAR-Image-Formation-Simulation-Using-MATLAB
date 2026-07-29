%% W7 SEED TEST — พิสูจน์ว่า "โครงกล่องจริง ไม่ขึ้นกับเลขสุ่ม (seed)"
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  คำถาม: เราใส่ rough (สุ่มเฟสรายจุด) เข้าไป -> จะรู้ได้ไงว่ากล่องที่เห็น
%          เป็นของจริง ไม่ใช่ artifact ของเลขสุ่ม?
%  วิธีพิสูจน์: รัน circular SAR ด้วยหลาย seed แล้วดูว่า
%     (1) ผนังตกตำแหน่งเดิมไหม   (2) ความแรงรายหน้าใกล้กันไหม
%     (3) ภาพ top-down 2 seed correlate กันสูงไหม
%  ถ้าโครงนิ่ง เปลี่ยนแค่ speckle -> กล่อง = ของจริง (ล็อกด้วยเรขาคณิต)
%
%  ใช้ pipeline เดียวกับ W7_CircularSAR.m (freq-domain + Taylor + multilook)
%  แต่ลดโหลด (12 arcs, 4 passes) ให้รันหลาย seed ไหว ~1-2 นาที

clear; clc; close all;

seeds = [1, 7, 42];                 % << เปลี่ยน/เพิ่ม seed ได้
Nseed = numel(seeds);

%% ===== RADAR / FREQ-DOMAIN MODEL (เหมือน circular) =====
c   = physconst('LightSpeed');
fc  = 4e9;   lambda = c/fc;
rngres = 3;  bw = c/(2*rngres);
Nfreq = 64;  fgrid = fc + linspace(-bw/2, bw/2, Nfreq);
df    = fgrid(2)-fgrid(1);
wFreq = taylorwin(Nfreq, 4, -35).';
Nzp   = 512;
drAxis = ((0:Nzp-1).' - Nzp/2) * c/(2*df*Nzp);
kc = 4*pi*fc/c;

%% ===== CIRCULAR GEOMETRY (ลดโหลดสำหรับเทสหลาย seed) =====
Rg = 1000;  baseH = 500;
Npass = 8;  dzPass = 2.5;           % baseline 17.5 m (เท่า circular จริง)
arcDeg = 15; azStep = 0.1;          % 24 arcs, 150 pulses/arc (คมขึ้น ~ลด fuzz)
arcs = 0 : arcDeg : 360-arcDeg;  Narc = numel(arcs);
NazArc = round(arcDeg/azStep);
fprintf('Seed test: %d seeds | %d arcs | %d passes | %d pulses/arc\n', ...
    Nseed, Narc, Npass, NazArc);

%% ===== BOX =====
boxCenter = [1000; 0; 0];  boxL = 30; boxW = 20; boxH = 10;
Nf = 20;
faceNormal = [ -1 0 0;  1 0 0;  0 -1 0;  0 1 0;  0 0 1;  0 0 -1];
faceLabels = {'Front','Back','Left','Right','Top','Bottom'};

%% ===== VOXEL GRID =====
voxStep = 1.0;
xScene = (boxCenter(2)-25):voxStep:(boxCenter(2)+25);
yScene = (boxCenter(1)-25):voxStep:(boxCenter(1)+25);
zScene = -5:voxStep:20;
[Xg,Yg,Zg] = meshgrid(xScene,yScene,zScene);
sz = size(Xg);  voxT = [Yg(:), Xg(:), Zg(:)];
[Xe,Ye,Ze] = boxEdges(boxCenter, boxL/2, boxW/2, boxH);

% จุดวัดกึ่งกลางแต่ละหน้า + จุดภายใน (ควรต่ำ)
faceCtr = [boxCenter(1)-boxL/2, boxCenter(2),        boxH/2;    % Front
           boxCenter(1)+boxL/2, boxCenter(2),        boxH/2;    % Back
           boxCenter(1),        boxCenter(2)-boxW/2, boxH/2;    % Left
           boxCenter(1),        boxCenter(2)+boxW/2, boxH/2;    % Right
           boxCenter(1),        boxCenter(2),        boxH;      % Top
           boxCenter(1),        boxCenter(2),        boxH/2];   % interior

%% ===== รันแต่ละ seed =====
magAll = cell(1,Nseed);  mipAll = cell(1,Nseed);  magSmAll = cell(1,Nseed);
faceLvl = zeros(6,Nseed);  contrast = zeros(1,Nseed);
tAll = tic;
for si = 1:Nseed
    [scatPosAll, scatAmpAll, scatFace] = buildBoxScatterers( ...
        boxCenter, boxL, boxW, boxH, Nf, seeds(si));

    Icomb = zeros(sz);
    for iArc = 1:Narc
        azv = arcs(iArc) + (0:NazArc-1)*azStep;
        platXY = boxCenter(1:2) + Rg*[cosd(azv); sind(azv)];
        Hm = baseH + (Npass-1)/2*dzPass;
        losA = [platXY; Hm*ones(1,NazArc)] - boxCenter;
        visF = false(1,6);
        for f = 1:6, visF(f) = any(faceNormal(f,:)*losA > 0); end
        ampA = scatAmpAll .* visF(scatFace);
        keep = ampA ~= 0;  posK = scatPosAll(:,keep);  ampK = ampA(keep);

        bpA = zeros(sz);
        for p = 1:Npass
            Hp = baseH + (p-1)*dzPass;
            plat = [platXY; Hp*ones(1,NazArc)];
            [rc, rref] = simPassFreq(posK, ampK, plat, boxCenter, fgrid, wFreq, Nzp, c);
            bpA = bpA + bpOnePassDiff(rc, rref, plat, drAxis, voxT, kc, sz);
        end
        Icomb = Icomb + abs(bpA).^2;          % multilook
    end
    magC = sqrt(Icomb/Narc);  magC = magC ./ max(magC(:));
    magAll{si} = magC;
    mipAll{si} = squeeze(max(magC,[],3));                 % top-down MIP (ใช้หา correlation)
    magSm = smooth3(magC,'gaussian',[5 5 5],1.5);         % smooth เพื่อโชว์ 3D (เหมือน circular Fig3)
    magSmAll{si} = magSm ./ max(magSm(:));

    for f = 1:6
        v = interp3(Xg,Yg,Zg, magC, faceCtr(f,2), faceCtr(f,1), faceCtr(f,3), 'linear', 0);
        faceLvl(f,si) = 20*log10(v + eps);
    end
    contrast(si) = (faceLvl(1,si)+faceLvl(3,si)+faceLvl(4,si))/3 - faceLvl(6,si); % ผนัง−ภายใน
    fprintf('  seed %2d done (%.0f s)\n', seeds(si), toc(tAll));
end

%% ===== CORRELATION ระหว่าง seed (ภาพ top-down) =====
corrM = eye(Nseed);
for i = 1:Nseed
    for j = i+1:Nseed
        r = corr2m(mipAll{i}, mipAll{j});     % ไม่ง้อ toolbox
        corrM(i,j) = r;  corrM(j,i) = r;
    end
end

%% ===== FIGURE: กล่อง 3D ขาวดำ ต่อ seed (point cloud เหมือน circular Fig3) =====
pcThr = -16;                                     % dB (ขาว=แรง บนพื้นดำ)
xlims = [-18 18]; ylims = [982 1018]; zlims = [-3 15];
figure(1); set(gcf,'Name','Seed test — 3D box per seed','Color','w', ...
    'Position',[40 80 560*Nseed 640]);
tiledlayout(1,Nseed,'TileSpacing','compact','Padding','compact');
for si = 1:Nseed
    ax = nexttile;
    mdB = 20*log10(magSmAll{si}+eps);  idx = find(mdB >= pcThr);
    scatter3(Xg(idx), Yg(idx), Zg(idx), 10, mdB(idx), 'filled', 'MarkerFaceAlpha',0.4);
    hold on; plot3(Xe, Ye, Ze, '-', 'Color',[0.1 0.3 0.9], 'LineWidth',1.6); hold off;
    colormap(jet); clim([pcThr 0]);
    xlabel('Cross (m)'); ylabel('Range (m)'); zlabel('Height (m)');
    title(sprintf('seed = %d', seeds(si)));
    xlim(xlims); ylim(ylims); zlim(zlims); daspect([1 1 1]);
    grid on; box on; view(-35,22);
end
cb = colorbar; ylabel(cb,'dB'); cb.Layout.Tile='east';
sgtitle('W7 seed test — กล่อง 3D นิ่งทุก seed (เปลี่ยนแค่ speckle)');
rotate3d on;

figDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end
exportgraphics(figure(1), fullfile(figDir,'fig_w7_seedtest.png'), 'Resolution',150);

%% ===== ผลลัพธ์ =====
fprintf('\n========== SEED TEST RESULT ==========\n');
fprintf('[Recovered level ต่อหน้า ต่อ seed]  (0 dB = peak)\n');
fprintf('  %-8s', 'Face');  for si=1:Nseed, fprintf('%8d ', seeds(si)); end;  fprintf('  spread\n');
for f = 1:6
    fprintf('  %-8s', faceLabels{f});
    fprintf('%9.1f', faceLvl(f,:));
    fprintf('   %.1f dB\n', max(faceLvl(f,:))-min(faceLvl(f,:)));
end
fprintf('\n[Wall-to-interior contrast]  '); fprintf('%.1f  ', contrast); fprintf('dB\n');
fprintf('[Correlation ภาพ top-down ระหว่าง seed]\n');
disp(round(corrM,3));
fprintf('-> ผนังตำแหน่งเดิม + contrast สูงทุก seed + corr ~1  =>  กล่องเป็นของจริง\n');
fprintf('   (เปลี่ยน seed = เปลี่ยนแค่ speckle, โครงกล่องล็อกด้วยเรขาคณิต)\n');


%%%% ===== LOCAL FUNCTIONS (ยกจาก W7_CircularSAR + เพิ่ม seed) =====

function [pos, amp, faceIdx] = buildBoxScatterers(bc, L, Wd, H, Nf, seed)
% จุดครบ 6 หน้า, rough (สุ่มเฟส) ด้วย seed ที่กำหนด
    rng(seed);
    hL=L/2; hW=Wd/2; hH=H/2; u=linspace(-1,1,Nf); pf=Nf^2;
    [V,W]=meshgrid(u*hW,u*hH);
    F=cell(1,6);
    F{1}=[(bc(1)-hL)*ones(1,pf); (bc(2)+V(:))'; (bc(3)+hH+W(:))'];
    F{2}=[(bc(1)+hL)*ones(1,pf); (bc(2)+V(:))'; (bc(3)+hH+W(:))'];
    [U,W]=meshgrid(u*hL,u*hH);
    F{3}=[(bc(1)+U(:))'; (bc(2)-hW)*ones(1,pf); (bc(3)+hH+W(:))'];
    F{4}=[(bc(1)+U(:))'; (bc(2)+hW)*ones(1,pf); (bc(3)+hH+W(:))'];
    [U,V]=meshgrid(u*hL,u*hW);
    F{5}=[(bc(1)+U(:))'; (bc(2)+V(:))'; (bc(3)+H)*ones(1,pf)];
    F{6}=[(bc(1)+U(:))'; (bc(2)+V(:))'; zeros(1,pf)];
    pos=[F{:}];  faceIdx = repelem(1:6, pf);
    amp = exp(1j*2*pi*rand(1,6*pf));            % rough: สุ่มเฟส (ขนาด=1)
end


function [rc, rref] = simPassFreq(pos, amp, plat, sceneCtr, fgrid, wFreq, Nzp, c)
    Naz  = size(plat,2);
    rref = sqrt(sum((plat - sceneCtr).^2, 1));
    dR   = sqrt( (pos(1,:)'-plat(1,:)).^2 + (pos(2,:)'-plat(2,:)).^2 ...
               + (pos(3,:)'-plat(3,:)).^2 ) - rref;
    Nf = numel(fgrid);  S = zeros(Nf, Naz);
    for n = 1:Nf
        S(n,:) = wFreq(n) * (amp * exp(-1j*4*pi*fgrid(n)/c * dR));
    end
    rc = fftshift(ifft(S, Nzp, 1), 1);
end


function img = bpOnePassDiff(rc, rref, plat, drAxis, voxT, kc, sz)
    Nzp = numel(drAxis); ddr = drAxis(2)-drAxis(1); d1 = drAxis(1);
    Nvox = size(voxT,1);  acc = zeros(Nvox,1);
    for a = 1:size(plat,2)
        d  = voxT - plat(:,a).';
        dv = sqrt(sum(d.^2, 2)) - rref(a);
        fi = (dv - d1)/ddr + 1;  i0 = floor(fi);
        v  = (i0 >= 1) & (i0 < Nzp);  w = fi - i0;
        Sa = rc(:,a);  iv = zeros(Nvox,1);
        iv(v) = (1-w(v)).*Sa(i0(v)) + w(v).*Sa(i0(v)+1);
        acc = acc + iv .* exp(1j*kc*dv);
    end
    img = reshape(acc, sz);
end


function r = corr2m(A,B)
% normalized cross-correlation (แทน corr2 ไม่ต้องมี Image Processing Toolbox)
    a=A(:)-mean(A(:)); b=B(:)-mean(B(:));
    r=(a'*b)/sqrt((a'*a)*(b'*b)+eps);
end


function [Xe,Ye,Ze] = boxEdges(bc, hl, hw, H)
    cx=bc(2); cy=bc(1);
    C=[cx-hw,cy-hl,0; cx+hw,cy-hl,0; cx+hw,cy+hl,0; cx-hw,cy+hl,0; ...
       cx-hw,cy-hl,H; cx+hw,cy-hl,H; cx+hw,cy+hl,H; cx-hw,cy+hl,H];
    E=[1 2;2 3;3 4;4 1; 5 6;6 7;7 8;8 5; 1 5;2 6;3 7;4 8];
    Xe=nan(1,3*size(E,1)); Ye=Xe; Ze=Xe;
    for e=1:size(E,1)
        j=3*(e-1)+1;
        Xe(j:j+1)=C(E(e,:),1); Ye(j:j+1)=C(E(e,:),2); Ze(j:j+1)=C(E(e,:),3);
    end
end
