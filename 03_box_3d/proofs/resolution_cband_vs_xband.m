%% CircularSAR_ResCompare — วัด RESOLUTION จริง: C-band vs X-band บน POINT target
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%
%  ทำไมต้องมีไฟล์นี้: W7_CircularSAR(_Xband) ใช้ฉากเป็น "กล่อง rough" (เฟสสุ่ม)
%  = distributed target มี SPECKLE เต็ม -> วัด -3dB ที่ยอด speckle ไม่ได้
%  (ยิ่ง X-band resolution ละเอียดกว่า voxel 1 m -> peakExtents คืน 0.00)
%
%  วิธีที่ถูก: ยิง "จุดเดียว" ที่ scene center แล้ววัด IRF (impulse response)
%  ด้วย grid ละเอียด (range/cross/height ~1-2 cm) ตามแกน -> mainlobe คมจริง
%  รันทั้ง C-band 4 GHz และ X-band 9.6 GHz ด้วย geometry เดียวกันเป๊ะ
%  -> โชว์ว่า band สั้นลงทำให้ cross + height คมขึ้น ~lambda ratio (7.5/3.1 = 2.4x)
%     ส่วน range เท่าเดิม (ขึ้นกับ bandwidth ไม่ใช่ band)
%
%  หมายเหตุ weighting: ใส่ Taylor -35 dB บนแกนความถี่ (range) เหมือนงานจริง
%     -> range -3dB จะกว้างกว่า c/2BW ประมาณ 1.3x | แกน cross/height ไม่ taper (uniform)
%  RUNTIME: ~30-60 วินาที (จุดเดียว, วัดเป็นเส้น 1D ไม่ใช่ทั้ง volume)

clear; clc; close all;
c = physconst('LightSpeed');

%% ===== GEOMETRY เดียวกันทั้งสอง band (fair) =====
Rg=1000; baseH=500; Npass=8; dzPass=2.5; arcDeg=15; azStep=0.05;
arcs=0:arcDeg:360-arcDeg; Narc=numel(arcs); NazArc=round(arcDeg/azStep);
B_elev=(Npass-1)*dzPass;                       % 17.5 m
sceneCtr=[1000;0;0];                           % [range; cross; height]
pt=[1000;0;0];  amp=1;                         % จุดเป้าเดี่ยวที่ center (z=0)
r0=sqrt(Rg^2+baseH^2); th=atan(baseH/Rg); cos2=cos(th)^2;
fprintf('Geometry: %d arcs x %d deg | %d passes (baseline %.1f m) | r0 %.0f m\n', ...
    Narc, arcDeg, Npass, B_elev, r0);

%% ===== PROCESSING (BW เท่ากันทั้งสอง band) =====
rngres=3; bw=c/(2*rngres); Nfreq=64; Nzp=512;  % BW 50 MHz

%% ===== เส้นวัด IRF ละเอียด (relative to จุดเป้า) =====
yLine=1000+(-6:0.02:6);  voxR=[yLine.',      zeros(numel(yLine),1), zeros(numel(yLine),1)];  % range (step 0.02 พอจับ X-band ~0.08 m)
xLine=(-3:0.02:3);       voxC=[1000*ones(numel(xLine),1), xLine.',  zeros(numel(xLine),1)];   % cross
zLine=(-3:0.02:3);       voxZ=[1000*ones(numel(zLine),1), zeros(numel(zLine),1), zLine.'];    % height

bands=[4e9, 9.6e9];  names={'C-band 4.0 GHz','X-band 9.6 GHz'};
res=zeros(2,3);  profs=cell(2,3);
figure('Color','w','Position',[80 80 1300 420]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for b=1:2
    fc=bands(b); lambda=c/fc; kc=4*pi*fc/c;
    fgrid=fc+linspace(-bw/2,bw/2,Nfreq); df=fgrid(2)-fgrid(1);
    wFreq=taylorwin(Nfreq,4,-35).';
    drAxis=((0:Nzp-1).'-Nzp/2)*c/(2*df*Nzp);

    accR=zeros(numel(yLine),1); accC=zeros(numel(xLine),1); accZ=zeros(numel(zLine),1);
    tB=tic;
    for iArc=1:Narc
        azv=arcs(iArc)+(0:NazArc-1)*azStep;
        platXY=sceneCtr(1:2)+Rg*[cosd(azv);sind(azv)];
        bpR=0; bpC=0; bpZ=0;
        for p=1:Npass
            Hp=baseH+(p-1)*dzPass;  plat=[platXY; Hp*ones(1,NazArc)];
            [rc,rref]=simPassFreq(pt,amp,plat,sceneCtr,fgrid,wFreq,Nzp,c);
            bpR=bpR+bpLine(rc,rref,plat,drAxis,voxR,kc);
            bpC=bpC+bpLine(rc,rref,plat,drAxis,voxC,kc);
            bpZ=bpZ+bpLine(rc,rref,plat,drAxis,voxZ,kc);
        end
        accR=accR+abs(bpR).^2;  accC=accC+abs(bpC).^2;  accZ=accZ+abs(bpZ).^2;  % multilook
    end
    mR=sqrt(accR/Narc); mC=sqrt(accC/Narc); mZ=sqrt(accZ/Narc);
    res(b,1)=w3(yLine-1000,mR);  res(b,2)=w3(xLine,mC);  res(b,3)=w3(zLine,mZ);
    profs{b,1}=mR; profs{b,2}=mC; profs{b,3}=mZ;
    fprintf('%s: lambda %.1f cm | -3dB range %.2f | cross %.3f | height %.2f m | %.0f s\n', ...
        names{b}, lambda*100, res(b,1), res(b,2), res(b,3), toc(tB));
end

%% ===== FIGURE: overlay IRF C vs X (range / cross / height) =====
axName={'range \Deltay (m)','cross \Deltax (m)','height \Deltaz (m)'};
axVec ={yLine-1000, xLine, zLine};
for k=1:3
    nexttile; hold on;
    plot(axVec{k}, 20*log10(profs{1,k}/max(profs{1,k})+eps),'b','LineWidth',1.5);
    plot(axVec{k}, 20*log10(profs{2,k}/max(profs{2,k})+eps),'r','LineWidth',1.5);
    yline(-3,'k:'); grid on; ylim([-30 1]);
    if k==1, xlim([-6 6]); else, xlim([-3 3]); end
    xlabel(axName{k}); ylabel('dB');
    title(sprintf('%s\\newlineC %.2f  vs  X %.2f m', axName{k}, res(1,k), res(2,k)));
    if k==1, legend(names,'Location','south'); end
end
sgtitle('Circular SAR — IRF บน point target: C-band (น้ำเงิน) vs X-band (แดง)');

%% ===== TABLE =====
fprintf('\n========== RESOLUTION COMPARISON (-3 dB, measured) ==========\n');
fprintf('  %-16s  %-9s  %-9s  %-9s\n','band','range','cross','height');
for b=1:2
    fprintf('  %-16s  %6.2f m   %6.3f m   %6.2f m\n', names{b}, res(b,1),res(b,2),res(b,3));
end
fprintf('  %-16s  %6.2fx   %6.2fx   %6.2fx  (C/X)\n','ratio', ...
    res(1,1)/max(res(2,1),eps), res(1,2)/max(res(2,2),eps), res(1,3)/max(res(2,3),eps));
fprintf('\n[Theory height res] C %.2f m | X %.2f m  (lambda*r0/(2*B*cos^2))\n', ...
    (c/4e9)*r0/(2*B_elev*cos2), (c/9.6e9)*r0/(2*B_elev*cos2));
fprintf('[อ่านผล] range ควรเท่ากัน (BW เท่ากัน) | cross+height X-band คมกว่า ~2.4x (lambda ratio 7.5/3.1)\n');

figDir=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'figure');
if ~exist(figDir,'dir'); mkdir(figDir); end
exportgraphics(gcf, fullfile(figDir,'fig_res_compare_CvsX.png'),'Resolution',150);
fprintf('\nFigure saved: %s\n', fullfile(figDir,'fig_res_compare_CvsX.png'));


%%%% ===== LOCAL FUNCTIONS =====

function [rc, rref] = simPassFreq(pos, amp, plat, sceneCtr, fgrid, wFreq, Nzp, c)
    Naz=size(plat,2);
    rref=sqrt(sum((plat-sceneCtr).^2,1));
    dR=sqrt((pos(1,:)'-plat(1,:)).^2+(pos(2,:)'-plat(2,:)).^2+(pos(3,:)'-plat(3,:)).^2)-rref;
    Nf=numel(fgrid); S=zeros(Nf,Naz);
    for n=1:Nf, S(n,:)=wFreq(n)*(amp*exp(-1j*4*pi*fgrid(n)/c*dR)); end
    rc=fftshift(ifft(S,Nzp,1),1);
end

function out = bpLine(rc, rref, plat, drAxis, voxT, kc)
% BP บนชุด voxel เป็น "เส้น" (voxT = Nx3 [range,cross,height]) คืน complex Nx1
    Nzp=numel(drAxis); ddr=drAxis(2)-drAxis(1); d1=drAxis(1);
    N=size(voxT,1); acc=zeros(N,1);
    for a=1:size(plat,2)
        d=voxT-plat(:,a).';
        dv=sqrt(sum(d.^2,2))-rref(a);
        fi=(dv-d1)/ddr+1; i0=floor(fi); w=fi-i0;
        v=(i0>=1)&(i0<Nzp); Sa=rc(:,a); iv=zeros(N,1);
        iv(v)=(1-w(v)).*Sa(i0(v))+w(v).*Sa(i0(v)+1);
        acc=acc+iv.*exp(1j*kc*dv);
    end
    out=acc;
end

function w = w3(ax, cut)
    ax=ax(:); m=abs(cut(:)); d=20*log10(m/max(m)+eps); a=ax(d>=-3);
    if isempty(a), w=0; else, w=max(a)-min(a); end
end
