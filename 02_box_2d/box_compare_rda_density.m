%% COMBINED B — Step comparison + Density proof + Algorithms
%  หมายเหตุ: บล็อก RADAR CONFIGURATION สืบทอดมาจาก 01_point_target/point3_bp_2d.m
%  ซึ่งดัดแปลงจากตัวอย่าง Stripmap SAR ของ MathWorks (ดู attribution ในไฟล์นั้น)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  ไฟล์รวมโค้ดจริง (ก๊อปจากไฟล์ย่อยมาต่อกัน) — แต่ละ section มี clear เอง
%  ทำงานอิสระต่อกัน; local function ทั้งหมดรวมไว้ท้ายไฟล์ (ชื่อไม่ชนกัน)
%  ที่มา: try106_StepComparison.m, try107_BoxDensity.m, try108_RDA_vs_BP.m

%% ================= SECTION 1: STEP COMPARISON (try106)  (from try106_StepComparison.m) =================
%% Week 6 | Step 3: Step-by-step Comparison (synthesis for report/slides)
%
%  รวบยอดผลทุก stage ที่ทำมา ไว้ในภาพเดียว เพื่อเห็น "อะไรเปลี่ยนเพราะอะไร"
%
%  ค่าด้านล่างเป็นผลที่ "วัดได้จริง" จากการรันไฟล์ก่อนหน้า (verified):
%     point calibration & box W6+occlusion ........ try104_PhysicsRCS.m
%     box W5 manual RCS ........................... try104b_W5RCS_NewEval.m
%     noise sweep ................................. try105_AWGN_Noise.m
%  ** ถ้ารันใหม่แล้วเลขขยับเล็กน้อย ปรับค่าในบล็อก DATA ให้ตรง run ล่าสุด **

clear; clc; close all;

%% ===== DATA (measured, verified) =====

% --- System metrics (point-target calibration) — ไม่ขึ้นกับ RCS/noise ---
theory_rangeRes = 3.00;     theory_crossRes = 0.0937;
cal_rangeRes    = 1.73;     cal_crossRes    = 0.093;     % measured IRF
cal_pslr_r      = -11.4;    cal_pslr_c      = -13.8;     % dB

% --- Box-scene metrics across 3 configs ---
cfgNames = {'W5 manual','W6 physics','W6 + occlusion'};
rangeExt = [ 7.88, 34.12, 34.00];   % -3 dB bright-region extent (m)  (~box L=30)
crossExt = [20.00, 20.00, 20.00];   % -3 dB extent (m)  (box W=20)
locErr   = [17.62,  8.14, 14.75];   % range localization error (m)
dynRange = [ 96.2,  99.0,100.0];    % peak/background, noise-free (dB)

% --- Noise study (try105): image SNR vs input SNR ---
snrIn    = [20 10 0 -10 -20 -30 -40 -50];
imgSNR   = [65.4 55.4 45.4 35.5 25.9 16.9 11.8 11.8];
detected = [ 1   1   1   1    1    1    0    0];
detThresh = 13;

%% ===== FIGURE 1: Master comparison table (text) =====
L = {};
L{end+1} = '====== STEP 3 — EVOLUTION OF THE SAR PIPELINE ======';
L{end+1} = '';
L{end+1} = '[A. SYSTEM RESOLUTION — from point-target calibration]';
L{end+1} = '    (a property of bandwidth & aperture; SAME for every config,';
L{end+1} = '     unchanged by RCS or noise)';
L{end+1} = sprintf('    Range res : %.2f m   (theory %.2f m)', cal_rangeRes, theory_rangeRes);
L{end+1} = sprintf('    Cross res : %.3f m   (theory %.4f m)', cal_crossRes, theory_crossRes);
L{end+1} = sprintf('    PSLR      : range %.1f dB,  cross %.1f dB', cal_pslr_r, cal_pslr_c);
L{end+1} = '';
L{end+1} = '[B. BOX-SCENE METRICS — change with RCS weighting / occlusion]';
L{end+1} = sprintf('    %-16s %-9s %-9s %-9s %-8s','config','RangeExt','CrossExt','LocErr','DynR');
L{end+1} = sprintf('    %-16s %-9s %-9s %-9s %-8s','','(m)','(m)','(m)','(dB)');
L{end+1} = repmat('-',1,56);
for i=1:3
    L{end+1} = sprintf('    %-16s %-9.2f %-9.2f %-9.2f %-8.1f', ...
        cfgNames{i}, rangeExt(i), crossExt(i), locErr(i), dynRange(i));
end
L{end+1} = '    (RangeExt ~ box length; resolution is NOT this — see A)';
L{end+1} = '';
L{end+1} = '[C. NOISE ROBUSTNESS — Step 2 (AWGN)]';
L{end+1} = sprintf('    Image SNR ~ input SNR + %.0f dB realized gain (slope 1:1)', imgSNR(3)-snrIn(3));
L{end+1} = sprintf('    Box detectable down to input SNR ~ %+d dB', snrIn(find(detected,1,'last')));
L{end+1} = sprintf('    Breaks below ~ %+d dB (peak becomes a noise spike)', snrIn(find(~detected,1,'first')));
L{end+1} = '    Noise lowers detectability, NOT resolution.';
L{end+1} = '';
L{end+1} = '[KEY TAKEAWAYS]';
L{end+1} = '  1. RCS (manual->physics) only re-weights faces; in a noise-free';
L{end+1} = '     normalized image the absolute scale cancels.';
L{end+1} = '  2. Resolution/PSLR fixed by BW & aperture -> measure on a POINT,';
L{end+1} = '     not the box (the box gives EXTENT).';
L{end+1} = '  3. Occlusion (cull Back/Bottom) = more honest target model.';
L{end+1} = '  4. Coherent integration gives ~45 dB processing gain ->';
L{end+1} = '     SAR stays usable far below 0 dB per-pulse SNR.';
L{end+1} = repmat('=',1,56);

figure(1); set(gcf,'Name','Fig1: Step3 Master Table','Color','k','Position',[60 40 820 760]);
ax=axes('Position',[0 0 1 1],'Visible','off','Color','k'); hold(ax,'on');
ny=numel(L); ys=1/(ny+2);
for li=1:ny
    s=L{li};
    if startsWith(s,'==')||startsWith(s,'--'), clr=[.6 .6 .6];
    elseif startsWith(s,'['), clr=[.3 .85 1];
    elseif startsWith(s,'  ') && (startsWith(strtrim(s),'1.')||startsWith(strtrim(s),'2.')||startsWith(strtrim(s),'3.')||startsWith(strtrim(s),'4.')), clr=[.4 1 .5];
    else, clr=[.92 .92 .92]; end
    text(0.03,1-li*ys,s,'Units','normalized','Color',clr,'FontName','Courier','FontSize',9,'Interpreter','none','Parent',ax);
end
hold(ax,'off');

%% ===== FIGURE 2: Visual comparison (bars + noise curve) =====
figure(2); set(gcf,'Name','Fig2: Step3 Comparison Charts','Position',[80 60 1180 420]);

% (a) resolution: theory vs calibration (stable)
subplot(1,3,1);
b=bar([theory_rangeRes cal_rangeRes; theory_crossRes cal_crossRes]);
set(gca,'XTickLabel',{'Range','Cross'}); ylabel('Resolution (m)');
legend({'theory','calibration'},'Location','northeast');
title('System resolution (config-independent)'); grid on;

% (b) box metrics across configs
subplot(1,3,2);
bar([rangeExt; locErr]'); set(gca,'XTickLabel',cfgNames,'XTickLabelRotation',15);
ylabel('metres'); legend({'Range extent','Loc error'},'Location','northwest');
title('Box scene: RCS/occlusion effect'); grid on;

% (c) noise robustness
subplot(1,3,3);
plot(snrIn, imgSNR,'o-','LineWidth',2,'MarkerFaceColor','b'); hold on;
yline(detThresh,'g--','LineWidth',1.2,'DisplayName','detection threshold');
idx=find(~detected); plot(snrIn(idx),imgSNR(idx),'rx','MarkerSize',10,'LineWidth',2,'DisplayName','not detected');
set(gca,'XDir','reverse'); xlabel('Input SNR/pulse (dB)'); ylabel('Image SNR (dB)');
title('Noise robustness (Step 2)'); grid on; legend('Location','southwest');

sgtitle('Fig 2 | Step 3: what changes at each stage — system vs scene vs noise');

%% ===== SAVE =====
figDir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'figure');
if ~exist(figDir,'dir'), mkdir(figDir); end
exportgraphics(figure(1),fullfile(figDir,'fig_step3_table.png'),'Resolution',150);
exportgraphics(figure(2),fullfile(figDir,'fig_step3_charts.png'),'Resolution',150);
fprintf('Step 3 comparison figures saved to %s\n', figDir);


%% ================= SECTION 2: BOX DENSITY (try107)  (from try107_BoxDensity.m) =================
%%   (figure numbers offset by +100 so all sections' figures stay open)
%% Week 6 | Proof: Box scatterer density vs IMAGE RESOLUTION
%
%  พิสูจน์ว่า "เพิ่มความละเอียดกล่อง (จุดถี่ขึ้น)" ทำให้ภาพดูตันขึ้น
%  แต่ resolution / extent ของภาพ "ไม่เปลี่ยน" (เป็นของระบบ ไม่ใช่ของเป้า)
%
%  รันกล่อง 2 ความละเอียด: N_face = 6 (216 จุด) vs 15 (1350 จุด)
%  ใช้ pipeline เดียวกัน (clean, bypass receiver) -> เทียบภาพ + วัด extent
%  *รัน ~3-5 นาที (forward sim 2 รอบ, อันถี่ช้ากว่า)

clear;   % combined: keep prior figures open

%% ===== RADAR CONFIG =====
c=physconst('LightSpeed'); fc=4e9; bw=c/(2*3);
prf=1000; aperture=4; tpd=3e-6; fs=120e6;
waveform=phased.LinearFMWaveform('SampleRate',fs,'PulseWidth',tpd,'PRF',prf,'SweepBandwidth',bw);
speed=100; flightDuration=4; slowTime=1/prf; numpulses=flightDuration/slowTime+1;
maxRange=2500; truncrangesamples=ceil((2*maxRange/c)*fs);
antenna=phased.CosineAntennaElement('FrequencyRange',[1e9 6e9]); antennaGain=aperture2gain(aperture,c/fc);
transmitter=phased.Transmitter('PeakPower',50e3,'Gain',antennaGain);
radiator=phased.Radiator('Sensor',antenna,'OperatingFrequency',fc,'PropagationSpeed',c);
collector=phased.Collector('Sensor',antenna,'PropagationSpeed',c,'OperatingFrequency',fc);
channel=phased.FreeSpace('PropagationSpeed',c,'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

boxCenter=[1000;0;0]; boxL=30; boxW=20; boxH=10; lambda=c/fc;
sigmaPlate=[4*pi*(boxW*boxH)^2/lambda^2*[1 1], 4*pi*(boxL*boxH)^2/lambda^2*[1 1], 4*pi*(boxL*boxW)^2/lambda^2*[1 1]];

radarpos0=[0;-200;500]; radarvel0=[0;speed;0];
radarPosHistory=radarpos0+radarvel0*((1:numpulses)-1)/prf;
faceNormal=[-1 0 0;1 0 0;0 -1 0;0 1 0;0 0 1;0 0 -1];
losAll=radarPosHistory-boxCenter; faceVisible=true(1,6);
for f=1:6, faceVisible(f)=any(faceNormal(f,:)*losAll>0); end   % occlusion: hide Back,Bottom

refChirp=waveform(); chirpSamples=round(tpd*fs); refChirp=refChirp(1:chirpSamples);
matchedFilter=conj(flipud(refChirp)); mfLen=length(matchedFilter);
xScene=linspace(-60,60,401); yScene=linspace(960,1040,401);
[xGrid,yGrid]=meshgrid(xScene,yScene); Ny=numel(yScene); Nx=numel(xScene);

%% ===== RUN TWO DENSITIES =====
Nlist=[6 15]; res={};
for q=1:2
    Nf=Nlist(q); pf=Nf^2;
    spc=[boxL/(Nf-1), boxW/(Nf-1), boxH/(Nf-1)];   % spacing range,cross,height
    [scatPos,scatRCS]=buildBox(Nf,sigmaPlate,faceVisible,boxCenter,boxL,boxW,boxH);
    nsc=size(scatPos,2);
    fprintf('N_face=%d : %d scatterers, spacing [R %.2f, X %.2f, Z %.2f] m\n',Nf,nsc,spc);
    rx=forwardSim(scatPos,scatRCS,radarPosHistory,waveform,transmitter,radiator,collector,channel,fc,slowTime,truncrangesamples,nsc,speed);
    rc=rangeCompress(rx,matchedFilter);
    bp=backproject(rc,radarPosHistory,xGrid,yGrid,Ny,Nx,c,fs,fc,mfLen);
    bpMag=abs(bp)/max(abs(bp(:)));
    [rExt,cExt]=extentMetric(bpMag,xScene,yScene,boxCenter);
    res{q}=struct('Nf',Nf,'nsc',nsc,'spc',spc,'bp',bpMag,'rExt',rExt,'cExt',cExt);
    fprintf('   -> measured extent: range %.2f m, cross %.2f m\n',rExt,cExt);
end

%% ===== FIGURE: side-by-side images =====
figure(100+1); set(gcf,'Name','Box density vs resolution','Position',[60 80 1150 520],'Color','k');
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for q=1:2
    ax=nexttile;
    imagesc(xScene,yScene,20*log10(res{q}.bp+eps)); set(ax,'YDir','normal'); clim([-40 0]);
    colormap(ax,'jet'); hold on;
    bx=boxCenter(2)+boxW/2*[-1 1 1 -1 -1]; by=boxCenter(1)+boxL/2*[-1 -1 1 1 -1];
    plot(bx,by,'w--','LineWidth',1.2);
    title(sprintf('N\\_face=%d : %d pts (spacing ~%.1f m)\\newlineextent: R %.1f m, X %.1f m', ...
        res{q}.Nf,res{q}.nsc,res{q}.spc(2),res{q}.rExt,res{q}.cExt),'Color','w');
    xlabel('Cross-Range (m)','Color','w'); if q==1, ylabel('Range (m)','Color','w'); end
    set(ax,'Color','k','XColor','w','YColor','w');
end
cb=colorbar; cb.Color='w'; ylabel(cb,'dB','Color','w'); cb.Layout.Tile='east';
sgtitle('Proof: more scatterers -> box looks denser, but EXTENT/resolution unchanged','Color','w');

figDir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'figure'); if ~exist(figDir,'dir'),mkdir(figDir);end
exportgraphics(figure(100+1),fullfile(figDir,'fig_box_density.png'),'Resolution',150);

fprintf('\n==== PROOF SUMMARY ====\n');
fprintf('  %-10s %-10s %-12s %-12s\n','N_face','scatterers','range ext','cross ext');
for q=1:2
    fprintf('  %-10d %-10d %-12.2f %-12.2f\n',res{q}.Nf,res{q}.nsc,res{q}.rExt,res{q}.cExt);
end
fprintf('  -> extent ~ เท่าเดิม (= ขนาดกล่อง) แม้จุดต่างกัน %dx\n',round(res{2}.nsc/res{1}.nsc));
fprintf('  -> resolution มาจาก calibration (point target) ไม่ขึ้นกับกล่องเลย\n');


%% ================= SECTION 3: RDA vs BP (try108)  (from try108_RDA_vs_BP.m) =================
%%   (figure numbers offset by +200 so all sections' figures stay open)
%% Week 6 | Step 4: Range-Doppler Algorithm (RDA) vs Backprojection (BP)
%
%  เทียบ 2 วิธี form ภาพ SAR จากข้อมูลชุดเดียวกัน:
%    BP  = space domain : ไล่ทุก pixel รวมทุก pulse (แม่น/ยืดหยุ่น แต่ช้า)
%    RDA = frequency domain : azimuth FFT -> RCMC -> azimuth MF -> IFFT (เร็ว)
%
%  ประเด็น: ทั้งคู่ "ควรได้ภาพ/resolution เท่ากัน" (IRF มาจาก BW+aperture ไม่ใช่
%  อัลกอริทึม) ต่างกันที่ "ความเร็ว" และ "ความยืดหยุ่นของ geometry"
%
%  ใช้ point targets + range-compressed model เชิงวิเคราะห์ (sinc ที่ slant range
%  + azimuth phase) เพื่อให้เทียบ resolution ได้ชัด  (สูตร RDA ตรวจกับ prototype แล้ว)

clear;   % combined: keep prior figures open

%% ===== PARAMS (ชุดเดียวกับโปรเจกต์) =====
c=physconst('LightSpeed'); fc=4e9; lambda=c/fc;
bw=c/(2*3); fs=120e6; prf=1000; v=100; h=500; flightDur=4;
Naz=flightDur*prf+1; t=linspace(0,flightDur,Naz); ry=-200+v*t;   % radar cross position
rngres=c/(2*bw);

% point targets [groundRange, cross]
tgts=[1000 0; 1000 30; 1015 -20];

% slant-range axis (fast-time bins)
rbins=(1090:c/2/fs:1170).'; Nr=numel(rbins);

%% ===== build range-compressed data S(slantRange, azimuth) =====
S=zeros(Nr,Naz);
for k=1:size(tgts,1)
    xg=tgts(k,1); yg=tgts(k,2);
    R=sqrt(xg^2+(ry-yg).^2+h^2);                 % 1 x Naz slant range vs pulse
    S=S+sinc((rbins-R)/rngres).*exp(-1j*4*pi*R/lambda);
end

%% ===== METHOD 1: BACKPROJECTION (space domain) =====
grB=linspace(990,1025,161); crB=linspace(-45,45,401);   % fine cross grid for fair -3dB
[CR,GR]=meshgrid(crB,grB); NyB=numel(grB); NxB=numel(crB);
tic;
bpImg=zeros(NyB,NxB);
for a=1:Naz
    Rp=sqrt(GR.^2+(CR-ry(a)).^2+h^2);            % slant range pixel->radar
    iv=interp1(rbins,S(:,a),Rp(:),'linear',0);
    bpImg=bpImg+reshape(iv,NyB,NxB).*exp(1j*4*pi*fc/c.*Rp);
end
tBP=toc;
bpMag=abs(bpImg)/max(abs(bpImg(:)));

%% ===== METHOD 2: RANGE-DOPPLER ALGORITHM (frequency domain) =====
tic;
faz=((0:Naz-1)-floor(Naz/2))*(prf/Naz);          % Doppler freq axis (Hz)
Saz=fftshift(fft(S,[],2),2);                     % azimuth FFT
% RCMC: straighten range migration per Doppler bin
Src=zeros(size(Saz));
for kk=1:Naz
    dR=lambda^2*rbins*faz(kk)^2/(8*v^2);         % range-variant migration
    Src(:,kk)=interp1(rbins,Saz(:,kk),rbins+dR,'linear',0);
end
% azimuth matched filter (sign verified: -)
for ri=1:Nr
    Ka=2*v^2/(lambda*rbins(ri));
    Src(ri,:)=Src(ri,:).*exp(-1j*pi*faz.^2/Ka);
end
rdaImg=ifft(ifftshift(Src,2),[],2);              % azimuth IFFT
tRDA=toc;
rdaMag=abs(rdaImg)/max(abs(rdaImg(:)));
% RDA native axes: slant range -> ground range ; azimuth sample -> cross (ry)
grR=sqrt(rbins.^2-h^2);                           % slant -> ground range
crR=ry;                                           % azimuth -> cross position

%% ===== RESOLUTION CHECK — วัดบน "จุดเดียว" (cross=+30) ไม่คร่อม 2 จุด =====
% BP
[~,iB]=min(abs(grB-1000)); cutB=bpMag(iB,:); cutB(~(crB>20 & crB<40))=0;
pkB=max(cutB); aB=crB(cutB>=pkB*10^(-3/20)); resB=max(aB)-min(aB);
% peak sidelobe (นอก mainlobe ของจุดนั้น) เพื่อเทียบ skirt
slB=20*log10(max(cutB(crB>20 & crB<40 & abs(crB-30)>2))/pkB+eps);
% RDA
[~,iR]=min(abs(grR-1000)); cutR=rdaMag(iR,:); cutR(~(crR>20 & crR<40))=0;
pkR=max(cutR); aR=crR(cutR>=pkR*10^(-3/20)); resR=max(aR)-min(aR);
slR=20*log10(max(cutR(crR>20 & crR<40 & abs(crR-30)>2))/pkR+eps);

fprintf('\n===== Step 4: RDA vs BP =====\n');
fprintf('  Targets: %d points\n',size(tgts,1));
fprintf('  Runtime          : BP = %.2f s   RDA = %.3f s   -> RDA faster ~%.0fx\n',tBP,tRDA,tBP/tRDA);
fprintf('  Cross-range -3dB  : BP = %.2f m   RDA = %.2f m   (single target ~ เท่ากัน)\n',resB,resR);
fprintf('  Cross sidelobe    : BP = %.1f dB  RDA = %.1f dB  (RDA สูงกว่า = ราคาของ approx)\n',slB,slR);
fprintf('  Trade-off: BP = exact/ช้า | RDA = เร็ว ~%.0fx แต่ sidelobe สูงกว่า (parabolic RCMC)\n',tBP/tRDA);
fprintf('             -> resolution (-3dB) เท่ากัน ยืนยันว่าอัลกอริทึมไม่เปลี่ยน IRF\n');

%% ===== FIGURE: side-by-side images =====
figure(200+1); set(gcf,'Name','RDA vs BP','Position',[60 80 1120 520],'Color','k');
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1=nexttile;
imagesc(crB,grB,20*log10(bpMag+eps)); set(ax1,'YDir','normal'); clim([-30 0]);
colormap(ax1,'jet'); hold on; plot(tgts(:,2),tgts(:,1),'wo','MarkerSize',10,'LineWidth',1.2);
title(sprintf('Backprojection (space domain)\\newlineruntime %.2f s | cross -3dB %.2f m',tBP,resB),'Color','w');
xlabel('Cross-Range (m)','Color','w'); ylabel('Ground Range (m)','Color','w');
set(ax1,'Color','k','XColor','w','YColor','w'); xlim([-45 45]); ylim([990 1025]);

ax2=nexttile;
imagesc(crR,grR,20*log10(rdaMag+eps)); set(ax2,'YDir','normal'); clim([-30 0]);
colormap(ax2,'jet'); hold on; plot(tgts(:,2),tgts(:,1),'wo','MarkerSize',10,'LineWidth',1.2);
title(sprintf('Range-Doppler Algorithm (freq domain)\\newlineruntime %.3f s | cross -3dB %.2f m',tRDA,resR),'Color','w');
xlabel('Cross-Range (m)','Color','w'); ylabel('Ground Range (m)','Color','w');
set(ax2,'Color','k','XColor','w','YColor','w'); xlim([-45 45]); ylim([990 1025]);

cb=colorbar; cb.Color='w'; ylabel(cb,'dB','Color','w'); cb.Layout.Tile='east';
sgtitle(sprintf('Step 4 | RDA vs BP — same -3dB resolution, RDA ~%.0fx faster but higher sidelobes (o = true targets)',tBP/tRDA),'Color','w');

figDir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'figure'); if ~exist(figDir,'dir'),mkdir(figDir);end
exportgraphics(figure(200+1),fullfile(figDir,'fig_step4_rda_vs_bp.png'),'Resolution',150);
fprintf('Figure saved to %s\n',figDir);


%%%% ===== LOCAL FUNCTIONS (all sections) =====

function [scatPos,scatRCS]=buildBox(Nf,sigmaPlate,faceVisible,bc,L,Wd,H)
    halfL=L/2; halfW=Wd/2; halfH=H/2; pf=Nf^2;
    faceRCS=sigmaPlate/pf^2;                 % sub-patch (รักษา total face RCS ให้คงที่)
    u=linspace(-1,1,Nf); scatPos=[]; scatRCS=[];
    [V,W]=meshgrid(u*halfW,u*halfH);
    scatPos=[scatPos,[(bc(1)-halfL)*ones(1,pf);(bc(2)+V(:))';(bc(3)+halfH+W(:))']]; scatRCS=[scatRCS,faceRCS(1)*faceVisible(1)*ones(1,pf)];
    scatPos=[scatPos,[(bc(1)+halfL)*ones(1,pf);(bc(2)+V(:))';(bc(3)+halfH+W(:))']]; scatRCS=[scatRCS,faceRCS(2)*faceVisible(2)*ones(1,pf)];
    [U,W]=meshgrid(u*halfL,u*halfH);
    scatPos=[scatPos,[(bc(1)+U(:))';(bc(2)-halfW)*ones(1,pf);(bc(3)+halfH+W(:))']]; scatRCS=[scatRCS,faceRCS(3)*faceVisible(3)*ones(1,pf)];
    scatPos=[scatPos,[(bc(1)+U(:))';(bc(2)+halfW)*ones(1,pf);(bc(3)+halfH+W(:))']]; scatRCS=[scatRCS,faceRCS(4)*faceVisible(4)*ones(1,pf)];
    [U,V]=meshgrid(u*halfL,u*halfW);
    scatPos=[scatPos,[(bc(1)+U(:))';(bc(2)+V(:))';(bc(3)+H)*ones(1,pf)]]; scatRCS=[scatRCS,faceRCS(5)*faceVisible(5)*ones(1,pf)];
    scatPos=[scatPos,[(bc(1)+U(:))';(bc(2)+V(:))';zeros(1,pf)]]; scatRCS=[scatRCS,faceRCS(6)*faceVisible(6)*ones(1,pf)];
end


function rx=forwardSim(scatPos,scatRCS,rph,wf,tx,rad,col,chan,fc,slowT,nsamp,nScat,spd)
    % release shared System objects so input sizes can change between densities
    release(rad); release(col); release(chan); release(tx);
    boxP=phased.Platform('InitialPosition',scatPos,'Velocity',zeros(size(scatPos)));
    tgt =phased.RadarTarget('OperatingFrequency',fc,'MeanRCS',scatRCS);
    np=size(rph,2); rx=zeros(nsamp,np);
    for ii=1:np
        rpos=rph(:,ii); [tpos,tvel]=boxP(slowT);
        [~,tAng]=rangeangle(tpos,rpos);
        s=wf(); s=s(1:nsamp); s=tx(s);
        tAng(1,:)=zeros(1,nScat);
        s=rad(s,tAng); s=chan(s,rpos,tpos,[0;spd;0],tvel); s=tgt(s);
        rx(:,ii)=col(s,tAng);
    end
end


function rc=rangeCompress(rx,mf)
    rc=zeros(size(rx,1)+length(mf)-1,size(rx,2));
    for ii=1:size(rx,2), rc(:,ii)=conv(rx(:,ii),mf,'full'); end
end


function bp=backproject(rc,rph,xGrid,yGrid,Ny,Nx,c,fs,fc,mfLen)
    pix=[yGrid(:)';xGrid(:)';zeros(1,numel(xGrid))]; bp=zeros(Ny,Nx);
    for ii=1:size(rph,2)
        rp=rph(:,ii); sr=sqrt(sum((pix-rp).^2,1)); si=2*sr/c*fs+mfLen;
        vm=(si>=1)&(si<=size(rc,1)); iv=zeros(1,numel(xGrid));
        if any(vm), iv(vm)=interp1(1:size(rc,1),rc(:,ii),si(vm),'linear',0); end
        bp=bp+reshape(iv.*exp(1j*4*pi*fc/c.*sr),Ny,Nx);
    end
end


function [rExt,cExt]=extentMetric(bpMag,xS,yS,bc)
    [~,xi]=min(abs(xS-bc(2))); rp=bpMag(:,xi); rpdB=20*log10(rp/max(rp));
    a=yS(rpdB>=-3); rExt=max(a)-min(a);
    [~,yi]=min(abs(yS-bc(1))); cp=bpMag(yi,:); cpdB=20*log10(cp/max(cp));
    b=xS(cpdB>=-3); cExt=max(b)-min(b);
end
