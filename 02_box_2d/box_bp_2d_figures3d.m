%% SAR Backprojection — 3D Box Target | ALL FIGURES IN 3D
%  หมายเหตุ: บล็อก RADAR CONFIGURATION สืบทอดมาจาก 01_point_target/point3_bp_2d.m
%  ซึ่งดัดแปลงจากตัวอย่าง Stripmap SAR ของ MathWorks (ดู attribution ในไฟล์นั้น)
%  path ทั้งหมดอิงจาก root ของ repo — ย้ายโฟลเดอร์ทั้งก้อนได้ ไม่ต้องแก้โค้ด
%  ทุก Figure ถูก upgrade เป็น 3D visualization
%  เลือก method ที่เหมาะกับแต่ละ Figure:
%
%  Fig 1  — scatter3   : Box geometry + scatterers แยกตามหน้า + solid box faces
%  Fig 2  — waterfall  : Raw echo (surf 3D waterfall per pulse)
%  Fig 3  — surf       : Raw echo magnitude เป็น 3D surface (range × pulse × amp)
%  Fig 4  — surf       : Range-compressed data เป็น 3D surface
%  Fig 5  — surf       : BP Image linear scale เป็น 3D surface + box wireframe
%  Fig 6  — surf+contour: BP Image dB scale เป็น 3D surface พร้อม contour ด้านล่าง
%  Fig 7  — plot3      : Range profile เป็น 3D ribbon ผ่าน scene
%  Fig 8  — plot3      : Cross-range profiles 3 slice ใน 3D space
%  Fig 9  — scatter3+patch: Box geometry สมบูรณ์ + radar flight path
%  Fig 10 — isosurface : BP Image volume rendering (NEW)
%  Fig 11 — SAR Scene  : Full 3D scene overview (radar + box + beam footprint)

clear; clc; close all;

%% ===== RADAR CONFIGURATION =====

c  = physconst('LightSpeed');
fc = 4e9;
rangeResolution = 3;
bw = c / (2*rangeResolution);

prf   = 1000;
aperture = 4;
tpd   = 3e-6;
fs    = 120e6;

waveform = phased.LinearFMWaveform( ...
    'SampleRate', fs, 'PulseWidth', tpd, ...
    'PRF', prf, 'SweepBandwidth', bw);

speed          = 100;
flightDuration = 4;

radarPlatform = phased.Platform( ...
    'InitialPosition', [0; -200; 500], 'Velocity', [0; speed; 0]);

slowTime   = 1/prf;
numpulses  = flightDuration/slowTime + 1;

maxRange          = 2500;
truncrangesamples = ceil((2*maxRange/c)*fs);
fastTime          = (0:1/fs:(truncrangesamples-1)/fs);

antenna     = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);
antennaGain = aperture2gain(aperture, c/fc);
transmitter = phased.Transmitter('PeakPower', 50e3, 'Gain', antennaGain);
radiator    = phased.Radiator('Sensor', antenna, 'OperatingFrequency', fc, 'PropagationSpeed', c);
collector   = phased.Collector('Sensor', antenna, 'PropagationSpeed', c, 'OperatingFrequency', fc);
receiver    = phased.ReceiverPreamp('SampleRate', fs, 'NoiseFigure', 30);
channel     = phased.FreeSpace('PropagationSpeed', c, 'OperatingFrequency', fc, ...
                               'SampleRate', fs, 'TwoWayPropagation', true);

%% ===== 3D BOX DEFINITION =====

boxCenter = [1000; 0; 0];
boxL = 30;   halfL = boxL/2;
boxW = 20;   halfW = boxW/2;
boxH = 10;   halfH = boxH/2;

N_face = 6;
u = linspace(-1, 1, N_face);

scatPos = [];
scatRCS = [];
faceColors = [1 0.2 0.2; 0.2 0.4 1; 0.2 0.8 0.2; 0.8 0.2 0.8; 0.2 0.8 0.8; 1 0.8 0.2];
faceLabels = {'Front','Back','Left','Right','Top','Bottom'};
faceRCS    = [2.0, 0.5, 1.0, 1.0, 0.8, 0.3];

[V, W] = meshgrid(u*halfW, u*halfH);

% Face 1: FRONT
xF = (boxCenter(1)-halfL)*ones(size(V));
scatPos = [scatPos, [xF(:)'; (boxCenter(2)+V(:))'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(1)*ones(1,N_face^2)];

% Face 2: BACK
xF = (boxCenter(1)+halfL)*ones(size(V));
scatPos = [scatPos, [xF(:)'; (boxCenter(2)+V(:))'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(2)*ones(1,N_face^2)];

[U, W] = meshgrid(u*halfL, u*halfH);

% Face 3: LEFT
yF = (boxCenter(2)-halfW)*ones(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; yF(:)'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(3)*ones(1,N_face^2)];

% Face 4: RIGHT
yF = (boxCenter(2)+halfW)*ones(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; yF(:)'; (boxCenter(3)+halfH+W(:))']];
scatRCS = [scatRCS, faceRCS(4)*ones(1,N_face^2)];

[U, V] = meshgrid(u*halfL, u*halfW);

% Face 5: TOP
zF = (boxCenter(3)+boxH)*ones(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; (boxCenter(2)+V(:))'; zF(:)']];
scatRCS = [scatRCS, faceRCS(5)*ones(1,N_face^2)];

% Face 6: BOTTOM
zF = zeros(size(U));
scatPos = [scatPos, [(boxCenter(1)+U(:))'; (boxCenter(2)+V(:))'; zF(:)']];
scatRCS = [scatRCS, faceRCS(6)*ones(1,N_face^2)];

numScatterers = size(scatPos, 2);
ptsPerFace = N_face^2;
fprintf('Box: %d scatterers total\n', numScatterers);

%% ===== PLATFORM HISTORY =====

radarpos0 = [0; -200; 500];
radarvel0 = [0; speed; 0];
pulseIndex = 1:numpulses;
timeAxis   = (pulseIndex-1) / prf;
radarPosHistory = radarpos0 + radarvel0 * timeAxis;

%% ===== SIMULATE RAW ECHO =====

fprintf('Simulating raw echo...\n');
rxsig = zeros(truncrangesamples, numpulses);

scatVel      = zeros(size(scatPos));
boxPlatform2 = phased.Platform('InitialPosition', scatPos, 'Velocity', scatVel);
boxTarget    = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', scatRCS);

radarPlatform2 = phased.Platform('InitialPosition', [0;-200;500], 'Velocity', [0;speed;0]);

for ii = 1:numpulses
    [radarpos, radarvel] = radarPlatform2(slowTime);
    [tpos, tvel]         = boxPlatform2(slowTime);
    [~, targetAngle]     = rangeangle(tpos, radarpos);
    sig = waveform();
    sig = sig(1:truncrangesamples);
    sig = transmitter(sig);
    targetAngle(1,:) = zeros(1, numScatterers);
    sig = radiator(sig, targetAngle);
    sig = channel(sig, radarpos, tpos, radarvel, tvel);
    sig = boxTarget(sig);
    sig = collector(sig, targetAngle);
    rxsig(:,ii) = receiver(sig);
end
fprintf('Raw echo done.\n');

rangeAxis = fastTime * c / 2;

%% ===== RANGE COMPRESSION =====

refPulse      = waveform();
chirpSamples  = round(tpd * fs);
refChirp      = refPulse(1:chirpSamples);
matchedFilter = conj(flipud(refChirp));

rxsigRC = zeros(size(rxsig,1)+length(matchedFilter)-1, numpulses);
for ii = 1:numpulses
    rxsigRC(:,ii) = conv(rxsig(:,ii), matchedFilter, 'full');
end

rangeAxisRC = ((0:size(rxsigRC,1)-1).' - (length(matchedFilter)-1)) / fs * c/2;
rxsigRCmag  = abs(rxsigRC);
rxsigRCmag  = rxsigRCmag ./ max(rxsigRCmag(:));

%% ===== BACKPROJECTION =====

fprintf('Running Backprojection...\n');

xScene = linspace(-100, 100, 401);
yScene = linspace(850,  1200, 401);
[xGrid, yGrid] = meshgrid(xScene, yScene);
Ny = length(yScene);  Nx = length(xScene);
pixelPos = [yGrid(:)'; xGrid(:)'; zeros(1,numel(xGrid))];

bpImage = zeros(Ny, Nx);
for ii = 1:numpulses
    rp         = radarPosHistory(:, ii);
    diffVec    = pixelPos - rp;
    slantRange = sqrt(sum(diffVec.^2, 1));
    sampleIdx  = 2 .* slantRange ./ c .* fs + length(matchedFilter);
    validMask  = (sampleIdx >= 1) & (sampleIdx <= size(rxsigRC,1));
    interpVals = zeros(1, numel(xGrid));
    if any(validMask)
        interpVals(validMask) = interp1(1:size(rxsigRC,1), rxsigRC(:,ii), ...
                                        sampleIdx(validMask), 'linear', 0);
    end
    phaseCorr    = exp(1j*4*pi*fc/c .* slantRange);
    bpImage      = bpImage + reshape(interpVals .* phaseCorr, Ny, Nx);
end
fprintf('Backprojection complete.\n');

bpImageMag = abs(bpImage);
bpImageMag = bpImageMag ./ max(bpImageMag(:));
bpImagedB  = 20*log10(bpImageMag + eps);

% box outline corners
boxX_cr = boxCenter(2) + halfW*[-1 1 1 -1 -1];
boxY_r  = boxCenter(1) + halfL*[-1 -1 1 1 -1];
boxZ_ln = zeros(1,5);

%% ===== FIGURE 1: 3D Box Geometry =====

figure(1);
set(gcf,'Name','Fig1: 3D Box Geometry');

[Verts, Faces] = makeBoxPatch(boxCenter(2), boxCenter(1), boxCenter(3), halfW, halfL, boxH);
patch('Vertices', Verts, 'Faces', Faces, ...
      'FaceColor', [0.7 0.85 1.0], 'FaceAlpha', 0.15, ...
      'EdgeColor', [0.3 0.3 0.8], 'LineWidth', 1.5);
hold on;
for f = 1:6
    idx = (f-1)*ptsPerFace + (1:ptsPerFace);
    scatter3(scatPos(2,idx), scatPos(1,idx), scatPos(3,idx), ...
        faceRCS(f)*25 + 10, ...
        repmat(faceColors(f,:), ptsPerFace, 1), 'filled', ...
        'DisplayName', sprintf('%s (RCS=%.1f)', faceLabels{f}, faceRCS(f)));
end
xlabel('Cross-Range (m)'); ylabel('Range (m)'); zlabel('Height (m)');
title('Fig 1 | 3D Box Target — Scatterer Geometry');
legend('Box (transparent)', faceLabels{:}, 'Location','northeast');
grid on; view(40, 25);
xlim([-60 60]); ylim([920 1080]); zlim([-5 35]);
hold off;

%% ===== FIGURE 2: SAR Acquisition Scene 3D =====

figure(2);
set(gcf,'Name','Fig2: SAR Acquisition Scene 3D');

bpImg2dB = max(bpImagedB, -30);
surf(xScene, yScene, zeros(Ny,Nx), bpImg2dB, 'EdgeAlpha',0, 'FaceAlpha',0.85);
colormap(gca,'jet'); hold on;

[Verts2, Faces2] = makeBoxPatch(boxCenter(2),boxCenter(1),boxCenter(3),halfW,halfL,boxH);
patch('Vertices',Verts2,'Faces',Faces2, ...
    'FaceColor',[0.3 0.6 1],'FaceAlpha',0.3,'EdgeColor',[1 0.9 0.1],'LineWidth',3, ...
    'DisplayName','Box');
scatter3(scatPos(2,:), scatPos(1,:), scatPos(3,:), ...
    scatRCS*15+5, 'w', 'filled', 'MarkerFaceAlpha',0.5, 'HandleVisibility','off');

rph2 = radarPosHistory;
plot3(rph2(2,1:50:end), rph2(1,1:50:end), rph2(3,1:50:end), ...
    'b-', 'LineWidth', 2.5, 'DisplayName','Radar path');
scatter3(rph2(2,1),   rph2(1,1),   rph2(3,1),   200, 'gs','filled','DisplayName','Start');
scatter3(rph2(2,end), rph2(1,end), rph2(3,end), 200, 'rs','filled','DisplayName','End');

midP2 = round(numpulses/2);
rp2   = rph2(:, midP2);
fpX2  = boxCenter(2) + halfW*2*[-1 1 1 -1];
fpY2  = boxCenter(1) + halfL*2*[-1 -1 1 1];
for cr = 1:4
    hb = plot3([rp2(2) fpX2(cr)],[rp2(1) fpY2(cr)],[rp2(3) 1],'c-','LineWidth',1.2);
    try hb.Color(4)=0.45; catch; hb.Color=[0 0.85 0.85]; end 
end
text(rp2(2)+5, rp2(1), rp2(3)+15,'Radar (mid)','Color',[0 0.3 1],'FontWeight','bold','FontSize',9);
text(boxCenter(2)+halfW+3, boxCenter(1), boxH+25,'BOX','Color','w','FontWeight','bold','FontSize',11);

cb2 = colorbar; clim([-30 0]); ylabel(cb2,'dB');
xlabel('Cross-Range (m)'); ylabel('Range (m)'); zlabel('Height (m)');
title('Fig 2 | SAR Scene — Radar Path + Beam + Box + BP Image (dB)');
legend({'BP Image','Box','Radar path','Start','End'},'Location','northwest');
grid on;
xlim([-150 150]); ylim([-280 1100]); zlim([-5 560]);
view(-20, 28); camlight left; lighting gouraud;
hold off;

%% ===== FIGURE 3: Raw Echo — imagesc (linear amplitude) =====

figure(3);
set(gcf,'Name','Fig3: Raw Echo imagesc (linear)');

inRange2  = (rangeAxis >= 1060) & (rangeAxis <= 1170);
pulseAxis = 1:numpulses;
radarY    = radarpos0(2) + radarvel0(2)*(pulseAxis-1)/prf;
expectedRange = sqrt((boxCenter(1)-radarpos0(1))^2 + (boxCenter(2)-radarY).^2 + radarpos0(3)^2);

imagesc(1:numpulses, rangeAxis(inRange2), abs(rxsig(inRange2, :)));
set(gca,'YDir','normal');
colormap(gca,'parula'); colorbar;
xlabel('Pulse Index (slow time)'); ylabel('Slant Range (m)');
title('Fig 3 | Raw Echo — Range × Pulse (linear amplitude)');
hold on;
plot(pulseAxis, expectedRange, 'r--', 'LineWidth', 1.5, 'DisplayName','Expected slant range (box center)');
ylim([1060 1170]);
legend('Location','northeast'); grid on;
hold off;

%% ===== FIGURE 4: Raw Echo — imagesc (dB scale) =====

figure(4);
set(gcf,'Name','Fig4: Raw Echo imagesc (dB)');

rawdB = 20*log10(abs(rxsig(inRange2, :)) + eps);
rawdB = rawdB - max(rawdB(:));
imagesc(1:numpulses, rangeAxis(inRange2), rawdB);
set(gca,'YDir','normal');
clim([-40 0]);
colormap(gca,'jet'); cb4 = colorbar; ylabel(cb4,'dB');
xlabel('Pulse Index (slow time)'); ylabel('Slant Range (m)');
title('Fig 4 | Raw Echo — Range × Pulse (dB, clip -40 dB)');
hold on;
plot(pulseAxis, expectedRange, 'w--', 'LineWidth', 1.5, 'DisplayName','Expected slant range (box center)');
ylim([1060 1170]);
legend('Location','northeast'); grid on;
hold off;

%% ===== FIGURE 5: Range Compressed — imagesc (dB) =====

figure(5);
set(gcf,'Name','Fig5: Range Compressed');

inR5    = (rangeAxisRC >= 1060) & (rangeAxisRC <= 1170);
step5   = 10;
pulses5 = 1:step5:numpulses;
imagesc(pulses5, rangeAxisRC(inR5), rxsigRCmag(inR5, pulses5));
set(gca,'YDir','normal');
colormap(gca,'jet'); cb5 = colorbar; ylabel(cb5,'Normalized Amp');
xlabel('Pulse Index'); ylabel('Slant Range (m)');
title('Fig 5 | Range-Compressed Data — peak คมขึ้นหลัง matched filter');
ylim([1060 1170]); grid on;

%% ===== FIGURE 6: BP Image — 3D Surface + 2D Views =====

figure(6);
set(gcf,'Name','Fig6: BP Image 3D+2D');
set(gcf,'Color','k');

bpClip = max(bpImagedB, -40);

tl = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');

% --- Left column: 3D Surface (span both rows) ---
ax1 = nexttile([2 1]);
surf(xScene, yScene, bpClip, 'EdgeAlpha', 0);
shading interp; colormap(ax1,'jet'); clim([-40 0]);
hold on;
contour3(xScene, yScene, bpClip, 12, 'k', 'LineWidth', 0.4);
plot3(boxX_cr, boxY_r, 2*ones(1,5), 'y--', 'LineWidth', 2);
xlabel('Cross-Range (m)','Color','w'); ylabel('Range (m)','Color','w'); zlabel('dB','Color','w');
title('3D Surface (dB)', 'Color','w');
colorbar; view(-40, 30); grid on; camlight left; lighting gouraud;
set(ax1,'Color','k','XColor','w','YColor','w','ZColor','w'); hold off;

% --- Top-right: Top-down imagesc ---
ax2 = nexttile;
imagesc(xScene, yScene, bpClip); axis xy;
colormap(ax2,'jet'); clim([-40 0]);
hold on;
plot(boxX_cr, boxY_r, 'y--', 'LineWidth', 1.8);
hold off;
xlabel('Cross-Range (m)','Color','w'); ylabel('Range (m)','Color','w');
title('Top-down View', 'Color','w');
colorbar; set(ax2,'Color','k','XColor','w','YColor','w');

% --- Bottom-right: Range Profile at cross-range = 0 ---
ax3 = nexttile;
[~, xMid] = min(abs(xScene - boxCenter(2)));
rangeProf = bpClip(:, xMid);
plot(rangeProf, yScene, 'c-', 'LineWidth', 1.8); hold on;
xline(boxCenter(1)-halfL, 'y--', 'Front', 'LineWidth',1.2);
xline(boxCenter(1)+halfL, 'y--', 'Back',  'LineWidth',1.2);
xlabel('dB','Color','w'); ylabel('Range (m)','Color','w');
title('Range Profile (cross=0)', 'Color','w');
grid on; set(ax3,'Color','k','XColor','w','YColor','w'); hold off;

sgtitle('Fig 6 | BP Image', 'Color','w', 'FontSize',11);

%% ===== FIGURE 7: Cross-Range Profiles — 3D Filled Slices =====

figure(7);
set(gcf,'Name','Fig7: Cross-Range Profiles 3D');

rangeSlices = [boxCenter(1)-halfL, boxCenter(1), boxCenter(1)+halfL];
sliceLabels = {'Front face','Box center','Back face'};
sliceColors3 = {[1 0.2 0.2], [0.1 0.1 0.1], [0.2 0.7 0.2]};

hold on;
for k = 1:3
    [~, rIdx] = min(abs(yScene - rangeSlices(k)));
    prof = bpImageMag(rIdx, :);
    % plot3: x=cross-range, y=fixed range, z=amplitude
    fill3(xScene([1:end end 1]), ...
          rangeSlices(k)*ones(1,Nx+2), ...
          [prof 0 0], ...
          sliceColors3{k}, 'FaceAlpha', 0.3, 'EdgeColor', sliceColors3{k}, ...
          'DisplayName', sliceLabels{k});
end

% วาดเส้น box boundary ใน 3D
plot3([boxX_cr; boxX_cr], [boxY_r; boxY_r], ...
      [zeros(1,5); 0.8*ones(1,5)], 'k--', 'LineWidth', 1);

xlabel('Cross-Range (m)'); ylabel('Ground Range (m)'); zlabel('Amplitude');
title('Fig 7 | Cross-Range Profiles — 3D Filled Slices');
legend show; grid on; view(-30, 25);
hold off;

%% ===== FIGURE 12: SAR Image 2D — Top-down (dB) =====
%  ภาพ SAR มาตรฐานที่คนในสายงานคาดหวัง

figure(8);
set(gcf,'Name','Fig8: SAR Image 2D (dB)');

imagesc(xScene, yScene, bpImagedB);
set(gca,'YDir','normal');
clim([-40 0]);
colormap(gca,'jet');
cb12 = colorbar; ylabel(cb12,'dB');
hold on;

% วาด box outline สีขาว
plot(boxX_cr, boxY_r, 'w--', 'LineWidth', 2, 'DisplayName','Box outline');
% cross-hair ที่ center
plot([0 0], [yScene(1) yScene(end)], 'w:', 'LineWidth', 1, 'HandleVisibility','off');
plot([xScene(1) xScene(end)], [boxCenter(1) boxCenter(1)], 'w:', 'LineWidth', 1, 'HandleVisibility','off');

xlabel('Cross-Range (m)'); ylabel('Ground Range (m)');
title('Fig 8 | SAR BP Image — Top-down 2D (dB, clip -40 dB)');
legend('Box outline','Location','northeast');
axis equal tight; grid on;
hold off;

%% ===== FIGURE 13: Range Profile 2D — dB with -3dB width =====
%  วัด range resolution: peak width ที่ -3 dB

figure(9);
set(gcf,'Name','Fig9: Range Profile 2D (dB)');

[~, xZeroIdx13] = min(abs(xScene));
rProf13    = bpImageMag(:, xZeroIdx13);
rProf13dB  = 20*log10(rProf13 + eps);
rProf13dB  = rProf13dB - max(rProf13dB);   % normalize peak to 0 dB

plot(yScene, rProf13dB, 'b-', 'LineWidth', 1.8); hold on;
yline(-3,  'r--', 'LineWidth', 1.5, 'DisplayName','-3 dB');
yline(-10, 'm--', 'LineWidth', 1,   'DisplayName','-10 dB');

% box face markers
faceR13  = [boxCenter(1)-halfL, boxCenter(1), boxCenter(1)+halfL];
fClr13   = {'r','k','g'};
fName13  = {'Front (985m)','Center (1000m)','Back (1015m)'};
for k = 1:3
    xline(faceR13(k), fClr13{k}, fName13{k}, 'LineWidth', 1.5, 'LabelVerticalAlignment','bottom');
end

% คำนวณ -3 dB width
above3dB = yScene(rProf13dB >= -3);
if ~isempty(above3dB)
    w3 = max(above3dB) - min(above3dB);
    xline(min(above3dB), 'r:', 'LineWidth', 1, 'HandleVisibility','off');
    xline(max(above3dB), 'r:', 'LineWidth', 1, 'HandleVisibility','off');
    title(sprintf('Fig 9 | Range Profile (Cross-Range=0) | -3dB width = %.1f m  (theory ≈ %.0f m)', w3, rangeResolution));
else
    title('Fig 9 | Range Profile (Cross-Range=0)');
end

xlim([850 1200]); ylim([-40 5]);
xlabel('Ground Range (m)'); ylabel('Normalized (dB)');
legend('Location','northeast'); grid on;
hold off;

%% ===== FIGURE 14: Cross-Range Profile 2D — dB with -3dB width =====
%  วัด azimuth/cross-range resolution

figure(10);
set(gcf,'Name','Fig10: Cross-Range Profile 2D (dB)');

[~, rCenterIdx] = min(abs(yScene - boxCenter(1)));
crProf14   = bpImageMag(rCenterIdx, :);
crProf14dB = 20*log10(crProf14 + eps);
crProf14dB = crProf14dB - max(crProf14dB);   % normalize peak to 0 dB

plot(xScene, crProf14dB, 'b-', 'LineWidth', 1.8); hold on;
yline(-3,  'r--', 'LineWidth', 1.5, 'DisplayName','-3 dB');
yline(-10, 'm--', 'LineWidth', 1,   'DisplayName','-10 dB');

% box cross-range edges
xline(-halfW, 'k--', 'LineWidth', 1.5, 'DisplayName','Box edges (±10m)');
xline(+halfW, 'k--', 'LineWidth', 1.5, 'HandleVisibility','off');

% คำนวณ -3 dB cross-range width
above3cr = xScene(crProf14dB >= -3);
if ~isempty(above3cr)
    w3cr = max(above3cr) - min(above3cr);
    xline(min(above3cr), 'r:', 'LineWidth', 1, 'HandleVisibility','off');
    xline(max(above3cr), 'r:', 'LineWidth', 1, 'HandleVisibility','off');
    % theoretical cross-range resolution: λ*R / (2*aperture_length)
    lambda = c/fc;
    apertureLen = speed * flightDuration;
    crRes_theory = lambda * boxCenter(1) / (2 * apertureLen);
    title(sprintf('Fig 10 | Cross-Range Profile (Range=1000m) | -3dB width = %.1f m  (theory ≈ %.2f m)', w3cr, crRes_theory));
else
    title('Fig 10 | Cross-Range Profile (Range=1000m)');
end

xlim([-100 100]); ylim([-40 5]);
xlabel('Cross-Range (m)'); ylabel('Normalized (dB)');
legend('Location','northeast'); grid on;
hold off;

%% ===== SAVE FIGURES =====

figDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figure');
if ~exist(figDir, 'dir'); mkdir(figDir); end
for fn = 1:10
    fh = figure(fn);
    exportgraphics(fh, fullfile(figDir, sprintf('fig%02d.png', fn)), 'Resolution', 150);
end
fprintf('Figures saved to %s\n', figDir);

%% ===== EVALUATION =====

fprintf('\n========== EVALUATION — SAR Backprojection Box Target ==========\n');

%% -- 1. Radar & Scene Parameters --
lambda      = c / fc;
apertureLen = speed * flightDuration;   % synthetic aperture length (m)
fprintf('\n[Radar Parameters]\n');
fprintf('  Center frequency : %.1f GHz\n', fc/1e9);
fprintf('  Bandwidth        : %.1f MHz\n', bw/1e6);
fprintf('  Wavelength       : %.4f m\n', lambda);
fprintf('  PRF              : %d Hz\n', prf);
fprintf('  Flight duration  : %.1f s  |  Aperture length: %.0f m\n', flightDuration, apertureLen);
fprintf('  Platform speed   : %.0f m/s  |  Height: %.0f m\n', speed, radarpos0(3));
fprintf('\n[Box Target]\n');
fprintf('  Center (range, cross, height): (%.0f, %.0f, %.0f) m\n', boxCenter(1), boxCenter(2), boxCenter(3));
fprintf('  Size: %.0f m (range) x %.0f m (cross) x %.0f m (height)\n', boxL, boxW, boxH);
fprintf('  Scatterers: %d total (%d per face, 6 faces)\n', numScatterers, ptsPerFace);
fprintf('  Face RCS: Front=%.1f, Back=%.1f, Left=%.1f, Right=%.1f, Top=%.1f, Bottom=%.1f\n', faceRCS(:)');

%% -- 2. Theoretical Resolution --
rangeRes_theory  = c / (2 * bw);
crRes_theory     = lambda * boxCenter(1) / (2 * apertureLen);
fprintf('\n[Theoretical Resolution]\n');
fprintf('  Range resolution     : %.2f m  (c/2B)\n', rangeRes_theory);
fprintf('  Cross-range resolution: %.4f m  (λR/2L)\n', crRes_theory);
fprintf('  NOTE: cross-range res is very fine (%.2f cm) — scatterers may be individually resolved\n', crRes_theory*100);

%% -- 3. Measured Range Resolution (-3 dB) --
[~, xZeroIdx2] = min(abs(xScene));
rProf2   = bpImageMag(:, xZeroIdx2);
rProfDB2 = 20*log10(rProf2 + eps);
rProfDB2 = rProfDB2 - max(rProfDB2);

above3r = yScene(rProfDB2 >= -3);
above6r = yScene(rProfDB2 >= -6);

fprintf('\n[Measured Range Resolution (at Cross-Range=0)]\n');
if ~isempty(above3r)
    fprintf('  -3 dB width : %.2f m  (theory: %.2f m)\n', max(above3r)-min(above3r), rangeRes_theory);
    fprintf('  -3 dB extent: %.1f – %.1f m\n', min(above3r), max(above3r));
end
if ~isempty(above6r)
    fprintf('  -6 dB width : %.2f m  (box length = %.0f m)\n', max(above6r)-min(above6r), boxL);
end

% Peak sidelobe ratio (range) — ดูบริเวณนอก main lobe
[peakVal, peakIdx] = max(rProf2);
sidelobeRegion = rProf2;
% mask out ±(boxL/2 * 2) m around peak
maskRange = abs(yScene - yScene(peakIdx)) > boxL;
if sum(maskRange) > 0
    slPeak = max(sidelobeRegion(maskRange));
    pslr_range = 20*log10(slPeak / peakVal);
    fprintf('  PSLR (range) : %.1f dB  (good SAR: < -13 dB)\n', pslr_range);
end

%% -- 4. Measured Cross-Range Resolution (-3 dB) --
[~, rCenterIdx2] = min(abs(yScene - boxCenter(1)));
crProf2   = bpImageMag(rCenterIdx2, :);
crProfDB2 = 20*log10(crProf2 + eps);
crProfDB2 = crProfDB2 - max(crProfDB2);

above3cr = xScene(crProfDB2 >= -3);

fprintf('\n[Measured Cross-Range Resolution (at Range=%.0f m)]\n', boxCenter(1));
if ~isempty(above3cr)
    fprintf('  -3 dB width : %.2f m  (theory: %.4f m)\n', max(above3cr)-min(above3cr), crRes_theory);
    fprintf('  -3 dB extent: %.1f – %.1f m\n', min(above3cr), max(above3cr));
end

% PSLR cross-range
[crPeakVal, crPeakIdx] = max(crProf2);
maskCR = abs(xScene - xScene(crPeakIdx)) > halfW*2;
if sum(maskCR) > 0
    slCRpeak = max(crProf2(maskCR));
    pslr_cr = 20*log10(slCRpeak / crPeakVal);
    fprintf('  PSLR (cross-range): %.1f dB\n', pslr_cr);
end

%% -- 5. Peak SNR & Dynamic Range --
bpPeak = max(bpImageMag(:));
% noise floor = mean of corner regions (far from box)
cornerMask = (abs(xGrid - boxCenter(2)) > halfW*5) & (abs(yGrid - boxCenter(1)) > halfL*5);
noiseFloor = mean(bpImageMag(cornerMask));
snr_dB = 20*log10(bpPeak / noiseFloor);
fprintf('\n[Dynamic Range & SNR]\n');
fprintf('  BP image peak    : %.4f (normalized)\n', bpPeak);
fprintf('  Noise floor (est): %.4f\n', noiseFloor);
fprintf('  Peak SNR (est)   : %.1f dB\n', snr_dB);

%% -- 6. Box Localization --
% หา centroid ของ bright region (> 0.5 * peak)
brightMask = bpImageMag > 0.5;
if sum(brightMask(:)) > 0
    xCentroid = mean(xGrid(brightMask));
    yCentroid = mean(yGrid(brightMask));
    fprintf('\n[Box Localization]\n');
    fprintf('  True center   : (cross=%.0f, range=%.0f) m\n', boxCenter(2), boxCenter(1));
    fprintf('  BP centroid   : (cross=%.1f, range=%.1f) m\n', xCentroid, yCentroid);
    fprintf('  Location error: cross=%.2f m,  range=%.2f m\n', ...
        abs(xCentroid-boxCenter(2)), abs(yCentroid-boxCenter(1)));
end

fprintf('\nDone. Figures 1–10 (pipeline order):\n');
fprintf('  1=BoxGeom  2=Scene3D  3=RawLinear  4=RawdB  5=RC\n');
fprintf('  6=BPdB3D   7=CRslices  8=SARimage  9=RangeProf  10=CRProf\n');
fprintf('=============================================================\n');

%% ===== LOCAL FUNCTIONS =====

function [Verts, Faces] = makeBoxPatch(cx, cy, cz, hl, hw, hh)
% makeBoxPatch — สร้าง vertices และ faces ของกล่องสำหรับ patch()
%   cx = cross-range center (X-axis)
%   cy = range center       (Y-axis)
%   cz = height base        (Z-axis, bottom of box)
%   hl = half-length along X (cross-range half-width)
%   hw = half-width  along Y (range half-length)
%   hh = full height along Z
    Verts = [cx-hl, cy-hw, cz;       % 1 bottom-front-left
             cx+hl, cy-hw, cz;       % 2 bottom-front-right
             cx+hl, cy+hw, cz;       % 3 bottom-back-right
             cx-hl, cy+hw, cz;       % 4 bottom-back-left
             cx-hl, cy-hw, cz+hh;    % 5 top-front-left
             cx+hl, cy-hw, cz+hh;    % 6 top-front-right
             cx+hl, cy+hw, cz+hh;    % 7 top-back-right
             cx-hl, cy+hw, cz+hh];   % 8 top-back-left
    Faces = [1 2 3 4;   % bottom
             5 6 7 8;   % top
             1 2 6 5;   % front  (near radar)
             3 4 8 7;   % back   (far from radar)
             1 4 8 5;   % left
             2 3 7 6];  % right
end