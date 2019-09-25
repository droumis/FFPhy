%load /data25/sjadhav/RippleInterruption/OthControlData/Con/concellinfo
%load /data25/sjadhav/RippleInterruption/OthControlData/Con/conripples04

% load /data25/sjadhav/RippleInterruption/REe_direct/REespikes05
% d=5;
% tet1=11; cell1=1;
% tet2=11; cell2=2;
 
load /data25/sjadhav/RippleInterruption/REf_direct/REfspikes03
d=3;
tet1=12; cell1=1;
tet2=1; cell2=3;


% Binsize, normalization etc like calcxcorrmeasures
binsize = 0.002;
winrun = 1;
winsleep = 0.4;
sw1 = 0.005; % Small gaussian window
sw2run = 0.25; % Big smoothing window
sw2sleep = 0.1; % Big smoothing window

% Run
ep1=2; ep2=4;
s1 = [spikes{d}{ep1}{tet1}{cell1}.data(:,1); spikes{d}{ep2}{tet1}{cell1}.data(:,1)];
s2 = [spikes{d}{ep1}{tet2}{cell2}.data(:,1); spikes{d}{ep2}{tet2}{cell2}.data(:,1)];
figure; hold on; redimscreen_halfvert; 
xc = spikexcorr(s1, s2, binsize, winrun);
subplot(2,1,1); hold on;
if ~isempty(xc.c1vsc2)
    %normcorr = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2);
    %g = gaussian(3, 18);
    %smoothcorr = smoothvect(xc.c1vsc2, g);
    [ec, normcorr, smoothcorr, basecorr] = excesscorr(xc.time, xc.c1vsc2,xc.nspikes1, xc.nspikes2, sw1, sw2run);
    plot(xc.time, normcorr, 'Color', [0.5 0.5 0.5]);
    plot(xc.time, smoothcorr,'k','LineWidth',3);
    plot(xc.time, basecorr,'k--','LineWidth',2);
    line([0 0], [0 max(normcorr)],'Color',[0.5 0.5 0.5],'LineWidth',2);
end
xlabel('Time (sec)');
ylabel('Probability')
title('Run Correlations (Ep 2 and 4)')
set(gca,'XLim',[-0.25 0.25 ]);


% Post- Sleep
ep1 = 3; ep2 = 5;
s1 = [spikes{d}{ep1}{tet1}{cell1}.data(:,1); spikes{d}{ep2}{tet1}{cell1}.data(:,1)];
s2 = [spikes{d}{ep1}{tet2}{cell2}.data(:,1); spikes{d}{ep2}{tet2}{cell2}.data(:,1)];
xcs = spikexcorr(s1, s2, binsize, winsleep);
subplot(2,1,2); hold on;
if ~isempty(xc.c1vsc2)
    %normcorr = xcs.c1vsc2 ./ sqrt(xcs.nspikes1 * xcs.nspikes2);
    %g = gaussian(3, 18);
    %smoothcorr = smoothvect(xcs.c1vsc2, g);
    [ec, normcorr, smoothcorr, basecorr] = excesscorr(xcs.time, xcs.c1vsc2,xcs.nspikes1, xcs.nspikes2, sw1, sw2sleep);
    plot(xcs.time, normcorr, 'Color', [0.5 0.5 0.5]);
    plot(xcs.time, smoothcorr,'k','LineWidth',3);
    plot(xcs.time, basecorr,'k--','LineWidth',2);
    line([0 0], [0 max(normcorr)],'Color',[0.5 0.5 0.5],'LineWidth',2);
end
xlabel('Time (sec)');
ylabel('Probability')
title('Sleep Correlations (Pre:Blue and Post-Beh)')
set(gca,'XLim',[-0.1 0.1]);

% Pre- Sleep
% ep1 = 1; ep2 = 1; % Double to count to mathc post-sleep
% s1 = [spikes{d}{ep1}{tet1}{cell1}.data(:,1); spikes{d}{ep2}{tet1}{cell1}.data(:,1)];
% s2 = [spikes{d}{ep1}{tet2}{cell2}.data(:,1); spikes{d}{ep2}{tet2}{cell2}.data(:,1)];
ep1=1;
s1 = [spikes{d}{ep1}{tet1}{cell1}.data(:,1)];
s2 = [spikes{d}{ep1}{tet2}{cell2}.data(:,1)];
xcs1 = spikexcorr(s1, s2, binsize, winsleep);
subplot(2,1,2); hold on;
if ~isempty(xc.c1vsc2)
    %normcorr = xcs1.c1vsc2 ./ sqrt(xcs1.nspikes1 * xcs1.nspikes2);
    %plot(xcs.time, normcorr, 'Color', [0.5 0.5 0.5]);
    %g = gaussian(3, 18);
    %smoothcorr = smoothvect(xcs.c1vsc2, g);
    [ec, normcorr, smoothcorr, basecorr] = excesscorr(xcs1.time, xcs1.c1vsc2,xcs1.nspikes1, xcs1.nspikes2, sw1, sw2sleep);
    plot(xcs1.time, smoothcorr,'b','LineWidth',3);
    plot(xcs1.time, basecorr,'b--','LineWidth',2);
end
%line([0 0], [0 max(normcorr)],'Color',[0.5 0.5 0.5],'LineWidth',2);
%xlabel('Time (sec)');
%ylabel('Probability')
%title('Sleep Correlations (Post-Beh)')
set(gca,'XLim',[-0.1 0.1]);



% r1 = ripples{4}{2}{17};
% r2 = ripples{4}{2}{23};
% %hack
% r1.maxthresh = ((r1.maxthresh * r1.std) - r1.baseline) / r1.std;
% r2.maxthresh = ((r2.maxthresh * r2.std) - r2.baseline) / r2.std;
% r1t = r1.midtime(r1.maxthresh > 1);
% r2t = r2.midtime(r2.maxthresh > 1);
% xc = spikexcorr(r1t, r2t, 0.050, 1);
% plot(xc.time, xc.c1vsc2/xc.nspikes2, 'k', 'LineWidth', 2)

% figure
% xc = spikexcorr(s1, r2t, 0.001, .1);
% plot(xc.time, xc.c1vsc2)


