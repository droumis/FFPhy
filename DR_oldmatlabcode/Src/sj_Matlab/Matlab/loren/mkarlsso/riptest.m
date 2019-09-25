load /data/loren/mkarlsso/Con/concellinfo
load /data/loren/mkarlsso/Con/conripples04
load /data/loren/mkarlsso/Con/conspikes04

%s = spikes{8}{2}{29}{1}.data(:,1);
r1 = ripples{4}{2}{17};
r2 = ripples{4}{2}{23};

%hack
r1.maxthresh = ((r1.maxthresh * r1.std) - r1.baseline) / r1.std;
r2.maxthresh = ((r2.maxthresh * r2.std) - r2.baseline) / r2.std;

r1t = r1.midtime(r1.maxthresh > 3);
r2t = r2.midtime(r2.maxthresh > 3);

xc = spikexcorr(r1t, r2t, 0.050, 1);
plot(xc.time, xc.c1vsc2/xc.nspikes2, 'k', 'LineWidth', 2)
xlabel('Time (sec)');
ylabel('Probability')

figure
xc = spikexcorr(s, r2t, 0.001, .1);
plot(xc.time, xc.c1vsc2)


