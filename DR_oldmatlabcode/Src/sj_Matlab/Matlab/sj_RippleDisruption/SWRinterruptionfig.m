%s1a, day 5

load /data/loren/shantanu/S1a/EEG/s1aeeg10-3-04;
load /data/loren/shantanu/S1a/EEG/s1aripple10-3-04;
e3 = eeg{10}{3}{4};
r3 = ripple{10}{3}{4};

load /data/loren/shantanu/S1a/EEG/s1aeeg10-4-04;
load /data/loren/shantanu/S1a/EEG/s1aripple10-4-04;
e4 = eeg{10}{4}{4};
r4 = ripple{10}{4}{4};

load /data/loren/shantanu/S1a/EEG/s1aeeg10-5-04;
load /data/loren/shantanu/S1a/EEG/s1aripple10-5-04;
e5 = eeg{10}{5}{4}
r5 = ripple{10}{5}{4};

%e3
figure
e3t = [4677.1 4677.6];
e3t = [6638.05 6638.45]
subplot(3,1,1)
t = geteegtimes(e3);
tind = find((t >= e3t(1)) & (t < e3t(2)));
plot(t(tind), e3.data(tind), 'k');
hold on
plot(t(tind), r3.data(tind,1), 'Color', [.8 0 0]);
set(gca, 'YLim', [-600 600])
set(gca, 'XLim', e3t)
h = title('Before disruption');
axis off

%e2
subplot(3,1,2)
e4t = [8405 8407.35];
t = geteegtimes(e4);
tind = find((t >= e4t(1)) & (t < e4t(2)));
plot(t(tind), e4.data(tind), 'k');
hold on
plot(t(tind), r4.data(tind,1), 'Color', [.8 0 0]);
set(gca, 'YLim', [-600 600])
set(gca, 'XLim', e4t)
h = title('During disruption');
axis off

%e5
subplot(3,1,3)
e5t = [9174.96 9175.36];
t = geteegtimes(e5);
tind = find((t >= e5t(1)) & (t < e5t(2)));
plot(t(tind), e5.data(tind), 'k');
hold on
plot(t(tind), r5.data(tind,1), 'Color', [.8 0 0]);
set(gca, 'XLim', e5t)
set(gca, 'YLim', [-600 600])
plot([e5t(2) - .1 e5t(2)], [-400 -400], 'k', 'LineWidth', 2);
plot([e5t(2) e5t(2)], [-200 -400], 'k', 'LineWidth', 2);
h = title('After disruption');
axis off



% size and duration figure
load /data/loren/shantanu/S1a/s1aripples10
r3 = ripples{10}{3}{4};
valid = r3.maxthresh > 3;
r3.maxthresh = r3.maxthresh(valid);
r3.duration = r3.endtime(valid) - r3.starttime(valid);;
r4 = ripples{10}{4}{4};
valid = r4.maxthresh > 3;
r4.maxthresh = r4.maxthresh(valid);
r4.duration = r4.endtime(valid) - r4.starttime(valid);;
r5 = ripples{10}{5}{4};
valid = r5.maxthresh > 3;
r5.maxthresh = r5.maxthresh(valid);
r5.duration = r5.endtime(valid) - r5.starttime(valid);;

figure
orient tall
subplot(4,2,1)
%[a b] = flhist(r3.maxthresh, 0:.5:20);
%plot(b,a./sum(a), 'LineWidth', 2);
[c3, x] = ecdf(r3.maxthresh);
plot(x, c3, 'b', 'LineWidth', 2);
hold on
[c4, x] = ecdf(r4.maxthresh);
plot(x, c4, 'r', 'LineWidth', 2);
[c5, x] = ecdf(r5.maxthresh);
plot(x, c5, 'k', 'LineWidth', 2);
%[a b] = flhist(r4.maxthresh, 0:.5:20);
%plot(b,a./sum(a), 'r', 'LineWidth', 2);
%[a b] = flhist(r5.maxthresh, 0:.5:20);
%plot(b,a./sum(a), 'k', 'LineWidth', 2);
set(gca, 'XLim', [0 20]);
xlabel('Ripple size (std)');
ylabel('Cumulative proportion');

subplot(4,2,2)
%[a b] = flhist((r3.endtime - r3.starttime), 0:.01:.35);
%plot(b,a./sum(a), 'LineWidth', 2);
[d3, x] = ecdf(r3.duration);
plot(x, d3, 'b', 'LineWidth', 2);
hold on
%[a b] = flhist((r4.endtime - r4.starttime), 0:.01:.35);
%plot(b,a./sum(a), 'r', 'LineWidth', 2);
%[a b] = flhist((r5.endtime - r5.starttime), 0:.01:.35);
%plot(b,a./sum(a), 'k', 'LineWidth', 2);
[d4, x] = ecdf(r4.duration);
plot(x, d4, 'r', 'LineWidth', 2);
[d5, x] = ecdf(r5.duration);
plot(x, d5, 'k', 'LineWidth', 2);

xlabel('Ripple duration (sec)');
set(gca, 'XLim', [0 0.25]);

legend({'Before', 'During', 'After'});

figure

e4t = [8407 8407.15];
t = geteegtimes(e4);
tind = find((t >= e4t(1)) & (t < e4t(2)));
plot(t(tind), -e4.data(tind), 'k');
hold on
%set(gca, 'YLim', [-600 600])
set(gca, 'XLim', e4t)
h = title('During disruption');

