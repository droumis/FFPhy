
windowLength = 0.3;
winOffset = 0.1;

t = [0:30000*windowLength-1]/30000 - winOffset;

%Load continuous data segment for samson
cd /data13/monster/samson/sam1025/
d = dir('*.cont');
if isempty(d)
    error('No CONT files found.');
end

% generateTimesFromCont
% stimdio = createstimstruct; save stimdio1 stimdio

%%% NOTE- pin "32" appears to correspond to "A" = CA3 = pin 16
%%% NOTE- pin "31" appears to correspond to "B" = EC = pin 17

load times
ne = size(ranges,1)-1;
for i = 1:length(d)
	e = 9;
    filename = d(i).name;
    load stimdio1;
    if (e <= length(stimdio)) && (~isempty(stimdio{e}))
      ts1 = stimdio{e}.pulsetimes(:,1);
      [eeg1{i}, times, samplingRate] = cont_window_c(filename, ts/10000 - winOffset, windowLength);
    end
    load stimdio2;
    if (e <= length(stimdio)) && (~isempty(stimdio{e}))
      ts2 = stimdio{e}.pulsetimes(:,1);
      [eeg2{i}, times, samplingRate] = cont_window_c(filename, ts/10000 - winOffset, windowLength);
    end
end

EEG1 = cat(3,eeg1{:});
tmp = EEG1(:,:,1:2:end);
EEG1(:,:,1:2:end) = EEG1(:,:,2:2:end);
EEG1(:,:,2:2:end) = tmp;

EEG2 = cat(3,eeg2{:});
tmp = EEG2(:,:,1:2:end);
EEG2(:,:,1:2:end) = EEG2(:,:,2:2:end);
EEG2(:,:,2:2:end) = tmp;


figure
for i = 1:16
    plot(t,-1000*(i-1) - mean(EEG1(:,:,i),2),'k');
    hold on
end

figure
for i = 1:16
    plot(t,-1000*(i-1) - mean(EEG2(:,:,i),2),'k');
    hold on
end

load sampos02.mat
p = pos{2}{4}.data;

% PLOT SPEED RELATIONSHIP FOR BOTH STIMULATORS

windows = ts1/10000 - winOffset;
windows = [lookup(windows,p(:,1)) lookup(windows + windowLength,p(:,1))];

speed = nan(length(ts1),1);
for s = 1:length(ts1)
    speed(s) = mean(p(windows(s,1):windows(s,2),8));
end

% Plot valid stims for fast>4 and slow<1
slowind = speed>1/4 & speed<1;

fastind = speed>10;

figure
plot(t,squeeze(-mean(EEG1(:,slowind,15),2)),'r')
hold on
plot(t,squeeze(-mean(EEG1(:,fastind,15),2)),'g')
title('EC Stimulation')
set(gca,'ylim',[-4000 1500],'xlim',[-0.025 0.1])
legend('Slow','Fast')
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/%d_%d_%d_samson_EC_stim.pdf', m, d, y);
print('-dpdf', savestring)

figure
for i = 1:16
    plot(t,-1000*(i-1) - mean(EEG1(:,fastind,i),2),'g');
    hold on
    plot(t,-1000*(i-1) - mean(EEG1(:,slowind,i),2),'r');
end

%NOW FOR CA3
windows = ts2/10000 - winOffset;
windows = [lookup(windows,p(:,1)) lookup(windows + windowLength,p(:,1))];

speed = nan(length(ts2),1);
for s = 1:length(ts2)
    speed(s) = mean(p(windows(s,1):windows(s,2),8));
end

% Plot valid stims for fast>4 and slow<1
slowind = speed>1/4 & speed<1;

fastind =speed>10;

figure
plot(t,squeeze(mean(EEG2(:,slowind,5),2)),'r')
hold on
plot(t,squeeze(mean(EEG2(:,fastind,5),2)),'g')
title('CA3 Stimulation')
set(gca,'ylim',[-750 400],'xlim',[-0.025 0.1])
% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/%d_%d_%d_samson_CA3_stim.pdf', m, d, y);
print('-dpdf', savestring)


figure
for i = 1:16
    plot(t,-1000*(i-1) - mean(EEG2(:,fastind,i),2),'g');
    hold on
    plot(t,-1000*(i-1) - mean(EEG2(:,slowind,i),2),'r');
end
