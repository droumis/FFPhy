
% 08/24/2012 - Plot eetgs from Probe 2

cd '/data25/sjadhav/HPExpt/Probe/PRb_direct/EEG';


load PRbeeg01-1-04.mat
eegstart = eeg{1}{1}{4}.starttime;
samprate = round(eeg{1}{1}{4}.samprate);
teeg = geteegtimes(eeg{1}{1}{4});

n=10000;

figure; hold on;
subplot(4,1,1); hold on;
plot(teeg(1:10000)-eegstart,eeg{1}{1}{4}.data(n+1:n+10000)); % Tet4Ch2 - Elec 8
set(gca,'YLim',[-450 250]); set(gca,'XLim',[-0.1 6]); axis off
load PRbeeg01-1-09.mat
subplot(4,1,2); hold on;
plot(teeg(1:10000)-eegstart,eeg{1}{1}{9}.data(n+1:n+10000)); % Tet9Ch3 - Elec 12
set(gca,'YLim',[-450 250]); set(gca,'XLim',[-0.1 6]); axis off
load PRbeeg01-1-08.mat
subplot(4,1,3); hold on;
plot(teeg(1:10000)-eegstart,eeg{1}{1}{8}.data(n+1:n+10000)); % Tet8Ch2 - Elec 13
set(gca,'YLim',[-450 250]); set(gca,'XLim',[-0.1 6]); axis off
load PRbeeg01-1-07.mat
subplot(4,1,4); hold on;
plot(teeg(1:10000)-eegstart,eeg{1}{1}{7}.data(n+1:n+10000)); % Tet7Ch0 - Elec 7
set(gca,'YLim',[-450 250]); set(gca,'XLim',[-0.1 6]);


