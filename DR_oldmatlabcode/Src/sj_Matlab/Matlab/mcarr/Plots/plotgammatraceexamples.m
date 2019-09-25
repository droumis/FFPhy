% Plot examples of raw CA1 eeg during fast and slow gamma periods

%Load raw eeg file
load '/data13/mcarr/Cor/EEG/Coreeg01-2-23.mat'
e = eeg{1}{2}{23};
clear eeg

%Load fast and slow gamma events
load '/data13/mcarr/Cor/Corgammal01.mat'
slow = gammal{1}{2}{23};
clear gammal
load '/data13/mcarr/Cor/Corgammah01.mat'
fast = gammah{1}{2}{23};
clear gammah

% Define slow event and scale bar
event = -1*e.data(slow.startind(273):slow.endind(273));
time = slow.starttime(16):1/slow.samprate:slow.endtime(16);
time = time(100:500) - slow.starttime(16);
scalex = [time(1) time(151)];
scaley = [min(event) min(event)+50];
figure
hold on
plot(time,event(100:500),'k')
line(scalex,[scaley(1) scaley(1)],'Color','k')
line([scalex(1) scalex(1)], scaley,'Color','k')
set(gca,'xLim',[time(1)-0.05 time(end)+0.05],'ylim',[min(event)-50 max(event)+50])
legend('scale bar: 100ms by 50uV','Location','NorthEast')
box off
axis off

% Save figure
% print('-dpsc', '/data13/mcarr/VelocityPaper/SlowGamma/slowgamma_example_raweeg.ps')

% Define fast event and scale bar
event = -1*e.data(fast.startind(50):fast.endind(50));
time = fast.starttime(16):1/fast.samprate:fast.endtime(16);
time = time(100:500) - fast.starttime(16);
scalex = [time(1) time(151)];
scaley = [min(event)-50 min(event)];
figure
plot(time,event(100:500),'k')
hold on
line(scalex,[scaley(1) scaley(1)],'Color','k')
line([scalex(1) scalex(1)], scaley,'Color','k')
set(gca,'xLim',[time(1)-0.05 time(end)+0.05],'ylim',[min(event)-50 max(event)+50])
legend('scale bar: 100ms by 50uV','Location','NorthEast')
box off
axis off

% Save figure
%print('-dpsc', '/data13/mcarr/VelocityPaper/FastGamma/fastgamma_example_raweeg.ps')
