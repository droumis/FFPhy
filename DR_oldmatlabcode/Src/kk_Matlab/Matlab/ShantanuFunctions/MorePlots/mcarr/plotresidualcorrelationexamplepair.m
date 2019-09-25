cell1_expected = nexpected(3,:);
cell1_resid = resid(3,:);
cell1_spikes = tmpspikes(3,:);

cell2_expected = nexpected(20,:);
cell2_resid = resid(20,:);
cell2_spikes = tmpspikes(20,:);

time = bins(1:end-1)+median(diff(bins))/2;


% Save cell 1 and cell 2 place fields
figure
plot(lf{pairsind(56,1)}.trajdata{2}(:,1),lf{pairsind(56,1)}.trajdata{j}(:,5),'b',lf{pairsind(56,1)}.trajdata{2}(:,1),lf{pairsind(56,2)}.trajdata{2}(:,5),'r')
xlabel('Linear distance (cm)')
ylabel('Firing rate (Hz)')

[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/ResidualCorrelation_placefields_example.pdf', m, d, y);
print('-dpdf', savestring)

% Save cell 1: expected and actual 100 sec
figure
hold on
plot(time,cell1_expected,'ko','MarkerFace','k')
plot(time,cell1_spikes,'bo','MarkerFace','b')
set(gca,'xlim',[1280 1380],'ylim',[0 6])
xlabel('Time (sec)')
ylabel('Number of spikes')

[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/ResidualCorrelation_cell1_100s.pdf', m, d, y);
print('-dpdf', savestring)

% Save cell 1: expected and actual 10 sec
set(gca,'xlim',[1310 1315],'ylim',[0 6])
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/ResidualCorrelation_cell1_10s.pdf', m, d, y);
print('-dpdf', savestring)

% Save cell 2: expected and actual 100 sec
figure
hold on
plot(time,cell2_expected,'ko','MarkerFace','k')
plot(time,cell2_spikes,'ro','MarkerFace','r')
set(gca,'xlim',[1280 1380],'ylim',[0 10])
xlabel('Time (sec)')
ylabel('Number of spikes')

[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/ResidualCorrelation_cell2_100s.pdf', m, d, y);
print('-dpdf', savestring)


% Save cell 2: expected and actual 10 sec
set(gca,'xlim',[1310 1315],'ylim',[0 6])
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/ResidualCorrelation_cell2_10s.pdf', m, d, y);
print('-dpdf', savestring)

% Save cell 1 & cell 2 residuals
timeind = lookup([1309:1316],time);
cell1_resid(timeind)
%         0         0    0.4794    0.7974    0.7817    0.0657

cell2_resid(timeind)
%     0.0129    0.0160    0.6981   -1.5298    1.5357    0.4314         0         0


% Plot Speed during time window
figure
plot(time,speed,'k')
set(gca,'xlim',[1280 1380])

[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/ResidualCorrelation_speed_100s.pdf', m, d, y);
print('-dpdf', savestring)

set(gca,'xlim',[1310 1315])
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/ResidualCorrelation_speed_10s.pdf', m, d, y);
print('-dpdf', savestring)

% This example pair correlation and speed correlation:
% correlation = 0.3558

% 1/4: 0.7899    1: 0.6674    4: 0.3851    16: 0.0626
