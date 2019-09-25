function stats= runStat

load(['dynamic_behav']);
load(['dynamic_data']);
load(['dynamic1d_result']);


%sopts{1}.name= 'MutualInfo';
%sopts{1}.nbinx= 20;
%sopts{1}.nbiny= 20;

%opt.name= 'LinCorr';
%opt.nbinx= 100;
%opt.nbiny= 100;
%sopts{2}= opt;


sopts{3}.name= 'Moments';
sopts{3}.sizex= 1;


lim.a{1}.traj=0;
lim.a{1}.linpos=[1; 149];
lim.a{1}.phase=[0;2*pi];
mintime=behav_info.startTime;
maxtime= behav_info.endTime-1;
dt=(maxtime-mintime)/10;
time= [mintime:dt:maxtime]';
lim.x{1}.time= time;

stats= calcStatistics(data, result, sopts, lim);

save(['dynamic1d_stats'], 'stats');

figure
ystr= {'area'; 'mean'; 'std'; 'skewness'};
for i=1:4
    subplot(2,2,i)
    plot(time, stats.Moments{1}(:,i))
    xlabel('time (s)');
    ylabel(ystr{i});
    hold on
end

subplot(2,2,2)
plot(data.time, genmodel.spatial.mx, 'r');
subplot(2,2,3)
plot(data.time, genmodel.spatial.Sx, 'r');
legend({'dyn. est.'; 'real'}, 'location', 'best');
subplot(2,2,4)
plot(data.time, zeros(1,length(data.time)), 'r');

clear functions
