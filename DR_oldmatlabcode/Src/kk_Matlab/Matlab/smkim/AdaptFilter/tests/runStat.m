function stats= runStat

load(['dynamic_behav']);
load(['dynamic_data']);
load(['dynamic_result']);


%sopts{1}.name= 'MutualInfo';
%sopts{1}.nbinx= 20;
%sopts{1}.nbiny= 20;

%opt.name= 'LinCorr';
%opt.nbinx= 100;
%opt.nbiny= 100;
%sopts{2}= opt;


sopts{3}.name= 'Moments2d';
sopts{3}.sizex= 1;
sopts{3}.sizey= .1;


lim.a{1}.traj=0;
lim.a{1}.linpos=[1; 149];
lim.a{1}.phase=[0;2*pi];
mintime=behav_info.startTime;
maxtime= behav_info.endTime-1;
dt=(maxtime-mintime)/10;
lim.x{1}.time= [mintime:dt:maxtime];


stats= calcStatistics(data, result, sopts, lim);

save(['dynamic_stats'], 'stats');
