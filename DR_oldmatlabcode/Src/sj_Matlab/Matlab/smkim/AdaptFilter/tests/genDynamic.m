function genDynamic

behav_info.startTime=0;
%behav_info.endTime=300;
behav_info.endTime=600;
behav_info.dT = 2/1000;
behav_info.minX = 0;
behav_info.maxX = 150;
behav_info.velocity= 25;
behav_info.thetaFreq = 8;

behav_info.minY = 0;
behav_info.maxY = 150;
behav_info.velX = 25;
behav_info.velY = 17;


data= genBehav(behav_info);
%data= genBehavXY(behav_info);
save dynamic_behav behav_info data

%showBehav(data,1)

T= length(data.time);

genmodel.name='PosPhase_Isi';
genmodel.max_isi= 800;
incr=[1:T]/T;
genmodel.spatial.name= 'Gauss';

%genmodel.spatial.alpha=  log(15 + 10 *incr);
%genmodel.spatial.mx=  25 + 30*incr;
%genmodel.spatial.my=  2.5 + 0.5*incr;
%genmodel.spatial.Sx=  12 + 4* incr;
%genmodel.spatial.Sy=  0.5 + 0.1*incr;
%genmodel.spatial.r=  tan(0.5*pi*-0.07*incr);

genmodel.spatial.alpha=  log(20 + 0 *incr);
genmodel.spatial.mx=  25 + 0*incr;
genmodel.spatial.my=  2.5 + 0.00*incr;
genmodel.spatial.Sx=  10 + 10.* incr;
genmodel.spatial.Sy=  0.5 + 0.0*incr;
genmodel.spatial.r=  tan(0.5*pi*(-0.7*ones(1,T)));
%genmodel.spatial.r=  tan(0.5*pi*(-0.7*incr));
%genmodel.spatial.r(1:T/2)=  tan(0.5*pi*(zeros(1,T/2)));
%genmodel.spatial.r(T/2+1:T)=  tan(0.5*pi*(-0.7*ones(1,T/2)));

genmodel.isi.name='Const';
genmodel.isi.x= 1;

gen.name='RescalingGen';

spikes = genSpikes(data, genmodel, gen);

data.spiketimes= spikes.spiketimes;
data.spikeindex= spikes.ispikes;

% [psth,x,y]= showHist(data, 1);
% [mx,my,Vx,Vy,R,I]=Hist2Stat(psth,x,y)

save('dynamic_data', 'data', 'genmodel', 'T')
