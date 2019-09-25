function genDynamic

behav_info.startTime=0;		% [s]
%behav_info.endTime=200;	% [s]
behav_info.endTime=800;		% [s]
behav_info.dT = 1/1000;		% [s]
behav_info.minX = 0;
behav_info.maxX = 150; 		% [cm]
behav_info.velocity= 25;	% [cm/s]
behav_info.thetaFreq = 8;	% [Hz]
behavdata{1}{2}= genBehav(behav_info);

%behav_info.minY = 0;
%behav_info.maxY = 150;
%behav_info.velX = 25;
%behav_info.velY = 17;
%data= genBehavXY(behav_info);

%showBehav(data,1)

save('behavdata01', 'behavdata');

T= length(behavdata{1}{2}.time);

genmodel.name='PosPhase_Isi';
incr=[1:T]*(behav_info.dT/60); % [min] , figure ca. 10 min total
genmodel.xp.name= 'Gauss';
genmodel.xp.alpha=  log(10 + 1.125 *incr);
genmodel.xp.mx=  25 + 7.5*incr;
genmodel.xp.my=  2.5 + 0.1*incr;
genmodel.xp.Sx=  12 + 0.45* incr;
genmodel.xp.Sy=  0.5 + 0.01*incr;
genmodel.xp.r=  tan(0.5*pi*-0.045*incr);
% genmodel.xp.alpha=  log(10 + 1.125 *incr);
% genmodel.xp.mx=  25 + 7.5*incr;
% genmodel.xp.my=  120 - 6*incr;
% genmodel.xp.Sx=  12 + 0.45* incr;
% genmodel.xp.Sy=  11 + 0.5*incr;
% genmodel.xp.r=  tan(0.5*pi*-0.045*incr);
%genmodel.xp.r= zeros(1,T);

genmodel.isi.name='Const';
genmodel.isi.a= 1;

gen.name='RescalingGen';

save('genmodel', 'genmodel', 'gen')

spikes = genSpikes(behavdata{1}{2}, genmodel, gen);

spikedata{1}{2}{1}{1}.time= spikes.spiketimes;
spikedata{1}{2}{1}{1}.index= spikes.ispikes;
spikedata{1}{2}{1}{1}.linpos= behavdata{1}{2}.linpos(spikes.ispikes);

save('spikedata01', 'spikedata');
