function genStatic

% generate behavioral data and spikes based on static model
% for cell [1 2 1 1]

behav_info.startTime=0;		% [s]
%behav_info.endTime=200;	% [s]
behav_info.endTime=800;		% [s]
behav_info.dT = 1/1000;		% [s]
behav_info.minX = 0;
behav_info.maxX = 150; 		% [cm]
behav_info.velocity= 25;	% [cm/s]
behav_info.thetaFreq = 8;	% [Hz]
behavdata{1}{2}= genBehav(behav_info);

save('behavdata01', 'behavdata');

T= length(behavdata{1}{2}.time);
inc=ones(T,1);
genmodel.name='PosPhase_Isi';
genmodel.xp.name= 'Gauss';
genmodel.xp.alpha=  log(20)*inc;
genmodel.xp.mx=  50*inc;
genmodel.xp.my=  3.5*inc;
genmodel.xp.Sx=  15*inc;
genmodel.xp.Sy=  0.7*inc;
genmodel.xp.r=  tan(0.5*pi*-0.7)*inc;

genmodel.isi.name='Const';
genmodel.isi.a= 1;

gen.name='RescalingGen';

save('genmodel', 'genmodel', 'gen')

spikes = genSpikes(behavdata{1}{2}, genmodel, gen);

spikedata{1}{2}{1}{1}.time= spikes.spiketimes;
spikedata{1}{2}{1}{1}.index= spikes.ispikes;
spikedata{1}{2}{1}{1}.linpos= behavdata{1}{2}.linpos(spikes.ispikes);

save('spikedata01', 'spikedata');
