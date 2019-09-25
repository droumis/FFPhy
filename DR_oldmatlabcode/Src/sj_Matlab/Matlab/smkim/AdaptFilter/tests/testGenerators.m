%function testGenerators
%function testGenerators
%
% Test the algorithm for simulating spike trains

% simulate animal's behavior
behav_info.startTime=0;
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
T= length(data.time);

% all tests with generator based on time-rescaling theorem
gen.name='RescalingGen';    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate spikes with a constant firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model.name='Const';
model.x=30;
spikes = genSpikes(data, model, gen);

data.spiketimes= spikes.spiketimes;
data.ispikes= spikes.ispikes;

% check generated spikes
hist(data.spiketimes,25);
[N,X]= hist(data.spiketimes,25);
N=N/(X(2)-X(1));

figure
bar(N);
hold on;
plot(0:length(X)+1,model.x*ones(1,length(X)+2), 'r', 'linew', 3);
xlabel('time bins')
ylabel('spike rate (Hz)')
set(gca, 'xlim', [0, 26]);
legend('empirical','theoretical', 'Location', 'SouthOutside', 'orientation', 'horiz');
title('constant firing rate')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate spikes with a two-factor model (spatial x isi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear genmodel n

genmodel.name='PosPhase_Isi';
%genmodel.max_isi= 0.08; % (sec)
genmodel.max_isi= behav_info.endTime;
incr=[1:T]/T;
genmodel.spatial.name= 'Gauss';

genmodel.spatial.alpha=  log(20 + 0 *incr);
genmodel.spatial.mx=  25 + 0*incr;
genmodel.spatial.my=  2.5 + 0.00*incr;
genmodel.spatial.Sx=  10 + 10.* incr;
genmodel.spatial.Sy=  0.5 + 0.0*incr;
genmodel.spatial.r=  tan(0.5*pi*(-0.7*ones(1,T)));

genmodel.isi.name='Const';

x= [1 2 4 8 16];
for i=1:length(x)
    genmodel.isi.x= x(i);
    spikes = genSpikes(data, genmodel, gen);
    n(i)= length(spikes.ispikes);
end

scale= sum(n)/sum(x);
figure
bar([n', scale*x'])
legend('empirical','theoretical', 'Location', 'northwest', 'orientation', 'horiz');

clear functions
