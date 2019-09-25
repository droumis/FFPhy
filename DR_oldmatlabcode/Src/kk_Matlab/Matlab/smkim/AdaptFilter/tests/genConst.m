function testConstGen()

behav_info.startTime=0;
behav_info.endTime=800;
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

gen.name='RescalingGen';

model.name='Const';
model.x=30;
spikes = genSpikes(data, model, gen);

data.spiketimes= spikes.spiketimes;
data.ispikes= spikes.ispikes;

save ConstTest data model
 
%test generated spikes
hist(data.spiketimes,25);
[N,X]= hist(data.spiketimes,25);
N=N/(X(2)-X(1));

figure(1); clf;
bar(N);
hold on;
plot(0:length(X)+1,model.x*ones(1,length(X)+2), 'r', 'linew', 3);
xlabel('time bins')
ylabel('spike rate (Hz)')
set(gca, 'xlim', [0, 26]);
legend('empirical','theoretical', 'Location', 'SouthOutside', 'orientation', 'horiz');
disp('Bars in histogram should be approx. around the red line.');
