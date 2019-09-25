function showRasters(rat, num1, num2)
%function showRasters(rat, num1, num2)
% compare raster plots of two cells
%  num= [d e t c]


if(num1(1)~= num2(1)) | (num1(2)~= num2(2))
    error('Cannot compare cells from different epochs.');
end

setRoot
load(sprintf('%s/%s/data2/spikedata%.2d.mat', root, rat, num1(1)));
sd{1}= spikedata{num1(1)}{num1(2)}{num1(3)}{num1(4)};
sd{2}= spikedata{num2(1)}{num2(2)}{num2(3)}{num2(4)};
load(sprintf('%s/%s/data2/behavdata%.2d.mat', root, rat, num1(1)));
bd= behavdata{num2(1)}{num2(2)};

auxPlotSpikePosTime(sd, bd)
%auxPlotSpikePhaseTime(sd, bd)
myprint([2 3], sprintf('raster-%s-%d-%d-%d-%d-vs-%d-%d', ...
    rat, num1,num2(3:4)), [], 0);

function auxPlotSpikePos(sd)
npart=20;
mint= 1e10; maxt= -1e10;
for i=1:length(sd)
    mint= min([sd{i}.time; mint]);
    maxt= max([sd{i}.time; maxt]);
end
tb= linspace(mint, maxt, npart+1);
y= [[0.15; 0.5], [0.5; 0.85]];
plotcol= {'k', 'r'};
for j=1:npart
    for i=1:length(sd)
        ind= find(tb(j)<sd{i}.time & sd{i}.time<tb(j+1));
        x= sd{i}.linpos(ind)';
        nsp= length(x);
        lh= line([x; x], j-1+y(:,i)*ones(1,nsp));
        set(lh, 'Color', plotcol{i});
        set(lh, 'LineWidth', 1);
    end
end
set(gca, 'YLim', [0 npart]);

function auxPlotSpikePosTime(sd, bd)
clf
hold on
bd.time= bd.time-bd.time(1);
sd{1}.time= sd{1}.time-bd.time(1);
sd{2}.time= sd{2}.time-bd.time(1);

%keyboard
%trajcol=[.8 .8 .8; 0 1 1; 1 1 0; 1 0 1;  0 1 0];
%for i=3:-1:-1
%    ind= find(bd.traj==i);
%    plot(bd.linpos(ind), bd.time(ind), '.', 'Color', trajcol(i+2,:));
%end
%plot(bd.linpos(1:100:end), bd.time(1:100:end), '-', ...
%    'LineWidth', 3, 'Color', .7*[1 1 1]);
plot(bd.linpos, bd.time, '.', 'MarkerSize', 2, 'Color', .7*[1 1 1]);
xlabel('position (cm)');
ylabel('time (s)');
set(gca, 'XLim', [min(bd.linpos)-2 max(bd.linpos)+2]);
set(gca, 'YLim', [min(bd.time)-15 max(bd.time)+15]);
%ind= find(bd.traj<0);
ind= find(bd.ripple);
plot(bd.linpos(ind), bd.time(ind), 'sy', 'MarkerSize', 4, 'MarkerFaceColor', 'y')
y= [[0.15; 0.5], [0.5; 0.85]];
%plot(sd{1}.linpos, sd{1}.time, 'k.');
%plot(sd{2}.linpos, sd{2}.time, 'ro', 'MarkerSize', 4);
spvalid= ismember(sd{2}.index, ind);
plot(sd{2}.linpos(spvalid), sd{2}.time(spvalid), 'ko', 'MarkerSize', 3);
spvalid= ismember(sd{1}.index, ind);
plot(sd{1}.linpos(spvalid), sd{1}.time(spvalid), 'r.', 'MarkerSize', 5);

function auxPlotSpikePhaseTime(sd, bd)
clf
hold on
%plot(bd.time(1:100:end), bd.phase(1:100:end), '-', ...
plot(bd.time, bd.phase, '-', ...
    'LineWidth', 3, 'Color', .7*[1 1 1]);
y= [[0.15; 0.5], [0.5; 0.85]];
plot(sd{1}.time, sd{1}.phase, 'k.');
plot(sd{2}.time, sd{2}.phase, 'ro', 'MarkerSize', 4);

%keyboard
