% set these variables before evaluating and run this in the animal's processed directory
% (example: '/data13/mcarr/thr'). Set epoch equal to the last linearized
% run.

epoch=2;
index = [2 epoch 3 1];
fileprefix = 'thr';
binsize = 2.5;
std = 2;
   
% load spikes, pos, linpos, task for day
eval(['load ',fileprefix,'spikes','0',num2str(index(1)), '.mat']);
eval(['load ',fileprefix,'pos','0',num2str(index(1)), '.mat']);
eval(['load ',fileprefix,'linpos','0',num2str(index(1)), '.mat']);
eval(['load ',fileprefix,'task','0',num2str(index(1)), '.mat']);

%run getbehave state, calclinfields, and twodoccupancy maps
[state lindist] = getbehavestate(linpos, index(1), 2, 1);
trajdata = calclinfields(spikes, state, lindist, linpos, [index(1) 2 index(3) index(4)]);
[output] = twodoccupancy(spikes,state, linpos, pos, [index(1) 2 index(3) index(4)]);

% plot the trajdata linearized firing rates and the twodoccupancy pictures
figure;
subplot(2,4,1); plot(trajdata{1}(:,5))
subplot(2,4,2); imagesc(output(1).smoothedspikerate)
subplot(2,4,3); plot(trajdata{2}(:,5))
subplot(2,4,4); imagesc(output(2).smoothedspikerate)
subplot(2,4,5); plot(trajdata{3}(:,5))
subplot(2,4,6); imagesc(output(3).smoothedspikerate)
subplot(2,4,7); plot(trajdata{4}(:,5))
subplot(2,4,8); imagesc(output(4).smoothedspikerate)


if max(epoch==4)
    figure;
    subplot(1,2,1); plot(pos{index(1)}{2}.data(:,2), pos{index(1)}{2}.data(:,3), 'x')
        hold on; plot(spikes{index(1)}{2}{index(3)}{index(4)}.data(:,2), spikes{index(1)}{2}{index(3)}{index(4)}.data(:,3), 'rx')
    subplot(1,2,2); plot(pos{index(1)}{4}.data(:,2), pos{index(1)}{4}.data(:,3), 'x')
        hold on; plot(spikes{index(1)}{4}{index(3)}{index(4)}.data(:,2), spikes{index(1)}{4}{index(3)}{index(4)}.data(:,3), 'rx')

    [state2 lindist2] = getbehavestate(linpos, index(1), 4, 1);
    trajdata2 = calclinfields(spikes, state2, lindist2, linpos, [index(1) 4 index(3) index(4)]);
    [output2] = twodoccupancy(spikes,state2, linpos, pos, [index(1) 4 index(3) index(4)]);
    [output3] = openfieldoccupancy(spikes, pos, [index(1) 6 index(3) index(4)], binsize, std);
    
    figure;
    subplot(1,2,1); imagesc(output3.smoothedspikerate)
    subplot(1,2,2); plot(pos{index(1)}{6}.data(:,2), pos{index(1)}{6}.data(:,3), 'x')
        hold on; plot(spikes{index(1)}{6}{index(3)}{index(4)}.data(:,2), spikes{index(1)}{6}{index(3)}{index(4)}.data(:,3), 'rx')
    
    figure;
        subplot(2,4,1); plot(trajdata2{1}(:,5))
        subplot(2,4,2); imagesc(output2(1).smoothedspikerate)
        subplot(2,4,3); plot(trajdata2{2}(:,5))
        subplot(2,4,4); imagesc(output2(2).smoothedspikerate) 
        subplot(2,4,5); plot(trajdata2{3}(:,5))
        subplot(2,4,6); imagesc(output2(3).smoothedspikerate)
        subplot(2,4,7); plot(trajdata2{4}(:,5))
        subplot(2,4,8); imagesc(output2(4).smoothedspikerate)
else
    % plot the open field occupancy normalized spike rate
    figure;
    plot(pos{index(1)}{2}.data(:,2), pos{index(1)}{2}.data(:,3), 'x')
        hold on; plot(spikes{index(1)}{2}{index(3)}{index(4)}.data(:,2), spikes{index(1)}{2}{index(3)}{index(4)}.data(:,3), 'rx')

    [output2] = openfieldoccupancy(spikes, pos, [index(1) 4 index(3) index(4)], binsize, std);
    figure;
    subplot(1,2,1)
    imagesc(output2.smoothedspikerate)
    subplot(1,2,2)
    plot(pos{index(1)}{4}.data(:,2),pos{index(1)}{4}.data(:,3),'x')
    hold on
    plot(spikes{index(1)}{4}{index(3)}{index(4)}.data(:,2), spikes{index(1)}{4}{index(3)}{index(4)}.data(:,3), 'rx')
end
