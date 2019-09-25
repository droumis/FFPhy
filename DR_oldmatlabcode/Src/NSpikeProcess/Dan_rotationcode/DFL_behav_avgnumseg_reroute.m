function [ output_args ] = behav_avgnumseg_reroute( directoryname,fileprefix,days, winsize, varargin )
% behav_avgnumseg_reroute(directoryname,fileprefix,days, winsize, varargin)
% function that will process the trajectory/trials from the linearized data
% (linpos), calculate the number of segments each trajectory took,
% and plot it smoothed with a moving average filter
% Input:
%   directoryname - full directory path the ripdata save is
%   fileprefix - animal name/prefix
%   days - days to process
%   winsize - window size of moving average
% Output:
%   Plots the moving average number of segments per trial and saves it

lowercasethree = '';

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'lowercasethree'
            lowercasethree = varargin{option+1};
            
    end
end

days = days(:)';


%Load all days to process
for day = days
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    eval(['load ',directoryname,fileprefix,'linpos',dsz,num2str(day), '.mat']);
    eval(['linpos = ',lowercasethree,'linpos;'])
    eval(['load ',directoryname,fileprefix,'task',dsz,num2str(day), '.mat']);
    eval(['task = ',lowercasethree,'task;'])
    
    linpos_all{day} = linpos{day};
    task_all{day} = task{day};
    
end

% setup variables
epoch_len = [];
day_len = [];
trialsegments_all = [];
day_totalseg = [];
numseg_mean = [];
numseg_std = [];
day_totaltime = [];
totaltime_mean = [];
totaltime_std = [];

% loop through days
for day = days
    
    % collecting trial segments across epochs and store which trial each day and
    % epoch ends so it can be used to plot later
    trialsegments_day = [];
    epochs = [];
    for i = 1:length(task_all{day})
        if ~isempty(task_all{day}{i}) && strcmp(task_all{day}{i}.type, 'reroute')
            trialsegments_day = [trialsegments_day linpos_all{day}{i}.trialsegments];
            epochs = [epochs i];
            epoch_len = [epoch_len size(linpos_all{day}{i}.trialsegments,2)];
        end
    end
    day_len = [day_len length(trialsegments_day)];
    
    % Extract length of each trial (# seg) and total time spent in each
    % trial
    num_seg = [cellfun(@(x) (length(x)), trialsegments_day(2,:))];
    trial_time = [cellfun(@(x) (x(end,1)-x(1,1)), trialsegments_day(2,:))];
    
    % Collect the lengths of each trial across all epochs in the day
    day_totalseg = [day_totalseg sum(num_seg)];
    numseg_mean = [numseg_mean movingstat(num_seg',winsize,@mean)'];
    numseg_std = [numseg_std movingstat(num_seg',winsize,@std)'];
    
    % Collect the time of each trial across all epochs in the day
    day_totaltime = [day_totaltime sum(trial_time)];
    totaltime_mean = [totaltime_mean movingstat(trial_time',winsize,@mean)'];
    totaltime_std = [totaltime_std movingstat(trial_time',winsize,@std)'];
    
    %trialsegments_all = [trialsegments_all trialsegments_day];
    
end

% Plotting number of segments per trial
figure(1);
set(gcf,'position',[0 0 800 250]);
set(gcf,'PaperPositionMode','auto');
set(gca,'FontSize',15);
hold off;
plot(numseg_mean,'LineWidth',3);
hold on;
plot(numseg_mean+numseg_std,'r:','LineWidth',2);
plot(numseg_mean-numseg_std,'r:','LineWidth',2);

% Annotating each day
for ii = 1:length(day_len)
    plot([sum(day_len(1:ii)) sum(day_len(1:ii))]-(winsize-1)*(ii), [0 15], 'k');
    text(sum(day_len(1:ii))-(winsize)*ii, 13, sprintf('Day %d\nn=%d\n#seg=%d',days(ii),day_len(ii),day_totalseg(ii)), ...
        'HorizontalAlignment','right','FontSize',10);
    
end
ylim([0 15]);
title(sprintf('Moving Average (win=%d) Num Seg per Trial (D%d-%d)',winsize,days(1),days(end)));
xlabel('Trials');
ylabel('Moving Average Num Seg');

% Saving plot
[s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
print(sprintf('%sPlot/behav/%s_numsegavg_d%d-%d_w%d',directoryname,fileprefix,days(1),days(end),winsize),'-depsc');

% Plotting time per trial
figure(2);
set(gcf,'position',[0 0 800 250]);
set(gcf,'PaperPositionMode','auto');
set(gca,'FontSize',15);
hold off;
plot(totaltime_mean,'LineWidth',3);
hold on;
plot(totaltime_mean+totaltime_std,'r:','LineWidth',2);
plot(totaltime_mean-totaltime_std,'r:','LineWidth',2);

% annotating
for ii = 1:length(day_len)
    plot([sum(day_len(1:ii)) sum(day_len(1:ii))]-(winsize-1)*(ii), [0 200], 'k');
    text(sum(day_len(1:ii))-(winsize)*ii, 85, sprintf('Day %d\nn=%d\ntime=%d',days(ii),day_len(ii),floor(day_totaltime(ii))), ...
        'HorizontalAlignment','right','FontSize',10);
    
end
ylim([0 100]);
title(sprintf('Moving Average (win=%d) Time per Trial (D%d-%d)',winsize,days(1),days(end)));
xlabel('Trials');
ylabel('Moving Average Trial Time');

% saving
[s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
print(sprintf('%sPlot/behav/%s_trialtimeavg_d%d-%d_w%d',directoryname,fileprefix,days(1),days(end),winsize),'-depsc');