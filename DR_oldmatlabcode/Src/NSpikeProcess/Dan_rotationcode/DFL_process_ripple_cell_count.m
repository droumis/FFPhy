function process_ripple_cell_count(directoryname,fileprefix,days, winsize, varargin)
% process_ripple_cell_count(directoryname,fileprefix,days, winsize, varargin)
% function that will process the data saved by the ripdata (unique cells
% firing per ripple) and plot it smoothed with a moving average filter
% Input:
%   directoryname - full directory path the ripdata save is
%   fileprefix - animal name/prefix
%   days - days to process
%   winsize - window size of moving average
% Output:
%   Plots the moving average ripple unique cell count and saves it


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
    eval(['load ',directoryname,fileprefix,'_',dsz,num2str(day),'task', '.mat']);
    eval(['task = ',lowercasethree,'task;'])
    eval(['load ',directoryname,fileprefix,'ripple_cell_count',dsz,num2str(day), '.mat']);
    
    linpos_all{day} = linpos{day};
    task_all{day} = task{day};
    ripcell_all{day} = ripple_cell_count{day};
end

% setup variables
epoch_len = [];
day_len = [];
day_totalcell = [];
ripcell_mean = [];
ripcell_std = [];

% loop through all days
for day = days
    ripcell_day = [];
    epochs = [];
    for i = 1:length(task_all{day})
        if ~isempty(task_all{day}{i}) && strcmp(task_all{day}{i}.type, 'reroute')
            ripcell_day = [ripcell_day ripcell_all{day}{i}];
            epochs = [epochs i];
            epoch_len = [epoch_len length(ripcell_all{day}{i})];
        end
    end
    day_len = [day_len length(ripcell_day)];
    
    
    % calculate total cells per day and moving averages
    day_totalcell = [day_totalcell sum(ripcell_day)];
    ripcell_mean = [ripcell_mean movingstat(ripcell_day',winsize,@mean)'];
    ripcell_std = [ripcell_std movingstat(ripcell_day',winsize,@std)'];
end


% Plotting

figure(1);
set(gcf,'position',[0 0 800 250]);
set(gcf,'PaperPositionMode','auto');
set(gca,'FontSize',15);
hold off;
plot(ripcell_mean,'LineWidth',3);
hold on;
plot(ripcell_mean+ripcell_std,'r:','LineWidth',2);
plot(ripcell_mean-ripcell_std,'r:','LineWidth',2);

for ii = 1:length(day_len)
    plot([sum(day_len(1:ii)) sum(day_len(1:ii))]-(winsize-1)*(ii), [0 6], 'k');
    text(sum(day_len(1:ii))-(winsize)*ii, 5, sprintf('Day %d\nn=%d\n#cells=%d',ii,day_len(ii),day_totalcell(ii)), ...
        'HorizontalAlignment','right','FontSize',10);
    
end
ylim([0 6]);
title(sprintf('Moving Average (win=%d) Unique Cells Per Ripple (D%d-%d)',winsize,days(1),days(end)));
xlabel('Ripples');
ylabel('Average Cells Per Ripple');

[s,mess,messid] = mkdir(sprintf('%sPlot/ripple/',directoryname));
print(sprintf('%sPlot/ripple/%s_ripcellavg_d%d-%d_w%d',directoryname,fileprefix,days(1),days(2),winsize),'-depsc');
