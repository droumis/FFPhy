function [ output_args ] = behav_learningstate_est_reroute( directoryname,fileprefix,days, varargin )
lowercasethree = '';
% behav_learningstate_est_reroute(directoryname, fileprefix, days, varargin)
% estimates the probability that the rat in a single day will consistently
% pick the most preferred or 'comfortable' trajectory.  The preferred
% trajectory is calculated by find the unique trajectory that the rat takes
% in a day the most.  These traj are set as a success '1' while all other
% trajectories are counted as a failure '0'. This creates a series of
% bernoulli trials that can be processed by the state-spaced trial
% estimator (Smith et. al 2004) to estimate the probability and confidence
% interval that the preferred trajectory will be taken by the rat.
% Input:
%   directoryname - full directory path the ripdata save is
%   fileprefix - animal name/prefix
%   days - days to process
% Output:
%   Plots the estimate of probability the rat will take its preferred
%   trajectory and saves the plot

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

epoch_len = [];
day_len = [];
for day = days
    trialsegments_day = [];
    epochs = [];
    
    
    % Collect all trajectories across epochs in a day
    for i = 1:length(task_all{day})
        if ~isempty(task_all{day}{i}) && strcmp(task_all{day}{i}.type, 'reroute')
            trialsegments_day = [trialsegments_day linpos_all{day}{i}.trialsegments];
            epochs = [epochs i];
            epoch_len = [epoch_len size(linpos_all{day}{i}.trialsegments,2)];
        end
    end
    day_len = [day_len length(trialsegments_day)];
    
    % differentiate between the direction each trajectory is going
    % (forward/backwards)
    % ASSUMING END AND START SEGMENTS/WELLS ARE THE SAME FOR ALL EPOCHS IN
    % A DAY
    endpoint = linpos_all{day}{epochs(1)}.wellSegmentInfo.segmentIndex;
    forward_mask1 = cellfun(@(x)(x(1,4) == endpoint(1)),trialsegments_day(2,:));
    forward_mask2 = cellfun(@(x)(x(end,4) == endpoint(2)),trialsegments_day(2,:));
    back_mask1 = cellfun(@(x)(x(1,4) == endpoint(2)),trialsegments_day(2,:));
    back_mask2 = cellfun(@(x)(x(end,4) == endpoint(1)),trialsegments_day(2,:));
    
    % Error checking to make sure endpoint wells are the same, if not then
    % these should not be included (and these misassignments should be
    % corrected in the preprocessing code)
    if ~isequal(forward_mask1, forward_mask2) || ~isequal(back_mask1, back_mask2) || ~isequal((length(find(xor(forward_mask1,back_mask1)))), length(trialsegments_day(2,:)))
        disp(sprintf('ERROR: some trajectories do not go between segments %d and %d (day %d)!!',endpoint(1),endpoint(2),day));
        
    end
    
    % Calculate most frequent forward trajectory
    [forward_correct_traj F C] = mode([trialsegments_day{1,forward_mask1}]);

    if(length(C) > 1)
        disp('WARNING: More than one forward traj is maximum, this is not handled properly yet so it will most likely be misclassified');
    end
    
    % Calculate most frequency backwards trajectory
    [back_correct_traj F C] = mode([trialsegments_day{1,back_mask1}]);
    if(length(C) > 1)
        disp('WARNING: More than one backwards traj is maximum, this is not handled properly yet so it will most likely be misclassified');
    end


    % find all traj that match the most frequent
    correct_ind1 = find([trialsegments_day{1,:}]==forward_correct_traj);
    correct_ind2 = find([trialsegments_day{1,:}]==back_correct_traj);
    correct_ind = sort([correct_ind1 correct_ind2]);

    % assign '1's to most freq, '0's to all others
    binary_behav_correct{day} = zeros([1 length([trialsegments_day{1,:}])]);
    binary_behav_correct{day}(correct_ind) = repmat(1,1,length(correct_ind));
    
    
end
% estimate state
getestprobcorrect([binary_behav_correct{:}]',0.01,1);

%plotting and saving
subplot(2,1,1)
hold on;
set(gcf,'position',[0 0 800 250]);
set(gcf,'PaperPositionMode','auto');
set(gca,'FontSize',15);
for ii = 1:length(day_len)
    plot([sum(day_len(1:ii)) sum(day_len(1:ii))], [0 1], 'k');
    text(sum(day_len(1:ii)), .8, sprintf('Day %d\nn=%d',days(ii),day_len(ii)), ...
        'HorizontalAlignment','right','FontSize',10);
    
end

end

