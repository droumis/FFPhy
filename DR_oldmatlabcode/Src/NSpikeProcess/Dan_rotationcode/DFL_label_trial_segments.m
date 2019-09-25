
function [labeled_trial_segments segmenttable trajwells] = label_trial_segments(trial_segments)
% [labeled_trial_segments segmenttable trajwells] = label_trial_segments(trial_segments)
% Takes all trial segments and labels each unique trajectory with a unique
% id.  Also generates the segmenttable and corresponding trajwells that
% stores all unique trajectories.
% Input:
%   trial_segments - must be a 1xn cell array that contains all n trials that you
%   want to assign unique ids to.
% Output:
%   labeled_trial_segments - a 2xn cell array with the first row is the
%   unique trajectory ID of a trial and the second is the same as
%   trial_segments
%   segmenttable - table of trajectories and their segments
%   trajwells - the end wells of each trajectory

max_trial_size = max(cellfun(@(x) size(x,1),trial_segments)); % count size of each cell (each trial) in cell array
% convert trajectory cell array into a matrix, easier to work with matrix
% padded with 0s than cell array
trials = cell2mat(cellfun(@(x,y) padarray(x(:,4),[y-size(x,1),0],'post'),trial_segments, ...
    num2cell(repmat(max_trial_size,1,length(trial_segments))),'UniformOutput',false))'; % turns it into a cell array

% find unique trials and create IDs
[unique_trials source_ind unique_ids] = unique(trials,'rows');

labeled_trial_segments = cell(2,length(trial_segments));
labeled_trial_segments(1,:) = num2cell(unique_ids');
labeled_trial_segments(2,:) = trial_segments;

% Calculate segmenttable
segmenttable = [];


trajwells = [];
for ii = 1:length(source_ind)
    seg_full = trial_segments{source_ind(ii)};
    seg = seg_full(:,4);
    segmenttable = [segmenttable; repmat(ii,length(seg),1), (1:length(seg))', seg];
    trajwells = [trajwells; seg_full(1,2), seg_full(end,3)];
end

