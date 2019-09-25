function [avecorrect] = getavecorrect(data, datasetnum, condition, ntrials)
% [avecorrect] = getavecorrect(data, datasetnum, condition)
%     goes through each condition for each dataset and get the running
%     average of the percentage correct using an ntrials moving average.
%     The first column of avecorrect is the trial number and the second column
%     is the percent correct

dataset = data{datasetnum};
correct = (dataset.cobj.response == 0);
% figure out what the number of this condition is */
[condnum tmp] = find(dataset.conds == condition)

if (isempty(condnum))
    % that specified condition is not present
    error('Incorrect condition number');
    avecorrect = [];
    return;
end

trials = dataset.condtrial{condnum};

ctrials = correct(trials);

% the "filter" here is 1/ntrials for all elements, producing a five trial
% moving average. We could also use a gaussian or any other function as long as
% it was normalized to have a total area (eg. sum) of one
filter = ones(1,ntrials) / ntrials;

% the next line looks a bit complicated, but it is just a way to get the right
% starting trial number for an arbitrary sized average
avecorrect(:,1) = [1:length(ctrials)]';

% to get the average correct, we convolve the filter with the data and remove
% the end points. 
avecorrect(:,2) = smoothvect(ctrials, filter);
