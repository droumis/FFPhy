function [logistcorrect] = getlogistcorrect(data, datasetnum, condition)
% [logistcorrect] = getlogistcorrect(data, datasetnum, condition)
%     goes through each condition for each dataset and fits a logistic model to the data. 
%     The first column of logistcorrect is the trial number and the second column
%     is the probability of a correct response

dataset = data{datasetnum};
correct = (dataset.cobj.response == 0);
% figure out what the number of this condition is */
[condnum tmp] = find(dataset.conds == condition)

if (isempty(condnum))
    % that specified condition is not present
    error('Incorrect condition number');
    logistcorrect = [];
    return;
end

trials = dataset.condtrial{condnum};

ctrials = correct(trials);

logistcorrect(:,1) = [1:length(ctrials)]';

% fit the logistic model  (1-p)/p = theta + x * beta
[beta, theta, dev, dl, d2l, p] = logist(ctrials, [1:length(ctrials)]', 0);


% this model is easier to work with if it fit p/1-p, so we flip the sign of theta
theta = theta * -1;

% the following produces the right fit, but I'm not quite sure why....
% more on this later 8-)
logistcorrect(:,2) = 1 - (1 ./ (1 + exp(theta + beta * logistcorrect(:,1))));
