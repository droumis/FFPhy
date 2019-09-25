
% we define a history of roi visits to be a list of the roi which the 
% animal has visited
% history(1) = current/most recent
% history(2) = last
% history(3) = the one before the last
% ... history is as long as we need it to be for the behavioral task
% for the W maze, length(history)==3
history = zeros(3,1);

% create a roi transition table which tells us which histories correspond to
% correct task performance (1) and which ones are incorrect (0). this table
% is a multi-dimensional array with 
% size(transition_table)==length(roi)*ones(length(history),1) 
contingency_table = zeros(length(visits.roi)*ones(1,length(history)));
contingency_table(2,1,:) = 1;
contingency_table(2,3,:) = 1;
contingency_table(3,2,1) = 1;
contingency_table(1,2,3) = 1;

% response bias:
% the bias of R turn response is given by 
% bias_Routbound = -1/2*(norminv(PRgivenCR)-norminv(PRgivenCL))
