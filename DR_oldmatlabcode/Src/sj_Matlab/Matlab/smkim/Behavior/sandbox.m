function correct = scoreCSAperformance(visits)
%
%   scoreCSAperformance(visits)
% identifies correct versus incorrect choices on inbound and
% outbound trials of the W-maze continuous spatial alternation
% task, tallies across all days and sessions

correct = struct;
correct.fields = 'day    epoch    correct?'

% both outbound and inbound are 3-column matrices
% column 1 is day, column 2 is epoch, column 3 is correct (1) or error (0)
% number of rows corresponds to the number of trials
correct.outbound = zeros(0,3);
correct.inbound = zeros(0,3);

% constant definitions of the roi indexes corresponding to choice arms
LEFT = 1;
CENTER = 2;
RIGHT = 3;

for i = 2:length(visits) % numdays
    for j = 1:length(visits{i}) % num sessions per day
        % grab a list of the roi visits made in that session
        session_visits = visits{i}{j}.data(:,1);
        % censor out repeated visits to the same arm
        session_visits = session_visits([session_visits(1) ...
             ; (find(diff(session_visits))+1)]);
        % given that judging performance on W-track CSA requires 
        % information from at most three consecutive trials, we
        % separately analyze the first 3 roi visits of the session
        % before scanning the remainder

        for k = 1:length(session_visits)
            % starting in center arm, examining outbound choice
            if session_visits(k) == CENTER
                % skip this iteration if we are at beginning or
                % end of session
                if (k == 1) || (k == length(session_visits))
                    continue % skip 
                end  
                % is the rat alternating correctly?
                if session_visits(k-1) ~= session_visits(k+1)
                    correct.outbound(end+1,:) = [i j 1];
                % or is the rat repeatedly visiting the same side arm?
                else
                    correct.outbound(end+1,:) = [i j 0];
                end
            % starting in left or right arm, examining inbound choice
            else
                % skip this iteration if we are at the end of session
                if (k == length(session_visits))
                    continue % skip 
                end  
                % is the rat returning to center arm correctly?
                if session_visits(k+1) == CENTER
                    correct.inbound(end+1,:) = [i j 1];
                % or is the rat making an error detour?
                else
                    correct.inbound(end+1,:) = [i j 0];
                end
            end
        end     
    end
end

%{

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
%
% this table is for the W maze:
contingency_table = zeros(length(visits.roi)*ones(1,length(history)));
contingency_table(2,1,:) = 1;
contingency_table(2,3,:) = 1;
contingency_table(3,2,1) = 1;
contingency_table(1,2,3) = 1;

% response bias:
% the bias of R turn response is given by 
bias_Routbound = -1/2*(norminv(PRgivenCR)-norminv(PRgivenCL))

%}
