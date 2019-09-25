function out = calcperformance(index, excludetimes, linpos, task, correctorder, varargin)
%%assumes animal is performing an alternation task, with a middle well and
%2, alternating outside wells, though it does not have to be a W-shaped
%track.  does not track trajectories, as animal may wander between time
%leaves one well and enters another.
%
% LINPOS is the output of linearizeposition for one day
% INDEX [day epoch]
% EXCLUDETIMES times to exclude from analysis
% CORRECTORDER is a vector of 3 well numbers, the middle number corresponds
% to the middle well, the first and third ones correspond to outside wells,
% for instance [2 1 3] would mean the middle well is well 1, the outside
% wells are 2 and 3.  these well numbers are determined in createtaskstruct
%
% out is proportion rewarded of included trajectories
%
% options
%   'includetrajbound', default includes all
%       [-1 0 1], matrix of all trajbounds to include, -1 for nontask-bound
%       trajectories, 0 for outbound trajectories (to outside task arms), 1
%       for inbound trajectories
%   'appendindex', 0 or 1, default 0

appendindex = 0;
incltrajb = [];
appendtaskinfo = 1;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'appendtaskinfo'
                appendtaskinfo = varargin{option+1};
            case 'includetrajbound'
                inctrajb = varargin{option+1};
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

if isempty(correctorder)
    correctorder = task{index(1,1)}{index(1,2)}.wellseq;
end

[trajbound rewarded time] = trackperformance(index, excludetimes, linpos, correctorder);

if ~isempty(incltrajb)
    rewarded = rewarded(ismember(trajbound,inctrajb));
end

tmpout = sum(rewarded)/length(rewarded);

if appendtaskinfo
    if isfield(task{index(1,1)}{index(1,2)}, 'switchday')
    tmpind = [task{index(1,1)}{index(1,2)}.exposure task{index(1,1)}{index(1,2)}.switchday task{index(1,1)}{index(1,2)}.tasknum];
    else
        tmpind = [task{index(1,1)}{index(1,2)}.exposure];
    end
    tmpout = [tmpout tmpind ];
end

if appendindex == 0
    out = tmpout; % perccorrect exposure experimentday exposureday
else
    out = [index tmpout];
end



end

