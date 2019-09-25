function [out] = getincludedtimes(excludeperiods, varargin)
%function [out] = getincludedtimes(excludeperiods)
%   This function takes excludeperiods and turns it into a n x 2 vector
%   which define the start and end time of included periods.

%Assign the options
starttime = -1;

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'starttime'
            starttime = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%Create a vector of included times
temp = excludeperiods;
temp = reshape(temp', 2*size(excludeperiods,1),1);
temp(1) = [];
temp(length(temp)) = [];
includeperiods = reshape(temp,2,size(excludeperiods,1)-1)';

if starttime ~= -1 && starttime <= excludeperiods(1,1)
    includeperiods = [starttime excludeperiods(1,1); includeperiods];
end

out = includeperiods;