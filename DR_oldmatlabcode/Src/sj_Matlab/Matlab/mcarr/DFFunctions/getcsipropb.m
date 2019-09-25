function out = calccsipropb(ind, excludeperiods, spikes, varargin)
% function out = calccsipropb(index, excludeperiods, spikes, varargin)
% 
% Calculates the csi and probability that a spike is part of a burst 
% Excluded time periods are not included in calculation.
%
% Options:
%   'burstlen', burstlength    specify the burst length to use in ms 
%                              (default 10 ms)
%   'appendindex', 1 or 0 -- set to 1 to append the cell index to the
%   output [tetrode cell value].  Default 0.
%

appendindex = 0;
burstlength = 10;
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'burstlength'
                burstlength = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% create a filtered set of spikes and get their amplitudes
s = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)};
goodspikes = ~isExcluded(s.data(:,1), excludeperiods);
[csi propb] = computesci(s.data(goodspikes(:,1), s.data(goodspikes(:,6),...
			burstlength);


if (appendindex)
    out = [index csi propb];
else
    out = [csi propb];
end
