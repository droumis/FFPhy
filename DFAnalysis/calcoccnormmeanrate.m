function out = calcoccnormmeanrate(index, excludeperiods, spikes, linpos, varargin)
% out = calcoccnormmeanrate(index, excludeperiods, spikes, linpos, options)
% Calculates the mean rate of a cell for a given epoch as the average occupancy-
% normalized linearized rate across spatial bins .  Excluded time periods are not 
% included in the total time.
%
% Options:
%   'appendindex', 1 or 0 -- set to 1 to append the cell index to the
%   output [tetrode cell value].  Default 0.
%

appendindex = 0;
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

linfields = filtercalclinfields(index, excludeperiods, spikes, linpos);
rates = [];
for i = 1:length(linfields.trajdata)
    rates = [rates; linfields.trajdata{i}(find(~isnan(linfields.trajdata{i}(:,5))),5)];
end

if (appendindex)
    out = [index mean(rates)];
else
    out = mean(rates);
end
        