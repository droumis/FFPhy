function out = calcpeakrate(index, excludeperiods, spikes, linpos, varargin)
% out = calcpeakrate(index, excludeperiods, spikes, linpos, options)
% Calculates the peak rate of a cell for a given epoch over all bins of the occupancy-
% normalized linearized rate.  Excluded time periods are not 
% included in calculation.
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
    if ~isempty(linfields.trajdata{i})
        rates = [rates; max(linfields.trajdata{i}(find(~isnan(linfields.trajdata{i}(:,5))),5))];
    end
end

if (appendindex)
    out = [index max(rates)];
else
    out = max(rates);
end