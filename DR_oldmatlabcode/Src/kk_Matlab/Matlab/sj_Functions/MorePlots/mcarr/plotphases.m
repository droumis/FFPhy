function out = plotphases(index, excludeperiods, theta, varargin)
% out = plotphases(index, excludeperiods, theta, options)
% plots theta phase for each epoch for each tetrode
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




if (appendindex)
    out = [index max(rates)];
else
    out = max(rates);
end
