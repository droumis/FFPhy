function [out] = getripdurations(index, excludeperiods, ripples, varargin)
% out = getripdurations(index, excludeperiods, ripples, options)
%
%   index [day epoc tetrode]
%
%   options are
%   'average' , 1 or 0, default 0
%           1: average all values together and report mean and standard error
%           0: report all values
%   'appendindex' , 1 or 0, default 0
%           set to 1 to append the cell index to the output [day epoch
%           value]
%
%   out is [durations]

% assign the options
appendindex = 0;
average = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'appendindex'
            appendindex = varargin{option+1};
     case 'average'
            average = varargin{option+1};
      otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

r = ripples{index(1)}{index(2)}{index(3)};
ripindex = ~isExcluded(r.starttime, excludeperiods);
durations = r.endtime(ripindex)-r.starttime(ripindex);


if appendindex == 0
    out = durations;
elseif appendindex ==1
    firstcolumns = repmat(index, length(durations), 1);
    out = [firstcolumns durations];
end

if (average)
    out = mean(out, 1);
    out = [ out stderr(durations)];

end

