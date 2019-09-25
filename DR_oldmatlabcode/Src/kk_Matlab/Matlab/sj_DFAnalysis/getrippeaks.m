function [out] = getrippeaks(index, excludeperiods, ripples, varargin)
% out = getriptimes(index, excludeperiods, ripples, options)
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
%   'excludethresh' 
%           exclude values > excludethresh
%
%   out is [peakrates]

% assign the options
appendindex = 0;
average = 0;
excludethresh =[];
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'appendindex'
            appendindex = varargin{option+1};
     case 'average'
            average = varargin{option+1};
     case 'excludethresh'
            excludethresh = varargin{option+1};
      otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

r = ripples{index(1)}{index(2)}{index(3)};
ripindex = ~isExcluded(r.starttime, excludeperiods);
peaks = r.peak(ripindex);
if ~isempty(excludethresh)
    peaks = peaks(find(peaks<excludethresh));
end
if appendindex == 0
    out = peaks;
elseif appendindex ==1
    firstcolumns = repmat(index, length(peaks), 1);
    out = [firstcolumns peaks];
end

if (average)
    out = mean(out, 1);
    out = [ out stderr(peaks)];

end

