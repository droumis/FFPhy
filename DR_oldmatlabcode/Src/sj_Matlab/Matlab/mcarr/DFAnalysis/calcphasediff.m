function out = calcphasediff(index, excludeperiods, theta, varargin)
% function out = calcphasediff(index, excludeperiods, theta, varargin)
%
%  Calculates the mean phase difference between theta traces on pairs of
%  tetrodes after excluded times have been removed.  


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

% assign a temporary variable for the theta
t1 = theta{index(1)}{index(2)}{index(3)};
t2 = theta{index(1)}{index(2)}{index(4)};
clear theta;

% the phase is the second field of the data element
p1 = t1.data(:,2);
p2 = t2.data(:,2);

% create a list of times for the theta data
tme = t1.starttime + [0:(length(p1)-1)] / t1.samprate;

% apply the exclude filter
p1 = p1(find(~isExcluded(tme, excludeperiods)));
p2 = p2(find(~isExcluded(tme, excludeperiods)));

out.bins = [-pi:pi/30:pi]';
out.phasediffhist = histc(double(p1-p2) / 10000, out.bins);

if (appendindex)
    out.index = index;
end
