function [out] = calcconcurrentrip(index, excludeperiods, ripples, cellinfo, varargin)
% out = getriptimes(index, excludeperiods, ripples, cellinfo, options)
%
%   index [day epoch tet]: list of tetrodes to include in the analysis
%
%   options are
%	'minthresh', a vector of min thresholds in std for ripples to include
%      Default: [3 5 7]

%   'tetfilter1' a string that is used to identify tetrodes in one group,
%      Default: 'isequal($area,''CA1'') & isequal($hemisphere,''right''')'

%   'tetfilter2' a string to identify tetrodes in the second group
%      Default: 'isequal($area,''CA1'') & isequal($hemisphere,''left'')'
%
%   out is Nx2 vector where N is the length of 'std'. First column is the
%   percent_concurrent, second column is the number of ripples detected at
%   that threshold
%   

% assign the options
minthresh = [3 5 7];
tetfilter1 = 'isequal($area,''CA1'') & isequal($hemisphere,''right'')';
tetfilter2 = 'isequal($area,''CA1'') & isequal($hemisphere,''left'')';

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'minthresh'
            minthresh = varargin{option+1};
        case 'tetfilter1'
            tetfilter1 = varargin{option+1};
        case 'tetfilter2'
            tetfilter2 = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%Use tetfilter1 and tetfilter2 to determine two groups of tetrodes
tet1 = evaluatefilter(cellinfo{index(1,1)}{index(1,2)},tetfilter1);
tet1 = unique(tet1(:,1));

tet2 = evaluatefilter(cellinfo{index(1,1)}{index(1,2)},tetfilter2);
tet2 = unique(tet2(:,2));

%Determine percent of SWRs detected in each group of tetrodes
r = zeros(size(minthresh,2),2);
count = 1;

for i = minthresh
    rip1 = getripples(index(1,[1 2]),ripples,cellinfo,'tetrodes',tet1,'minstd',i,'excludeperiods',excludeperiods);
    rip2 = getripples(index(1,[1 2]),ripples,cellinfo,'tetrodes',tet2,'minstd',i,'excludeperiods',excludeperiods);
    
    if isempty(rip1) | isempty(rip2)
        r(count,:) = NaN;
    elseif size(rip1,1) < size(rip2,1)
        tmp = periodAssign(rip1(:,3),rip2(:,[1 2]));
        r(count,:) = [sum(tmp>0)./length(tmp) size(rip1,1)];
    elseif size(rip2,1)<size(rip1,1)
        tmp = periodAssign(rip2(:,3),rip1(:,[1 2]));
        r(count,:) = [sum(tmp>0)./length(tmp) size(rip2,1)];
    else
        tmp1 = periodAssign(rip1(:,3),rip2(:,[1 2]));
        tmp2 = periodAssign(rip2(:,3),rip1(:,[1 2]));
        tmp = [sum(tmp1>0)./length(tmp1) sum(tmp2>0)./length(tmp2)];
        r(count,:) = [max(tmp) size(rip1,1)];
    end
    count = count+1;
end

out.concurrent = r(:,1);
out.nrips = r(:,2);
out.minthresh = minthresh;
end

