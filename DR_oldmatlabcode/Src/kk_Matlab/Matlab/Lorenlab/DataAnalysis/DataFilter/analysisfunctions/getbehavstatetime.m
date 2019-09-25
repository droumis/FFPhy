function [out] = getbehavstatetime(index, times, excludeperiods, pos, linpos, task, track_regions, varargin)
% function [out] = getbehavstatetime(index, times, excludeperiods, pos, linpos, task, track_regions, varargin)
%
%   index [day epoch]
%
%   out is behavior organized as trial x time x tetrode

% assign the options

behavVars = {'position', 'velocity'};
t_offset = 0;

[otherArgs] = procOptions(varargin);

out = [];
if isempty(times) | isempty(index)
  return;
end

if size(index,1)>1
  error('Multiple epochs sent to getbehavstatetime.');
end

if (nargin>4) && ~isempty(linpos) && length(linpos)>=index(1) && ~isempty(linpos{index(1)}) && isfield(linpos{index(1)}{index(2)},'statematrix')
  lp = linpos{index(1)}{index(2)}.statematrix.lindist;
else
  lp = [];
end

if isfield(pos{index(1)}{index(2)},'data')
  p = pos{index(1)}{index(2)}.data;
else
  p = [];
end

if (nargin>5)
  tk = task{index(1)}{index(2)};
else
  tk = [];
end

if ~(nargin>6)
  track_regions = [];
end

goodTimes = ~isExcluded(times(:,1), excludeperiods);
if size(times,2) > 1
  for tt = 1:size(times,2)
    goodTimes = goodTimes & ~isExcluded(times(:,1), excludeperiods);
  end
end
goodTimes = find(goodTimes);

times = times + t_offset;

if ~isempty(p)
  posinds = lookup(times(goodTimes,1),p(:,1));
end

for i = 1:length(behavVars)
  switch behavVars{i}
    case 'posinds'
      if ~isempty(p)
        out = [out posinds];
      else
        out = [out nan(length(goodTimes),1)];
      end
    case 'position'
      if ~isempty(p)
        out = [out p(posinds,2:3)];
      else
        out = [out nan(length(goodTimes),2)];
      end
    case 'velocity'
      if ~isempty(p)
        out = [out p(posinds,5)];
      else
        out = [out nan(length(goodTimes),1)];
      end
    case 'smoothvelocity'
      if ~isempty(p)
        out = [out p(posinds,11)];
      else
        out = [out nan(length(goodTimes),1)];
      end
    case 'direction'
      if ~isempty(p)
        out = [out p(posinds,4)];
      else
        out = [out nan(length(goodTimes),1)];
      end
    case 'quiescence'
      if ~isempty(p)
        out = [out p(posinds,12)];
      else
        out = [out nan(length(goodTimes),1)];
      end
    case 'novelty'
      if ~isempty(p)
        out = [out findNovelty(track_regions,tk,p(posinds,2:3))];
      else
        out = [out nan(length(goodTimes),2)];
      end
    case 'linpos'
      if ~isempty(lp)
        out = [out lp(posinds)];
      else
        out = [out nan(length(goodTimes),1)];
      end
  end
end
