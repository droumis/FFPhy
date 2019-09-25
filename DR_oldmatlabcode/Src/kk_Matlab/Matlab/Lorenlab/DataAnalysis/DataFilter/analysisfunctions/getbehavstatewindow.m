function [out] = getbehavstatetime(index, times, excludeperiods, pos, linpos, task, track_regions, varargin)
% function [out] = getbehavstatetime(index, times, excludeperiods, pos, linpos, task, track_regions, varargin)
%
%   index [day epoch]
%
%   out is behavior organized as trial x time x tetrode

% assign the options

behavVar = 'smoothvelocity';
behavwindow = [0 0];
FS = 60/1.001/2; % standard NTSC frame rate

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

goodTidx = ~isExcluded(times(:,1), excludeperiods);
if size(times,2) > 1
  for tt = 1:size(times,2)
    goodTidx = goodTidx & ~isExcluded(times(:,1), excludeperiods);
  end
end
goodTidx = find(goodTidx);


winInds = floor(behavwindow(1)*FS):ceil(behavwindow(2)*FS);
out = nan(length(goodTidx),length(winInds));


if isempty(p)
  return;
end

for i = 1:length(goodTidx)
  posinds = lookup(times(goodTidx(i),1),p(:,1));
  posinds = posinds + winInds;

  idx = find(posinds > 0);
  posinds = posinds(idx);

  switch behavVar
    case 'posinds'
        out(i,idx) = posinds;
    case 'position-x'
        out(i,idx) = p(posinds,2);
    case 'position-y'
        out(i,idx) = p(posinds,3);
    case 'velocity'
        out(i,idx) = p(posinds,5);
    case 'headdir'
        out(i,idx) = p(posinds,4);
    case 'smoothvelocity'
        out(i,idx) = p(posinds,11);
    case 'direction'
        out(i,idx) = p(posinds,4);
    case 'quiescence'
        out(i,idx) = p(posinds,12);
    case 'novelty'
        out(i,idx) = findNovelty(track_regions,tk,p(posinds,2:3));
    case 'linpos'
      if ~isempty(lp)
        out(i,idx) = lp(posinds);
      end
  end
end

