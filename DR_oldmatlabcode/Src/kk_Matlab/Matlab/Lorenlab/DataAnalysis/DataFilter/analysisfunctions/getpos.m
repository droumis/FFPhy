function [out] = getpos(index, excludeperiods, pos, varargin)
% function [out] = getpos(index, excludeperiods, pos, varargin)
% 
%   index [day epoc tetrode]
%
%   out is struct array of ripples

% assign the options

[otherArgs] = procOptions(varargin);

out = [];
if isempty(index)
  return;
end

p = pos{index(1)}{index(2)};

excluded = isExcluded(p.data(:,1), excludeperiods);

out = p.data(~excluded,:);

