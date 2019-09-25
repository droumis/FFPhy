function [out] = findpreceedingepoch(index, task, varargin)
% function [out] = findpreceedingsleep(index, task, varargin)

epochfield = 'type';
epochvalue = 'sleep';
[otherArgs] = procOptions(varargin);

out = [];
if isempty(index)
  return;
end

out = nan(size(index));

for i = 1:size(index,1)
  for j = index(i,2)-1:-1:1
    if strcmp(task{index(i,1)}{j}.(epochfield),epochvalue)
      out(i,:) = [index(i,1) j];
      break;
    end
  end
end

