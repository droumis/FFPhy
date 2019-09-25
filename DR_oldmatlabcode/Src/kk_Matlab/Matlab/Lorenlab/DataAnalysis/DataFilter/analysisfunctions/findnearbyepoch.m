function [out] = findpreceedingepoch(index, task, varargin)
% function [out] = findpreceedingsleep(index, task, varargin)

epochfield = 'type';
epochvalue = 'sleep';
direction = 'before';
[otherArgs] = procOptions(varargin);

out = [];
if isempty(index)
  return;
end

if strcmp(direction,'fullybefore')
  index = index(1,:);
  direction = 'before';
elseif strcmp(direction,'fullyafter')
  index = index(end,:);
  direction = 'after';
end

out = [];
k = 1;
for i = 1:size(index,1)
  flag = 0;
  if strcmp(direction,'before')
    for j = index(i,2)-1:-1:1
      if strcmp(task{index(i,1)}{j}.(epochfield),epochvalue)
        out(k,:) = [index(i,1) j];
        flag = 1;
        k = k + 1;
      elseif flag==1
        break;
      end
    end
  elseif strcmp(direction,'after')
    for j = index(i,2)+1:1:length(task{index(i,1)})
      if strcmp(task{index(i,1)}{j}.(epochfield),epochvalue)
        out(k,:) = [index(i,1) j];
        flag = 1;
        k = k + 1;
      elseif flag==1
        break;
      end
    end
  end
end

