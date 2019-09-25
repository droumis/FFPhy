
function out = DFA_loadfields_DR(index, excludetimes, linfields, mapfields, cellinfo)
%just load in the premade linfield and mapfield data for each cell
%create linfields and mapfields with DFS_DR_SaveLinMapFields.m

% if ~~isempty(excludetimes)
%     excludetimes = [];
% end

out.area = cellinfo{index(1)}{index(2)}{index(3)}{index(4)}.area;%get tag from cell info
out.index = index;
try
    out.trajdata = linfields{index(1)}{index(2)}{index(3)}{index(4)};
catch
    out.trajdata = NaN;
end
try
    out.mapdata = mapfields{index(1)}{index(2)}{index(3)}{index(4)};
catch
    out.mapdata = NaN;
end


%this was stupid of me 

% %have to do this incrementallly bc error if structure before the 'cell' does not exist 
% if ~isempty(linfields{index(1)});
%     if ~isempty(linfields{index(1)}{index(2)});
%         if ~isempty(linfields{index(1)}{index(2)}{index(3)});
%             if ~isempty(linfields{index(1)}{index(2)}{index(3)}{index(4)});
%                 out.trajdata = linfields{index(1)}{index(2)}{index(3)}{index(4)};
%             else
%                 out.trajdata = NaN;
%             end
%         else
%             out.trajdata = NaN;
%         end
%     else
%         out.trajdata = NaN;
%     end
% else
%     out.trajdata = NaN;
% end
% 
% 
% if ~isempty(mapfields{index(1)});
%     if ~isempty(mapfields{index(1)}{index(2)});
%         if ~isempty(mapfields{index(1)}{index(2)}{index(3)});
%             if ~isempty(mapfields{index(1)}{index(2)}{index(3)}{index(4)});
%                 out.mapdata = mapfields{index(1)}{index(2)}{index(3)}{index(4)};
%             else
%                 out.mapdata = NaN;
%             end
%         else
%             out.mapdata = NaN;
%         end
%     else
%         out.mapdata = NaN;
%     end
% else
%     out.mapdata = NaN;
% end