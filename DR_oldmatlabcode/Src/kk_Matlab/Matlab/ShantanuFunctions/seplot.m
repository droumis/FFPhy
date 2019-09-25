function h = seplot(x, varargin)
% function h = seplot(x, varargin)
%
% Either called with
% seplot(x, y, sey)
% or
% seplot(y, sey)

if ( (nargin > 1) & ~ischar(varargin{1}) )
   y = varargin{1};
   hh = plot(x,mean(y));
   flag2 = 1;
   if (length(varargin) > 1)
      varargin = {varargin{2:end}};
   else
      varargin = {};
   end
else
   y = x;
   hh = plot(mean(y));
   flag2 = 0;
end


holdState = ishold;
hold on
[lineStyle, lineWidth, color, marker, confalpha, faceAlpha, fcolor, vartype, otherArgs] = ...
process_options(varargin, ...
   'lineStyle', '-', 'lineWidth', 1, 'Color', [], 'Marker', [], ...
   'confalpha', 0.05, 'FaceAlpha', [], 'fcolor', [], 'VarType', 'se');


if isempty(vartype) | strcmpi(vartype,'se')
   se = transpose(std(y)./sqrt(size(y,1) - 1));
elseif strcmpi(vartype,'std')
   se = transpose(std(y));
elseif strcmpi(vartype,'conf')
   se = [];
   se1 = prctile(y,100*(1-confalpha));
   se2 = prctile(y,100*confalpha);
end

for i = 1:length(hh)
   xdata = get(hh(i),'XData');
   ydata = get(hh(i),'YData')';
   longXdata = [xdata, fliplr(xdata)];
   if ~isempty(se)
     longSE = [ydata + se(:,i); flipud(ydata - se(:,i))];
   else
     longSE = [se1(:); flipud(se2(:))];
   end
   if isempty(fcolor)
     fcolor{i} = get(hh(i),'Color');
   end
   if ~isempty(color)
       ff(i) = fill(longXdata,longSE,color,'EdgeColor','none');
   else
        ff(i) = fill(longXdata,longSE,fcolor{i},'EdgeColor','none');
   end
   delete(hh(i));
end

my = mean(y);

if flag2
   if (length(otherArgs) > 0)
      hh = plot(x,my,'linestyle',lineStyle,'linewidth',lineWidth, ...
         otherArgs{:});
   else
      hh = plot(x,my,'linestyle',lineStyle,'linewidth',lineWidth);
   end
else
   if (length(otherArgs) > 0)
      hh = plot(my,'linestyle',lineStyle,'linewidth',lineWidth, ...
         otherArgs{:});
   else
      hh = plot(my,'linestyle',lineStyle,'linewidth',lineWidth);
   end
end

for i = 1:length(hh)
   if ~isempty(color)
      set(hh(i),'Color',color*0.5);
     
   else
      set(hh(i),'Color',fcolor{i}*0.5);
   end
   if ~isempty(marker)
      set(hh(i),'Marker',marker);
   end
   if ~isempty(faceAlpha)
      set(ff(i),'FaceAlpha',faceAlpha);
   end
end

h = [hh; ff'];

if (holdState == 0)
   hold off
end
