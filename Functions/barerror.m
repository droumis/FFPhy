% [h he] = barerror(x, y, e, errorwidth, width, 'grouped')
%       Plots a bar graph of the data in y with error bounds as in e.
%       errorwidth is the width of the error bar tops
%       width and 'grouped" are optional and are passed directly into bar();
%       x should be 1 x N
%       y should be N x M
%       e should be N x M
%           where N is the number of groups (for a single set of bars)
%                   or the number of levels of each factor (for a grouped chart)
%                 M is 1 for a single set of bars
%		    or the number of groups for a grouped chart
%       h is a handle for the bars
%       he is a handle for the errobars
function [h, he] = barerror(x, y, e, ew, varargin)

% check to see if this is a 1 dimensional bar graph
if (min(size(y)) == 1) 
    h = bar(x, y, varargin{:});
    hold on
    he = errorbar2(x, y, e, ew, 'k.');
else
    % multidimensional graph
    h = bar(x, y, varargin{:});
    hold on
    % get the centers of the bars for each group
    for i = 1:length(h)
	xloc = mean(get(h(i), 'XData'));
        he{i} = errorbar2(xloc, y(:,i), e(:,i), ew, 'k.');
    end
end
% set up reasonable limits
% check to see if we have a vector y and a matrix x
if ((min(size(y)) == 1) & (min(size(e)) == 2))
    miny = min([min(y - e(1,:)) min(y + e(2,:))]);
    maxy = max([max(y - e(1,:)) max(y + e(2,:))]);
else
    miny = min(min(y - e));
    maxy = max(max(y + e));
end

if (miny < 0) 
    ymin = miny * 1.1;
elseif (miny >= 0)
    ymin = 0;
end

if (maxy > 0) 
    ymax = maxy * 1.1;
elseif (miny <= 0)
    ymax = maxy;
end

if (length(x) > 1)
    xmin = x(1) - mean(diff(x))/2;
    xmax = x(end) + mean(diff(x))/2;
else
    xmin = x(1) - 1;
    xmax = x(1) + 1;
end

set(gca, 'XTick', x, 'YLim', [ymin ymax], 'XLim', [xmin xmax]);
