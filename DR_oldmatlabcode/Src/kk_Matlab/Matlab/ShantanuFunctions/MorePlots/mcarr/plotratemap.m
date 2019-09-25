%    [h] = plotratemap(ratemap, options)
%           options are
%            'fontsize', n    n is the fontsize is point
%            'minrate', r    do not display ratemaps with peak < r (default 0)
%	     'showmax', 0 or 1  1 indicates that only the maximum rate should
%	     			be labeled
function [h ch] = plotratemap(rmap, varargin)

minrate = 0;
fontsize = 12;
showmax = 0;
if (~isempty(varargin))
    assign(varargin{:});
end


% create a single behav
cmap = jet(1024);
%cmap = hot(1024) ./ 1.5;
%cmap = cmap(100:920,:);
%cmap = hot(1024);
cmap(1,:) = 1;
colormap(cmap);
% set up the bounds to make it look good
mx = max(rmap(:));
if (mx > 20)
    minbound = -1;
elseif (mx > 10)
    minbound = -.5;
    rmap(find(rmap == -1)) = -.5;
elseif (mx > 3)
    minbound = -.1;
    rmap(find(rmap == -1)) = -.1;
else
    minbound = -.01;
    rmap(find(rmap == -1)) = -.01;
end

bounds = [minbound mx*.65];
if (mx < .1)
    bounds(2) = 1;
end
% set up the colormap to be white at some negative values
if (mx >= minrate)
    h = imagesc(rmap, bounds);
    ch = colorbar;
    set(ch, 'FontSize', fontsize);
    if (showmax)
	set(ch, 'YTick', floor(mx * .65));
    end
end

axis off;
