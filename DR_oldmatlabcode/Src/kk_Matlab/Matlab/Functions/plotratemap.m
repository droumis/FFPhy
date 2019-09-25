%    [h] = plotratemap(ratemap, options)
%           options are
%            'fontsize', n    n is the fontsize is point
%            'peakrate', r    set the peak of the colormap to r
%            'minrate', r    do not display ratemaps with peak < r (default 0)
%	     'showmax', 0 or 1  1 indicates that only the maximum rate should
%	     			be labeled
function [h ch] = plotratemap(rmap, varargin)

minrate = 0;
peakrate = [];
fontsize = 12;
showmax = 0;
if (~isempty(varargin))
    assign(varargin{:});
end

if (isempty(rmap))
    return;
end


% create a single behav
%cmap = jet(1024) ./ 1.5;
%cmap = hot(1024) ./ 1.5;
%cmap = cmap(100:920,:);
cmap = jet(1024);
cmap(1,:) = 1;
colormap(cmap);
% set up the bounds to make it look good
bounds = [0 0];
if (isempty(peakrate))
    peakrate = max(rmap(:));
    bounds(2) = peakrate * 0.65;
else
    bounds(2) = peakrate;
end
if (peakrate > 20)
    minbound = -1;
elseif (peakrate > 10)
    minbound = -.5;
    rmap(find(rmap == -1)) = -.5;
elseif (peakrate > 3)
    minbound = -.1;
    rmap(find(rmap == -1)) = -.1;
else
    minbound = -.01;
    rmap(find(rmap == -1)) = -.01;
end

bounds(1) = minbound;
if (peakrate < .1)
    bounds(2) = 1;
end
% set up the colormap to be white at some negative values
if (peakrate >= minrate)
    h = imagesc(rmap, bounds);
    ch = colorbar;
    set(ch, 'FontSize', fontsize);
    if (showmax)
	set(ch, 'YTick', floor(peakrate * .65));
    end
end

axis off;
axis equal
