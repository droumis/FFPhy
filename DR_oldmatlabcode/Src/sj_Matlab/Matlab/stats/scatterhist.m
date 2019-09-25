function h = scatterhist(x,y,nbins)
%SCATTERHIST 2D scatter plot with marginal histograms.
%   SCATTERHIST(X,Y) creates a 2D scatterplot of the data in the vectors X
%   and Y, and puts a univariate histogram on the horizontal and vertical
%   axes of the plot.  X and Y must be the same length.
%
%   SCATTERHIST(X,Y,NBINS) also accepts a two-element vector specifying
%   the number of bins for the X and Y histograms.  The default is to
%   compute the number of bins using Scott's rule based on the sample
%   standard deviation.
%
%   Any NaN values in either X or Y are treated as missing data, and are
%   removed from both X and Y.  Therefore the plots reflect points for
%   which neither X nor Y has a missing value.
%
%   H = SCATTERHIST(...) returns a vector of three axes handles for the
%   scatterplot, the histogram along the horizontal axis, and the histogram
%   along the vertical axis, respectively.
%
%   Example:
%      Independent normal and lognormal random samples
%         x = randn(1000,1);
%         y = exp(.5*randn(1000,1));
%         scatterhist(x,y)
%      Marginal uniform samples that are not independent
%         u = copularnd('Gaussian',.8,1000);
%         scatterhist(u(:,1),u(:,2))
%      Mixed discrete and continuous data
%         cars = load('carsmall');
%         scatterhist(cars.Weight,cars.Cylinders,[10 3])

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.1.2.1 $  $Date: 2007/01/19 03:16:22 $

% Check inputs
error(nargchk(2, 3, nargin, 'struct'))

if ~isvector(x) || ~isnumeric(x) || ~isvector(y) || ~isnumeric(y)
    error('stats:scatterhist:BadXY', ...
          'Both X and Y must be numeric vectors.');
end
if numel(x)~=numel(y)
    error('stats:scatterhist:BadXY','X and Y must have the same length.');
end
x = x(:);
y = y(:);
t = isnan(x) | isnan(y);
if any(t)
    x(t) = [];
    y(t) = [];
end

if nargin < 3 || isempty(nbins)
    % By default use the number of bins given by Scott's rule
    xctrs = dfhistbins(x);
    yctrs = dfhistbins(y);
    if length(xctrs)<2
        xctrs = 1;  % bin count 1 for constant data
    end
    if length(yctrs)<2
        yctrs = 1;
    end
elseif ~isnumeric(nbins) || numel(nbins)~=2 || ...
       any(nbins<=0)     || any(nbins~=round(nbins))
    error('stats:scatterhist:BadBins',...
          'NBINS must be a vector of two positive integers.');
else
    xctrs = nbins(1); % use nbins in place of bin centers
    yctrs = nbins(2);
end

% Create the histogram information
[nx,cx] = hist(x,xctrs);
if length(cx)>1
    dx = diff(cx(1:2));
else
    dx = 1;
end
xlim = [cx(1)-dx cx(end)+dx];

[ny,cy] = hist(y,yctrs);
if length(cy)>1
    dy = diff(cy(1:2));
else
    dy = 1;
end
ylim = [cy(1)-dy cy(end)+dy];

yoff = 0;
if prod(ylim)<0, yoff = min(y)*2; end

% Put up the plots in preliminary positions
clf
hScatter = subplot(2,2,2);
plot(x,y,'o'); 
axis([xlim ylim]);
xlabel('x'); ylabel('y');

hHistX = subplot(2,2,4);
bar(cx,-nx,1);
if nx==0
    nx = 1;
end
axis([xlim, -max(nx), 0]);
axis('off');

hHistY = subplot(2,2,1);
barh(cy-yoff,-ny,1); 
if ny==0
    ny = 1;
end
axis([-max(ny), 0, ylim-yoff]); 
axis('off');
line([0 0],ylim-yoff,'Color','k'); % 

% Make scatter plot bigger, histograms smaller
set(hScatter,'Position',[0.35 0.35 0.55 0.55],'tag','scatter');
set(hHistX,'Position',[.35 .1 .55 .15],'tag','xhist');
set(hHistY,'Position',[.1 .35 .15 .55],'tag','yhist');

colormap([.8 .8 1]); % more pleasing histogram fill color

% Leave scatter plot as current axes
set(get(hScatter,'parent'),'CurrentAxes',hScatter);

if nargout>0
    h = [hScatter hHistX hHistY];
end
