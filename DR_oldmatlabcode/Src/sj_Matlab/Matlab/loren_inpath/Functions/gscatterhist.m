function h = gscatterhist(x,y,g)
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
%   Use the data cursor to read precise values and observation numbers 
%   from the plot.
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

%   Copyright 2006-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2008/02/29 13:12:27 $

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
obsinds = 1:numel(x);
t = isnan(x) | isnan(y);
if any(t)
    x(t) = [];
    y(t) = [];
    obsinds(t) = [];
end





%% Create the histogram information

grpID = grp2idx(g); 
grp   = unique(grpID) ;
nGrp  = length(unique(g)) ;

cxmax = max(x) ;
cymax = max(y) ;

cxmin = min(x) ;
cymin = min(y) ;

for i = 1:nGrp

     dfx{i} = fitdist(x(grpID == grp(i) ), 'kernel') ;
     dfy{i} = fitdist(y(grpID == grp(i) ), 'kernel') ; 

end

dx = 0.1*range(x) ; 
dy = 0.1*range(y) ; 


xrange = [cxmin - dx:0.01*dx: cxmax + dx];
yrange = [cymin - dy:0.01*dy: cymax + dy];

for i = 1:nGrp

     px(:, i) = dfx{i}.pdf(xrange);
     py(:, i) = dfy{i}.pdf(yrange); 

end

%% Put up the plots in preliminary positions
clf
hScatter = subplot(2,2,2);
hScatterline = gscatter(x,y,g, 'bgrcmyk'); 
axis([xlim ylim]);
xlabel('x'); ylabel('y');
grid on;
box on
legend('Location', 'Best')

hHistX = subplot(2,2,4);


hX = plot(xrange ,-px) ;

set(hX, 'LineWidth', 1.5)
axis([min(xrange) max(xrange), -max(px(:)), 0]);
set(hHistX, 'yTick', [])
set(hHistX, 'xTick', [])

hHistY = subplot(2,2,1);

hY = plot( -py , yrange); 

set(hY, 'LineWidth', 1.5)
axis([-max(py(:)), 0, min(yrange), max(yrange)]); 
set(hHistY, 'yTick', [])
set(hHistY, 'xTick', [])

% Make scatter plot bigger, histograms smaller
set(hScatter,'Position',[0.35 0.35 0.55 0.55],'tag','scatter');
set(hHistX,'Position',[.35 .1 .55 .15],'tag','xhist');
set(hHistY,'Position',[.1 .35 .15 .55],'tag','yhist');

% Leave scatter plot as current axes
set(get(hScatter,'parent'),'CurrentAxes',hScatter);

if nargout>0
    h = [hScatter hHistX hHistY];
end

legend('Location', 'NorthEast')

