function h = gscatter(x,y,g,clr,sym,siz,doleg,xnam,ynam)
%GSCATTER   Scatter plot with grouping variable
%   GSCATTER(X,Y,G) creates a scatter plot of the vectors X and Y grouped
%   by G.  Points with the same value of G are shown with the same color
%   and marker.  G is a grouping variable defined as a categorical
%   variable, vector, cell array of strings, or string matrix, and it must
%   have the same number of rows as X and Y.  Alternatively G can be a cell
%   array of grouping variables (such as {G1 G2 G3}) to group the values in
%   X by each unique combination of grouping variable values.
%
%   GSCATTER(X,Y,G,CLR,SYM,SIZ) specifies the colors, markers, and
%   size to use.  CLR is either a string of color specifications or
%   a three-column matrix of color specifications.  SYM is a string
%   of marker specifications.  Type "help plot" for more information.
%   For example, if SYM='o+x', the first group will be plotted with a
%   circle, the second with plus, and the third with x.  SIZ is a
%   marker size to use for all plots.  By default, the marker is '.'.
%
%   GSCATTER(X,Y,G,CLR,SYM,SIZ,DOLEG) lets you control whether legends
%   are created.  Set DOLEG to 'on' (default) or 'off'.
%
%   GSCATTER(X,Y,G,CLR,SYM,SIZ,DOLEG,XNAM,YNAM) specifies XNAM and
%   YNAM as the names of the X and Y variables.  Each must be a
%   character string.  If you omit XNAM and YNAM, GSCATTER attempts to
%   determine the names of the variables passed in as the first and
%   second arguments.
%
%   H = GSCATTER(...) returns an array of handles to the objects
%   created.
%
%   Example:  Scatter plot of car data coded by country.
%      load carsmall
%      gscatter(Weight, MPG, Origin)
%
%   See also GRPSTATS, GRP2IDX.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 1.7.2.5 $  $Date: 2006/11/11 22:55:13 $

error(nargchk(2,9,nargin,'struct'));

% Default colors, markers, etc.
if (nargin < 4), clr = ''; end
if (nargin < 5) || isempty(sym), sym = '.'; end
if (nargin < 6), siz = []; end
if (nargin < 7), doleg = 'on'; end
if (nargin < 8), xnam = inputname(1); end
if (nargin < 9), ynam = inputname(2); end

% What should go into the plot matrix?
doleg = strcmp(doleg, 'on');

% Don't plot anything if either x or y is empty
if isempty(x) || isempty(y),
   if nargout>0
       h = [];
   end
   return
end

if (ndims(x)==2) && any(size(x)==1), x = x(:); end
if (ndims(y)==2) && any(size(y)==1), y = y(:); end

if ndims(x)>2 || ndims(y)>2
   error('stats:gscatter:MatrixRequired',...
         'X and Y must be 2-D.');
end
if size(x,1)~=size(y,1)
   error('stats:gscatter:InputSizeMismatch',...
         'X and Y must have the same length.');
end

if (nargin > 2)
   [g,gn,ignore1,ignore2,maxgrp] = mgrp2idx(g,size(x,1),',');
   ng = max(g);
else
   g = [];
   gn = [];
   ng = 1;
   maxgrp = 1;
end

if (length(g)>0) && (length(g) ~= size(x,1)),
   error('stats:gscatter:InputSizeMismatch',...
         'There must be one value of G for each row of X.');
end

if (isempty(siz))
   siz = repmat(get(0, 'defaultlinemarkersize'), size(sym));
   if any(sym=='.'),
      units = get(gcf,'units');
      set(gcf,'units','pixels');
      pos = get(gcf,'Position');
      set(gcf,'units',units);
      siz(sym=='.') = max(1,min(15, round(15*min(pos(3:4))/size(x,1))));
   end
end

newplot;
if isempty(clr)
    clr = hsv(maxgrp);
end
hh = iscatter(x, y, g, clr, sym, siz);

% Label plots
if (~isempty(xnam)), xlabel(deblank(xnam)); end
if (~isempty(ynam)), ylabel(deblank(ynam)); end

% Create legend if requested
if (doleg && ~isempty(gn))
   t = find(ismember(1:size(gn,1),g));
   legend(hh(t), gn(t,:));
end

% Nudge X axis limits if points are too close
xlim = get(gca, 'XLim');
d = diff(xlim);
xlim(1) = min(xlim(1), min(min(x))-0.05*d);
xlim(2) = max(xlim(2), max(max(x))+0.05*d);
set(gca, 'XLim', xlim);

if (nargout>0), h = hh; end

% Store information for gname
set(gca, 'UserData', {'gscatter' x y g});
