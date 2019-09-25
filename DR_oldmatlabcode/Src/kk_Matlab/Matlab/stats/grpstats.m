function [varargout] = grpstats(x,group,whichstats,alpha)
%GRPSTATS Summary statistics by group.
%   MEANS = GRPSTATS(X,GROUP) returns the MEANS of each column of X by
%   GROUP. X is a matrix of observations.  GROUP is a grouping variable
%   defined as a categorical variable, vector, string array, or cell array
%   of strings.  GROUP can also be a cell array of several grouping
%   variables (such as {G1 G2 G3}) to group the values in X by each unique
%   combination of grouping variable values.  GROUP can be [] or omitted to
%   compute the mean of the entire sample without grouping.  When there is
%   a single grouping variable, groups are sorted by order of appearance
%   (if GROUP is character), sorted numeric value (if GROUP is numeric), or
%   order of the levels property (if GROUP is categorical).
%
%   GRPSTATS(X,GROUP,ALPHA) displays a plot of the means versus index with
%   100(1 - ALPHA)%  confidence intervals around each mean.
%
%   [A,B,...] = GRPSTATS(X,GROUP,WHICHSTATS} returns the statistics specified
%   in WHICHSTATS.  WHICHSTATS can be a single function handle or name, or a
%   cell array containing multiple function handles or names.  The number of
%   outputs (A,B,...) must match the number function handles and names in
%   WHICHSTATS.  The names can be chosen from among the following:
%
%      'mean'     mean
%      'sem'      standard error of the mean
%      'numel'    count, or number of elements
%      'gname'    group name
%      'std'      standard deviation
%      'var'      variance
%      'meanci'   95% confidence interval for the mean
%      'predci'   95% prediction interval for a new observation
%
%   Each function included in WHICHSTATS must accept a vector of data and
%   compute a descriptive statistic for it.  For example, @median and
%   @skewness are suitable functions.  The size of the output is
%   NGROUPS-by-NVALS, where NGROUPS is the number of groups and NVALS is
%   the number of values returned by the function for a single group.  The
%   function may also accept a matrix of data and compute column
%   statistics.  In this case the output size is NGROUPS-by-NCOLS-by-NVALS,
%   where NCOLS is the number of columns of X.
%
%   [...] = GRPSTATS(X,GROUP,WHICHSTATS,ALPHA) specifies the confidence
%   level as 100(1-ALPHA)% for the 'meanci' and 'predci' options.  It does
%   not display a plot.
%
%   Example:  
%      load carsmall
%      [m,p,g] = grpstats(Weight,Model_Year,{'mean','predci','gname'})
%      n = length(m)
%      errorbar((1:n)',m,p(:,2)-m)
%      set(gca,'xtick',1:n,'xticklabel',g)
%      title('95% prediction intervals for mean weight by year')
%
%   See also GSCATTER, GRP2IDX.

%   Older syntax still supported:
%   [MEANS,SEM,COUNTS,GNAME] = GRPSTATS(X,GROUP) returns the standard error
%   of the mean in SEM, the number of elements in each group in COUNTS,
%   and the name of each group in GNAME.

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 2.15.2.12 $  $Date: 2006/11/11 22:55:12 $

if (nargin<1)
   error('stats:grpstats:TooFewInputs',...
         'GRPSTATS requires at least one argument.')
end
if ndims(x)>2 || ~isreal(x)
    error('stats:grpstats:BadData',...
          'X must be a vector or matrix of real numbers.');
elseif isvector(x)
    x = x(:);
end
[row,cols] = size(x);

% Recognize plotting syntax with alpha in 3rd position
doplot = false;
if nargin==3
    if isnumeric(whichstats) && isscalar(whichstats)
        alpha = whichstats;
        whichstats = {};
        doplot = true;
    else
        alpha = 0.05;
    end
elseif nargin<=2
    alpha = 0.05;   % used in nested functions
    whichstats = {};
end

if alpha<=0 || alpha>=1
    error('stats:grpstats:BadAlpha',...
          'ALPHA must be a number larger than 0 and smaller than 1.');
end

% Get list of statistics functions to call
if isempty(whichstats)
    % Default list
    whichstats = {@(x)localmean(x), @sem, @(x)size(x,1), 'gname'};
    if doplot
        minargs = 3;
    else
        minargs = 1;
    end
    whichstats = whichstats(1:max(minargs,nargout));
else
    if ~iscell(whichstats)
        whichstats = {whichstats};
    end

    % Convert keywords to function handles
    for j=1:numel(whichstats)
        hfun = whichstats{j};
        if ischar(hfun)
            switch(hfun)
              case 'mean',  hfun = @localmean;
              case 'sem',   hfun = @sem;
              case 'numel', hfun = @(x)size(x,1);
              case 'std',   hfun = @localstd;
              case 'var',   hfun = @localvar;
              case 'meanci',hfun = @meanci;
              case 'predci',hfun = @predci;
              %otherwise, may be a function name or 'gname'
            end
        whichstats{j} = hfun;
        end
    end

    if max(1,nargout)~=numel(whichstats)
        warning('stats:grpstats:ArgumentMismatch',...
                'Number of outputs does not match number of statistics.')
    end
end

% Get grouping variable information
if (nargin<2)
   group = [];
end
if isempty(group)
   group = ones(row,1);
end
[group,glabel,groupname,multigroup,ngroups] = mgrp2idx(group,row);
if length(group) ~= size(x,1)
    error('stats:grpstats:InputSizeMismatch',...
          'Must have one GROUP for each row of X.');
end

% Collect group information
groups = cell(1,ngroups);
for gnum = 1:ngroups
    groups{gnum} = find(group==gnum);
end

nfuns = numel(whichstats);
varargout = cell(1,max(1,nfuns));

for nfun = 1:nfuns
    hfun = whichstats{nfun};   % get function handle or name
    if isequal(hfun,'gname')
        % special case for gname, not applied separately to each column
        varargout{nfun} = groupname;
        continue
    end

    % Should we try to apply the function to an entire matrix or just a column?
    trymatrix = (cols~=1) && ~any(isnan(x(:)));

    % Test the function to see what we get
    if isempty(groups)
        rowidx = [];
    else
        rowidx = groups{1};
    end

    if trymatrix
        % Attempt to call the function on a data matrix
        try
            t = feval(hfun,x(rowidx,:));
            if size(t,2)~=cols
                trymatrix = false;
            end
        catch
            trymatrix = false;
        end
    end
    if trymatrix
        % Success, put results for this group into an array
        nstatvals = size(t,1);
        t1 = reshape(t',[1,cols,nstatvals]);   % 1st dim for groups
        z = repmat(t1,[ngroups,1,1]);          % one per group
        tsize = [nstatvals,cols];
    else    
        % Call the function on one column
        if size(x,2)>=1
            y = x(rowidx,1);
        else
            y = x(rowidx,[]);
        end
        t = tryeval(hfun,y(~isnan(y),:));
        nstatvals = size(t,1);
        t1 = reshape(t,[1,1,nstatvals]);       % dims 1-2 for group,col
        z = repmat(t1,[ngroups,cols,1]);       % one per group and col
        tsize = size(t);
        if ngroups>0 && cols>0
            for colnum = 2:cols
                % Now do the rest of the columns
                y = x(rowidx,colnum);
                z(1,colnum,:) = tryeval(hfun,y(~isnan(y),:),tsize);
            end
        end
    end

    % Now do the rest of the groups
    for gnum = 2:ngroups
        idx = groups{gnum};

        if trymatrix
            z(gnum,:,:) = tryeval(hfun,x(idx,:),tsize)';
        else
            for colnum = 1:cols
                y = x(idx,colnum);
                z(gnum,colnum,:) = tryeval(hfun,y(~isnan(y),:),tsize);
            end
        end
    end
    
    % Special case:  don't add 3rd dimension of there is just one column
    if cols==1
        z = reshape(z,ngroups,nstatvals);
    end
    varargout{nfun} = z;
end

if doplot
   means = varargout{1};
   sems = varargout{2};
   counts = varargout{3};
   p = 1 - alpha/2;
   xd = repmat((1:ngroups)',1,cols);
   h = errorbar(xd,means,tinv(p,counts-1) .* sems);
   set(h,'Marker','o','MarkerSize',2);
   set(gca,'Xlim',[0.5 ngroups+0.5],'Xtick',(1:ngroups));
   xlabel('Group');
   ylabel('Mean');
   if (multigroup)
      % Turn off tick labels and axis label
      set(gca, 'XTickLabel','','UserData',size(groupname,2));
      xlabel('');
      ylim = get(gca, 'YLim');
      
      % Place multi-line text approximately where tick labels belong
      for j=1:ngroups
         text(j,ylim(1),glabel{j,1},'HorizontalAlignment','center',...
              'VerticalAlignment','top', 'UserData','xtick');
      end

      % Resize function will position text more accurately
      set(gcf, 'ResizeFcn', @resizefcn, 'Interruptible','off');
      doresize(gcf);
   else
      set(gca, 'XTickLabel',glabel);
   end
   title('Means and Confidence Intervals for Each Group');
   set(gca, 'YGrid', 'on');
end

% Nested functions below here; they use alpha from caller
    function ci = meanci(y,m,s,n,d) % m,s,n,d are local variables
    n = size(y,1);
    if n<=1
        ci = nan(2,size(y,2));
    else
        m = mean(y,1);
        s = std(y,0,1) / sqrt(n);
        d = s * -tinv(alpha/2, max(0,n-1));
        ci = [m-d; m+d];
    end
    end

    % ----------------------------
    function ci = predci(y,m,s,n,d) % m,s,n,d are local variables
    n = size(y,1);
    if n<=1
        ci = nan(2,size(y,2));
    else
        m = mean(y,1);
        s = std(y,0,1) * sqrt(1 + 1/n);
        d = s * -tinv(alpha/2, max(0,n-1));
        ci = [m-d; m+d];
    end
    end
end

% Other local descriptive statistics functions


function t = localmean(y)
n = size(y,1);
if n==0
    t = nan(1,size(y,2));
else
    t = mean(y,1);
end
end

function t = localstd(y)
t = sqrt(localvar(y));
end


function t=localvar(y)
n = size(y,1);
if n<=1
    t = nan(1,size(y,2));
else
    t = var(y,0,1);
end
end


function t = sem(y)
n = size(y,1);
if n==0
    t = nan(1,size(y,2));
else
    t = std(y,0,1) / sqrt(n);
end
end

% ----------------------------
function t = tryeval(f,y,tsize)
errtype = 0;
try
    t = feval(f,y);
catch
    errtype = 1;
end
if nargin>=3 && ~isempty(tsize) && ~isequal(tsize,size(t))
    errtype = 2;
end
if errtype>0
    if ischar(f)
        fname = f;
    else
        fname = func2str(f);
    end
    if errtype==1    
        error('stats:grpstats:FunctionError', ...
              'Error evaluating function %s:\n%s',fname,lasterr);
    else
        error('stats:grpstats:BadFunctionResult', ...
            'Function %s returned result of size [%s], expected size [%s].',...
              fname,num2str(size(t)),num2str(tsize));
    end
end

end

% ----------------
function resizefcn(varargin)
% Resize callback
doresize(gcbf);
end

% -------------------------
function doresize(f)
% Adjust figure layout to make sure labels remain visible
h = findobj(f, 'UserData','xtick');
if (isempty(h))
   set(f, 'ResizeFcn', '');
   return;
end
ax = get(f, 'CurrentAxes');
nlines = get(ax, 'UserData');

% Position the axes so that the fake X tick labels have room to display
set(ax, 'Units', 'characters');
p = get(ax, 'Position');
ptop = p(2) + p(4);
if (p(4) < nlines+1.5)
   p(2) = ptop/2;
else
   p(2) = nlines + 1;
end
p(4) = ptop - p(2);
set(ax, 'Position', p);
set(ax, 'Units', 'normalized');

% Position the labels at the proper place
xl = get(gca, 'XLabel');
set(xl, 'Units', 'data');
p = get(xl, 'Position');
ylim = get(gca, 'YLim');
p2 = (p(2)+ylim(1))/2;
for j=1:length(h)
   p = get(h(j), 'Position') ;
   p(2) = p2;
   set(h(j), 'Position', p);
end
end