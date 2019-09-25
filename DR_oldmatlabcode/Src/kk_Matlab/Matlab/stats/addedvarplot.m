function addedvarplot(x,y,vnum,in,stats,f)
%ADDEDVARPLOT Create added-variable plot for stepwise regression
%   ADDEDVARPLOT(X,Y,VNUM,INMODEL) produces an added variable plot for the
%   response Y and the predictor in column VNUM of X.  This plot illustrates
%   the incremental effect of this predictor in a regression model in which
%   the columns listed in INMODEL are used as predictors.  X is an N-by-P
%   matrix of predictor values.  Y is vector of N response values.  VNUM is a
%   scalar index specifying the column of X to use in the plot.  INMODEL is
%   either a vector of column numbers or a logical vector of P elements,
%   specifying the columns of X to use in the base model.  By default, all
%   elements of INMODEL are false (the model has no predictors).  ADDEDVARPLOT
%   automatically includes a constant term in the model.
%
%   ADDEDVARPLOT(X,Y,VNUM,INMODEL,STATS) uses the structure STATS containing
%   fitted model results created by the STEPWISEFIT function.  If STATS is
%   omitted, this function computes it.
%
%   An added variable plot contains data and fitted lines.  Suppose X1 is
%   column VNUM of X.  The data curve plots Y versus X1 after removing the
%   effects of the other predictors specified by INMODEL.  The solid line is
%   a least squares fit to the data curve, and its slope is the coefficient
%   that X1 would have if it were included in the model.  The dotted lines
%   are 95% confidence bounds for the fitted line, and they can be used to
%   judge the significance of X1.
%
%   If VNUM also appears in INMODEL, the plot that ADDEDVARPLOT produces is
%   sometimes known as a partial regression leverage plot.
%
%   Example:  Perform a stepwise regression on the Hald data, and create
%             an added variable plot for the predictor in column 2.
%      load hald
%      [b,se,p,in,stats] = stepwisefit(ingredients,heat);
%      addedvarplot(ingredients,heat,2,in,stats)

%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.3.4.2 $  $Date: 2005/04/18 22:20:25 $

error(nargchk(3,6,nargin,'struct'));

P = size(x, 2);

% Check for valid inputs
if numel(y)~=length(y)
   error('stats:addedvarplot:VectorRequired','Y must be a vector.');
end
if (nargin < 4)
    in = false(1,P);
elseif islogical(in)
   if length(in)~=P
      error('stats:addedvarplot:BadInModel',...
            'INMODEL must have one value for each column of X.');
   end
else
   if any(~ismember(in,1:P))
      error('stats:addedvarplot:BadInModel',...
            'INMODEL must be a list of X column numbers.');
   end
   in = ismember((1:P),in);
end

% Perform fit if fit results are not done; otherwise retrieve some results
if (nargin < 5)
    [B,SE,PVAL,in,stats] = stepwisefit(x,y, 'maxiter', 0, 'display', 'off',...
                                       'inmodel',in);
else
    if ~isequal(stats.source, 'stepwisefit')
      error('stats:addedvarplot:BadStats',...
            'STATS must be from the output of STEPWISEFIT.');
    end
    B = stats.B;
    SE = stats.SE;
end

% Argument 6 processed below

N = length(y);
bk = B(vnum);
se = SE(vnum);
alpha = 0.05;

if in(vnum)
   % Create a partial regression leverage plot for a column that is in
   r = y - x(:,in)*B(in);
   r = r - mean(r);
   xnotk = [ones(N,1) x(:,in & ((1:P)~=vnum))];
   xk = x(:,vnum);
   bnotk = xnotk \ xk;
   xr = xk - xnotk*bnotk;
   yr = B(vnum)*xr + r;
   ttl = sprintf('Partial regression leverage plot for X%d',vnum);
else
   % Created added-variable plot for an X column that is now out
   if stats.dfe==0
      error('stats:addedvarplot:NotEnoughData',...
            'Cannot do added variable plot with no error degrees of freedom.')
   end
   varlist = find(~in);
   outnum = find(varlist==vnum);
   xr = stats.xr(:,outnum);
   yr = stats.yr;
   ttl = sprintf('Added variable plot for X%d',vnum);
end

t = -tinv(alpha/2,stats.dfe);

% Restore NaNs to get back to original row numbers
[xr,yr] = statinsertnan(stats.wasnan,xr,yr);

% Create informative title
in2 = in;
in2(vnum) = 0;
runstart = find([in2(1), in2(2:end)&~in2(1:end-1)]);
runend   = find([in2(1:end-1)&~in2(2:end), in2(end)]);
txt = '';
for j=1:length(runstart)
   if runstart(j)==runend(j)
      txt = sprintf('%s,X%d',txt,runstart(j));
   elseif runstart(j)==runend(j)-1
      txt = sprintf('%s,X%d,X%d',txt,runstart(j),runend(j));
   else
      txt = sprintf('%s,X%d-X%d',txt,runstart(j),runend(j));
   end
end
if ~isempty(txt)     % add to title after removing extra comma
   ttl = sprintf('%s\nAdjusted for %s',ttl,txt(2:end));
end

% Create plot
xx = linspace(min(xr),max(xr))';
if (nargin < 6)
    f = gcf;
end

ax = get(f,'CurrentAxes');
if isempty(ax)
    ax = axes('Parent', f);
    set(ax,'Position',[.13 .11 .78 0.78]);
end
plot(xr,yr,'x',xx,bk*xx,'r-',...
     [xx;NaN;xx],[(bk+t*se)*xx;NaN;(bk-t*se)*xx],'r:', 'Parent', ax);
title(ttl,'Parent', ax);
xlabel(sprintf('X%d residuals',vnum), 'Parent', ax);
ylabel('Y residuals', 'Parent', ax);

legend(ax,'Adjusted data',sprintf('Fit: y=%g*x',bk),...
       sprintf('%g%% conf. bounds',100*(1-alpha)),0);

