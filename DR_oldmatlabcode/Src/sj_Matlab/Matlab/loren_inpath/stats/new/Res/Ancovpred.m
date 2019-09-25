% ANCOVPRED:  Find predicted values of given points at the grand-mean abscissa 
%             (or other specified value), assuming a common slope (either 
%             among-group or pooled within-group) for multiple groups.  Uses the 
%             among-group slope if the pooled within-group slope is not 
%             statistically significant.
%
%     Usage: [ygpred,yipred,b1,b0,predval,useag] = ...
%                         ancovpred(x,y,grps,{getse},{predval},{doplot},{alpha})
%
%             x,y =     [n x 1] corresponding vectors of independent- and 
%                         dependent-variable values.
%             grps =    [n x 1] vector of group identifiers.
%             getse =   optional boolean flag indicating, if true, that standard 
%                         errors of estimates are to be provided [default = 1].
%             predval = optional abscissa value at which DV values are to 
%                         predicted [default = grand mean X].
%             doplot =  optional boolean flag indicating, if true, that
%                         scatterplots are to be produced [default = 0].
%             alpha =   optional alpha-level for determining significance of 
%                         pooled within-group slope [default = 0.05].
%             -------------------------------------------------------------------
%             ygpred =  [k x 1] predicted values at x=predval for group means; or
%                       [k x 2] matrix of predicted values and corresponding 
%                         standard errors.
%             yipred =  [n x 1] predicted values of individual observations at 
%                         x=predval + the within-group residual; or 
%                       [n x 2] matrix of predicted values and corresponding 
%                         standard errors.
%             b1 =      common regression slope; or [1 x 2] vector of slope and 
%                         standard error.
%             b0 =      [k x 1] regression intercepts; or [k x 2] matrix of 
%                         intercepts and standard errors.
%             predval = value at which DV values were predicted.
%             useag =   boolean flag indicating, if true, that the among-group 
%                         slope was used rather than the pooled within-group slope.
%

% RE Strauss, 10/31/00
%   4/22/01 - produce two plots rather than one.

function [ygpred,yipred,b1,b0,predval,useag] = ...
                                ancovpred(x,y,grps,getse,predval,doplot,alpha)
  if (nargin < 3)
    error('  ANCOVPRED: expected input values not passed.');
  end;

  if (nargin < 4) getse = []; end;
  if (nargin < 5) predval = []; end;
  if (nargin < 6) doplot = []; end;
  if (nargin < 7) alpha = []; end;

  err = 0;
  if (~isvector(x) | ~isvector(y) | ~isvector(grps))
    err = 1;
  end;

  x = x(:);
  y = y(:);
  grps = grps(:);
  n = length(x);
  
  if (length(y)~=n | length(grps)~=n | err)
    error('  ANCOVPRED: data and group identifiers must be vectors of equal length.');
  end;

  if (isempty(getse))
    getse = 1;
  end;
  if (isempty(predval))
    predval = mean(x);
  end;
  if (isempty(doplot))
    doplot = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;
  if (alpha > 1)
    alpha = alpha/100;
  end;

  [ug,fg] = uniquef(grps);              % Group identifiers and freqs
  ngrps = length(ug);                   % Number of groups

  if (min(fg)<3 & getse)
    disp('  ANCOVPRED warning: need at least 3 obs/grp for standard errors.');
  end;

  if (getse)                            % Allocate output matrices
    b0 = zeros(ngrps,2);
    ygpred = zeros(ngrps,2);
    yipred = zeros(n,2);
  else
    b0 = zeros(ngrps,1);
    ygpred = zeros(ngrps,1);
    yipred = zeros(n,1);
  end;

  xc = grpcentr(x,grps);                % Center groups on x & y
  yc = grpcentr(y,grps);

  useag = 0;
  [b,stats] = linregr(xc,yc);           % Pooled within-grp regression
  p = stats(6);
  if (p > alpha)                        % If not significant,
    [b,stats] = linregr(x,y);           %   use among-group regression
    useag = 1;
  end;

  b1 = b(2);                            % Stash slope
  if (getse)                            %   and standard error            
    mse = stats(2);
    seb1 = sqrt(mse/(x'*x));
  end;

  for g = 1:ngrps                       % Cycle thru groups
    i = find(grps==ug(g));                % Get obs in current grp
    ni = length(i);                       % Sample size for current grp
    rni = 1./ni;                          % Reciprocal of sample size

    ymean = mean(y(i));                   % Within-grp means
    xmean = mean(x(i));
    ssx = sum((x(i)-xmean).^2);           % Within-grp ssq of x

    b0i = ymean - b1*xmean;               % Get and stash intercept

    ei = y(i) - (b1*x(i)+b0i);            % Within-grp residuals
    if (ni>=3)                            % Within-grp mse
      msei = (ei'*ei)/(ni-2);                
    else
      msei = NaN;
    end;

    gpred = b1*predval + b0i;             % Predicted y at x=predval

    if (getse)                            % Stash output
      seb0i = sqrt(msei*(rni + (xmean^2)/ssx));
      b0(g,:) = [b0i, seb0i];

      sepredi = sqrt(msei*(rni + (predval-xmean).^2/ssx));
      ygpred(g,:) = [gpred, sepredi];
      yipred(i,:) = [gpred+ei, sepredi*ones(ni,1)];
    else
      b0(g) = b0i;
      ygpred(g) = gpred;
      yipred(i) = gpred+ei;
    end;

  end;

  if (getse)                            % Stash standard error of pooled slope        
    b1 = [b1 seb1];
  end;

  if (doplot)                           % Scatterplot of groups and regressions
    figure;
    plotgrps(x,y,grps,tostr(uniquef(grps,1))); % Scatterplot by group

    figure;
    plotgrps(x,y,grps,tostr(uniquef(grps,1)),[],[0 0 0 1]); % Scatterplot by group
    xx = x;
    yy = y;

    hold on;
    for g = 1:ngrps                       % Cycle thru groups
      i = find(grps==ug(g));                % Get obs in current grp

      buffer = 0.05*range(x(i));            % Bounds on x
      xmin = min(x(i)) - buffer;    
      xmin = min([xmin, predval]);        

      xmax = max(x(i)) + buffer;
      xmax = max([xmax, predval]);

      yxmin = b1(1)*xmin + b0(g,1);         % Corresponding y values
      yxmax = b1(1)*xmax + b0(g,1);
      plot([xmin xmax],[yxmin yxmax],'k');  % Plot regression line

      xx = [xx; xmin; xmax];
      yy = [yy; yxmin; yxmax];
    end;

    buffer = 0.03*range(yy);                % Vertical dotted line for predval
    plot([predval predval],[min(yy)-buffer max(yy)+buffer],'k:'); 
    putbnd(xx,yy);
    hold off;
  end;

  return;
