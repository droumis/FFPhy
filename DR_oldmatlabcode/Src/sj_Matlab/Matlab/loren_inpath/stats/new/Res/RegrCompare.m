% RegrCompare: Randomized permutation comparison of two or more linear regression 
%              lines for homogeneity.  The improvement-of-fit statistics test for
%              the improvement of fit of separate regressions with respect to the
%              pooled regression (assuming that all groups are identical).
%              If the response variable is censored (at the same level for all groups),
%              uses Tobit regression to estimate parameters and statistics.
%              Removes missing values in the predictor variable.
%
%     Usage: [b,iof_stats,slope_stats,intcpt_stats,limit] = ...
%                                       regrcompare(x,y,g,{se_iter},{pr_iter},{limit})
%
%         x =       [N x 1] vector for independent variable.
%         y =       matching vector for dependent variable.
%         g =       matching group-identification vector for k groups.
%         se_iter = number of iterations to estimate standard errors of regression
%                     coefficients for Tobit regression; required in the case of
%                     censored data, otherwise ignored.
%         pr_iter = optional number of random permutation iterations for regression-
%                     comparison statistics.
%         limit =   optional censoring limit value for dependent variable 
%                     [default = minimum non-censored value - eps].
%         --------------------------------------------------------------------------
%         b = [2 x k] matrix of regression coefficients for each of k groups.
%                The first row contains the intercepts and the second the slopes.
%         iof_stats = vector of improvement-of-fit statistics for regression lines:
%                F       - F statistic;
%                pr_asym - asymptotic significant level;
%                pr_boot - bootstrapped significant level;
%                sse(r)  - sum-of-squared error for reduced (pooled-regr) model;
%                sse(f)  - sse for full (separate-regr) model;
%                df(r)   - degrees of freedom for reduced model;
%                df(f)   - degrees of freedom for full model.
%         slope_stats = vector of anova statistics for heterogeneity of slopes:
%                F       - F statistic;
%                pr_asym - asymptotic significance level;
%                pr_boot - optional randomized significance level;
%                df1     - numerator degrees of freedom;
%                df2     - denominator degrees of freedom.
%         intcpt_stats = vector of anova statistics for heterogeneity of intercepts:
%                F       - F statistic;
%                pr_asym - asymptotic significance level;
%                pr_boot - optional randomized significance level;
%                df1     - numerator degrees of freedom;
%                df2     - denominator degrees of freedom.
%         limit = censoring limit value used in analysis.
%               

% RE Strauss, 10/14/02
%   3/19/03 - implements Tobit regression if y is censored.
%   8/11/03 - add 'limit' to output arguments.

function [b,iof_stats,slope_stats,intcpt_stats,limit] = regrcompare(x,y,g,se_iter,pr_iter,limit)
  if (nargin < 4) se_iter = []; end;
  if (nargin < 5) pr_iter = []; end;
  if (nargin < 6) limit = []; end;
  
  if (isempty(se_iter)) se_iter = 0; end;
  if (isempty(pr_iter)) pr_iter = 0; end;
  
  i = find(~isfinite(x));                     % Remove missing values in predictor
  if (~isempty(i))
    x(i) = [];
    y(i) = [];
    g(i) = [];
  end;
  
  do_tobit = 0;
  i = find(~isfinite(y));                     % If have missing data in response variable, 
  if (~isempty(i))                            %   try Tobit regression
    do_tobit = 1;
    limit = min(y(isfinite(y))) - eps;
    if (se_iter==0)
      error('  RegrCompare: Tobit regression standard errors must be bootstrapped (se_iter>0)');
    end;
  end;
  
  [grp_id,n] = uniquef(g);                    % Groups
  ngrps = length(grp_id);
  N = length(x);                              % Sample size
  
  if (do_tobit)
    [b_pooled,stats_pooled] = censoredregr(x,y,limit,[],0);  % Pooled regression
    [b_sep,stats_sep] = censoredregr([x,g],y,limit,[],0);    % Separate regressions    
  else
    [b_pooled,stats_pooled] = linregr(x,y);     % Pooled regression
    [b_sep,stats_sep] = linregr([x,g],y);       % Separate regressions
  end;
  
  b = zeros(2,ngrps);
  s = zeros(2,ngrps);
  ssef = 0;
  for i = 1:ngrps                             % Regressions by group
    ig = find(g==grp_id(i));
    if (do_tobit)
      [b(:,i),stats,s(:,i)] = censoredregr(x(ig),y(ig),limit,se_iter,0);     
    else
      [b(:,i),stats,p,r,s(:,i)] = linregr(x(ig),y(ig));      
    end;
    s(:,i) = s(:,i)*sqrt(n(i));
    ssef = ssef + stats(2)*stats(5);
  end;
  
  dfr = N-2;                                           % Improvement-of-fit statistics
  dff = N-2*ngrps;
  sser = stats_pooled(2)*stats_pooled(5);     
  F_iof = (((sser-ssef)-(dfr-dff))/(dfr-dff))/(ssef/dff);
  F_iof = max([0,F_iof]);
  pr_asym = 1-fcdf(F_iof,dfr-dff,dff);
  pr_boot = NaN;
  iof_stats = [F_iof,pr_asym,pr_boot,sser,ssef,dfr,dff]';
  
  [F_slopes,df,pr_asym] = anovaparam(b(2,:),s(2,:),n); % Slope statistics
  pr_boot = NaN;
  slope_stats = [F_slopes,pr_asym,pr_boot,df(1),df(2)]';
  
  [F_intcpt,df,pr_asym] = anovaparam(b(1,:),s(1,:),n); % Intercept statistics
  pr_boot = NaN;
  intcpt_stats = [F_intcpt,pr_asym,pr_boot,df(1),df(2)]';
  
  if (pr_iter)                                % Random permutations 
    F_iof_null = zeros(pr_iter,1);
    F_slopes_null = zeros(pr_iter,1);
    F_intcpt_null = zeros(pr_iter,1);
    
    for it = 1:pr_iter
      g = g(randperm(length(g)));               % Randomly permute group membership vector
  
      ssef = 0;
      for i = 1:ngrps                           % Regressions by group
        ig = find(g==grp_id(i));
        if (do_tobit)
          [b(:,i),stats,s(:,i)] = censoredregr(x(ig),y(ig),limit,se_iter,0);
        else
          [b(:,i),stats,p,r,s(:,i)] = linregr(x(ig),y(ig));
        end;
        s(:,i) = s(:,i)*sqrt(n(i));
        ssef = ssef + stats(2)*stats(5);
      end;
      
      F_iof_null(it) = (((sser-ssef)-(dfr-dff))/(dfr-dff))/(ssef/dff);
      F_iof_null(it) = max([0,F_iof_null(it)]);
      F_slopes_null(it) = anovaparam(b(2,:),s(2,:),n); % Slope statistics
      F_intcpt_null(it) = anovaparam(b(1,:),s(1,:),n); % Intercept statistics
    end;
    iof_stats(3) = randprob(F_iof,sort(F_iof_null));
    slope_stats(3) = randprob(F_slopes,sort(F_slopes_null));
    intcpt_stats(3) = randprob(F_intcpt,sort(F_intcpt_null));
  end;
  
  return;
  