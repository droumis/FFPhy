% R2ADJ:  Returns the adjusted squared multiple correlation, given SSR and SSTO, 
%         based on Wherry's (1931) correction.  Also returns the results of the 
%         F-test for overall significance of the regression model.
%           Input values may be scalars or column vectors or some combination of 
%         the two.  If at least one column vector is present, scalars are 
%         expanded to vectors and statistics are returned for each row of input.
%
%     Usage:  [R2a,F,pr,df] = r2adj(ssr,ssto,n,p)
%
%         ssr =   sum of squared deviations between observed and predicted values.
%         ssto =  sum of squared deviations from the mean for the dependent 
%                   variable.
%         n =     number of observations.
%         p =     number of parameters in model.
%         ------------------------------------------------------------------------
%         R2a =   adjusted R^2.
%         F =     F-statistic value.
%         pr =    significance level.
%         df =    numerator and denominator degrees of freedom.
%

% Note: Should implement the small-sample correction of Cattlin (1980)?

% Wherry RJ Sr. 1931. A new formula for predicting the shrinkage of the 
%   coefficient of multiple correlation.  Annals of Mathematical Statistics 
%   2:440-457.
% Catlin, P. 1980. Note on the estimation of the squared cross-validated 
%   multiple correlation of a regression model.  Psychological Bulletin 87:63-65.
% Neter J, W Wasserman, MH Kutner. 1985. Applied linear statistical models: 
%   regression, analysis of variance, and experimental designs.  2nd ed.  
%   Richard D. Irwin, Homewood IL.

% RE Strauss, 1/1/01
%   1/5/01 -    make F-test conditional, depending on 'nargout'.
%   1/6/01 -    added error messages.
%   11/9/01 -   corrected problem with ssto=0.
%   8/6/02 -    allow adjusted R2s to be negative.
%   11/6/02 -   add error checking for ssr==ssto and p>=n.
%   11/14/02 -  improve help documentation.

function [R2a,F,pr,df] = R2adj(ssr,ssto,n,p)
  get_F = 0;
  if (nargout > 1)
    get_F = 1;
  end;

  ssr = ssr(:);                             % Convert input to col vectors
  ssto = ssto(:);
  n = n(:);
  p = p(:);

  [ok,ssr,ssto,n,p] = samelength(ssr,ssto,n,p);  % Check for same length
  if (~ok)
    error('  R2ADJ: input matrices not consistent in length.');
  end;

  R2a = NaN;
  F = NaN;
  pr = NaN;
  df = NaN;
  
  if (any(p>=n))
    disp('  R2ADJ warning: number of parameters must be less than sample size.');
    return;
  end;
  if (any(ssr>=ssto))
    disp('  R2ADJ warning: SSR must be less than SSTO');
    return;
  end;
  if (any(ssto<eps))
    i = find(ssto<eps);
    ssto(i) = NaN*ones(length(i),1);
  end;
  
  sse = ssto-ssr;
  R2a = 1 - (((n-1)./(n-p)).*(sse./ssto));

  if (get_F)                                % F test
    df1 = max([ones(length(p),1), p'-1])';
    df2 = n-p;
    df = [df1 df2];
    F = (ssr./df1)./(sse./df2);
    pr = 1-fcdf(F,df1,df2);
  end;

  return;
