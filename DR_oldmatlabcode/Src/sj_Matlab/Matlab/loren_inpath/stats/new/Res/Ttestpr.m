% TTESTPR:  Paired t-test of T-C, for matched treatment and control groups, 
%           using a specified combination of T and C as a potential covariate to 
%           adjust for magnitude.  If the correlation between (T-C) and the 
%           covariate is significant, (T-C) is regressed on the covariate and the 
%           residuals are added to the predicted value of (T-C) corresponding to 
%           the mean covariate for the t-test.
%
%     Usage: [meandiff,pr,t,df,r,rpr,b] = ...
%                       ttestpr(T,C,{dir},{wts|covar},{r_alpha},{doplot})
%
%           T =         vector of values for the control group.
%           C =         matching vector of values for the treatment group.
%           dir =       value (pos, zero, neg) indicating direction of test:
%                           + = right-tailed test for T > C;
%                           0 = two-tailed test for T = C [default];
%                           - = left-tailed test for T < C.
%           wts =       vector of two weights indicating the linear combination 
%                         of T & C, respectively, to be used as the covariate.  
%                         If both weights are zero, the covariate-correction 
%                         procedure is not applied.  Examples:
%                           [1,0] = treatment variable (T) used as covariate;
%                           [0,1] = control variable (C) used as covariate [default];
%                           [1,1] = sum of treatment and control (T+C) used as 
%                               covariate;
%                           [0.5,0.5] = mean of treatment and control ((T+C)/2) 
%                               used as covariate;
%                           [0,0] = no covariate.
%           covar =     vector (same length as T & C) of values for an external 
%                         covariate (used in place of wts).
%           r_alpha =   alpha-level used to judge the correlation to be 
%                         significant; set to 1 for the correction always to be 
%                         supplied [default = 1].
%           doplot =    boolean flag indicating, if true, that a scatter of 
%                         (T-C) vs the covariate is to be plotted, with 
%                         regression line [default = 0].
%           ----------------------------------------------------------------------
%           meandiff =  mean difference, T-C.
%           pr =        significance level for the test.
%           t =         t-value for the test.
%           df =        degrees of freedom for the test.
%           r =         correlation between (T-C) and C.
%           rpr =       significance of correlation between (T-C) and C.
%           b =         intercept and slope of regression of (T-C) on C.
%

function [meandiff,pr,t,df,r,rpr,b] = ttestpr(T,C,dir,wts,r_alpha,doplot)
  if (nargin < 3) dir = []; end;
  if (nargin < 4) wts = []; end;
  if (nargin < 5) r_alpha = []; end;
  if (nargin < 6) doplot = []; end;

  if (~isvector(T) | ~isvector(C))
    error('  TTESTPR: treatment and control data must be vectors.');
  end;
  if (length(T) ~= length(C))
    error('  TTESTPR: treatment and control vectors not of same length.');
  end;
  T = T(:);
  C = C(:);

  if (isempty(dir))
    dir = 0;
  end;
  if (isempty(wts))
    wts = [0 1];
  end;
  if (isempty(r_alpha))
    r_alpha = 0.05;
  end;
  if (r_alpha > 1)
    r_alpha = r_alpha/100;
  end;

  if (length(wts)~=2 & length(wts)~=length(T))
    error('  TTESTPR: covariate/weights vector not of proper length.');
  end;

  use_covar = 0;
  if (any(abs(wts)>eps))
    use_covar = 1;
  end;

  diff = T-C;                             % Get paired differences
  origdiff = diff;                        %   and save them

  if (use_covar)                          % Get resids from regr on covariate
    if (length(wts) > 2)                    % Get covariate
      covar = wts;
    else
      covar = [T C]*wts(:);
    end;

    [r,rpr] = corr(diff,covar);                 % Corr between diff & control
    if (rpr<=r_alpha | doplot)              % Regress out covariate, use resids
      cmin = min(covar);
      cmax = max(covar);
      [b,stats,pred,resid] = linregr(covar,diff,[cmin;mean(covar);cmax]);
      if (rpr<=r_alpha)
        diff = pred(2)+resid;
      end;
    end;
  end;

  n = length(diff);
  df = n-1;
  meandiff = mean(diff);
  t = meandiff/sqrt(var(diff)/n);

  if (dir>eps)                            % Right-tailed
    pr = 1-tcdf(t,df);
  elseif (dir<-eps)                       % Left-tailed
    pr = tcdf(t,df);
  else                                    % Two-tailed
    pr = 2*(1-tcdf(abs(t),df));
  end;

  if (doplot & use_covar)
    plot(covar,origdiff,'ko');
    hold on;
    plot([cmin cmax],[pred([1,3])],'k');
    plot([cmin cmax],[0 0],'k:')
    hold off;
    putbnd([covar;cmin;cmax],[origdiff;pred([1,3])]);
    putxlab('Covariate');
    putylab('T-C');
  end;

  return;
