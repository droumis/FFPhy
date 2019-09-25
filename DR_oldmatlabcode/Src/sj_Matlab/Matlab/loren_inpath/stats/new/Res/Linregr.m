% LINREGR: Predictive multivariate linear regression of Y on X, optionally 
%          weighted.  Does not allow for missing data.
%
%     Syntax: [b,stats,pred,resid,stderrs,B,jpred,jresid,d,Fd,sr] = linregr(X,Y,{w},{Xp},{Yp})
%
%        X =       [n x p] matrix of independent variables.
%        Y =       [n x q] matrix of dependent variables, regressed separately
%                    by column.
%        w =       optional [n x 1] vector of weights, for weighted regression.
%        Xp =      optional [m x p] matrix of independent variables for auxiliary 
%                    data for which predicted values are to be estimated.
%        Yp =      optional [m x q] matching matrix of dependent variables for 
%                    auxiliary data, from which to calculate residuals.
%        -------------------------------------------------------------------------
%        b =       [p+1 x q] matrix of regression coefficients; the first row
%                    gives the intercepts (b0), the second row gives the
%                    coefficients for X1 (b1), etc.  
%        stats =   [6 x q] matrix of regression statistics; each column gives the
%                    intercepts for a single dependent variable:
%                      row 1: adjusted coefficient of determination (R_a^2);
%                          2: MSE (residual variance) or PMSE (predicted
%                               residual variance from jackknifed residuals,
%                               if requested;
%                          3: F-statistic value;
%                          4: df1 (numerator degrees of freedom, =dfr);
%                          5: df2 (denominator degrees of freedom = dfe);
%                          6: pr >= F (probability under the null hypothesis).
%        pred =    [n|m x q] matrix of predicted values.  If a non-null input
%                    matrix 'Xp' is provided, predicted values are estimated from it;
%                    otherwise predicted values are estimated from X.
%        resid =   [n|m x q] matrix of residuals (Y-Yhat).  If 'Xp' and 'Yp' are 
%                    provided, the [m x q] matrix of predicted residuals is returned 
%                    instead.
%        stderrs = matching matrix of standard errors of regression coefficients.
%        B =       matching matrix of stardized (beta) coefficients.
%        jpred =   [n x q] matrix of jackknifed predicted values.
%        jresid =  [n x q] matrix of jackknifed residuals.
%        d =       [n x q] matrix of Cook's-D influence values, which measure
%                    the influence of observations on the regression.
%        Fd =      [n x q] matrix of F percentiles corresponding to Cook's D values;
%                    if the percentiles are 50% or more, the observations can be 
%                    considered to have a considerable influence.
%        sr =      [n x q] matrix of jackknifed studentized residuals, distributed as 
%                    t with n-3 df.  When testing the hypothesis that all n obs
%                    are not outliers, the Bonferroni-adjusted type I error
%                    is >= n*alpha.  If 'Xp' and 'Yp' are provided, the [m x q] 
%                    matrix of these (non-jackknifed) residuals is returned instead,
%                    but based on the jackknifed residual variances.
%

% Neter, J., W. Wasserman, M.H. Kutner.  1985.  Applied Linear Statistical Models: 
%   Regression, Analysis of Variance, and Experimental Designs (2nd ed.).
%   Richard D. Irwin, Homewood IL.
% Jobson, J.D. 1991. Applied Multivariate Data Analysis. Vol.1: Regression
%   and Experimental Design.  Springer-Verlag.

% RE Strauss, 5/9/95
%   4/30/98 -   added F values, df's, and probabilities; organized into single
%                 'stats' matrix along with r2 and s2.
%   9/3/99 -    changed usage of I/O flags.
%   2/24/00 -   if Xp is row vector, transpose to column vector.
%   10/26/00 -  check for zero or negative dfe;
%                 if Xp is a row vector, transpose to column vector only if 
%                 length not equal to number of variables.
%   1/1/01 -    corrected degrees of freedom; 
%                 call r2adj() to calc coeff of multiple correlation.
%   2/2/01 -    added option for weighted regression.
%   4/13/01 -   added variance among jackknifed predicted residuals.
%   7/27/01 -   corrected degrees of freedom.
%   7/29/01 -   removed no-intercept option.
%   9/27/01 -   produce error message for missing data.
%   12/7/01 -   create weighting matrix only for weighted regression.
%   12/13/01 -  added F-percentiles for Cook's D values.
%   12/20/01 -  isolate linregrf(); calculate 'resid', 'sr' and 'jrv' for Xp 
%                 if provided.
%   3/2/02 -    corrected number of parameters for regression model.
%   9/12/02 -   corrected problem with scale of studentized residuals when Xp,Yp provided.
%   9/17/02 -   added calculation of jackknifed residuals (unadjusted); removed output 'jvar'.
%   9/23/02 -   added calculation of jackknifed predicted values.
%   9/29/02 -   added estimation of standard errors and standardized coefficients.
%   10/1/02 -   removed 'calc_stats' flag in pass to linregrf().
%   10/8/02 -   changed sequence of output arguments to be consistent with previous versions.
%   11/14/02 -  removed vector check for Xp.
%   12/6/02 -   altered vector check for Xp based on number of parameters of model.

function [b,stats,pred,resid,stderrs,B,jpred,jresid,d,Fd,sr] = linregr(X,Y,w,Xp,Yp)
  if (~nargin) help linregr; return; end;
  
  if (nargin < 3) w = []; end;
  if (nargin < 4) Xp = []; end;
  if (nargin < 5) Yp = []; end;

  nstats = 6;                     % Number of regression statistics calculated
  
  if (isvector(X))
    X = X(:);
  end;
  if (isvector(Y))
    Y = Y(:);
  end;

  [n,p] = size(X);
  [ny,q] = size(Y);

  if (any(~isfinite(X(:))) | any(~isfinite(Y(:))))
    error('  LINREGR: input matrices contain missing data.');
  end;

  if (n ~= ny)
    error('  LINREGR: X,Y must have same number of observations');
  end;

  X = [ones(n,1) X];                    % Augment the X-matrix
  b = zeros(p+1,q);                     % Allocate parameter matrices
  stderrs = zeros(p+1,q);
  B = zeros(p+1,q);

  if (isempty(w))                       % Weight vector supplied?
    weighted = 0;
    W = [];
  else
    weighted = 1;
    W = diag(w,0);                      % Diagonal weight matrix
  end;

  given_Xp = 0;                         % Prediction matrices supplied?
  given_Yp = 0;
  if (~isempty(Xp))
    given_Xp = 1;
    if (isvector(Xp))
      Xp = Xp(:);
      if (length(Xp)==p)
        Xp = Xp';
      end;
    end;
    [m,pp] = size(Xp);

    if (pp~=p)
      error('  LINREGR: matrix of IVs for auxiliary has incorrect number of columns');
    end;
    Xp = [ones(m,1) Xp];                % Augment
  end;
  if (~isempty(Yp))
    if (~given_Xp)
      error('  LINREGR: Y-prediction vector must have matching X-prediction matrix.');
    end;
    if (size(Yp,1)~=m)
      error('  LINREGR: X- & Y-prediction matrices must have same number of observations.');
    end;
    if (size(Yp,2)~=q)
      error('  LINREGR: Y-prediction matrix must have same number of columns as Y.');
    end;
    given_Yp = 1;
  end;

  calc_stats = 0;                       % Set processing flags
  calc_pred = 0;
  calc_resid = 0;
  calc_jpred = 0;
  calc_jresid = 0;
  calc_d = 0;
  calc_sr = 0;
  
  if (nargout >= 2)
    calc_stats = 1;
    stats = zeros(nstats,q);
  end;
  if (nargout >= 3)
    calc_pred = 1;
  end;
  if (nargout >= 4)
    calc_resid = 1;
  end;
  if (nargout >= 7)
    calc_jpred = 1;
    jpred = zeros(n,q);                   
  end;
  if (nargout >= 8)
    calc_jresid = 1;
    jresid = zeros(n,q);
    s2i = zeros(n,1);
  end;
  if (nargout >= 9)
    calc_d = 1;
    d = zeros(n,q);
    Fd = zeros(n,q);
  end;
  if (nargout >= 11)
    calc_sr = 1;
    if (given_Yp) 
      sr = zeros(m,q);
    else
      sr = zeros(n,q);
    end;
  end;
  
  if (calc_d)                           % Residual leverage values (h-hat) for Cook's D
    if (weighted)
      h = diag(X*inv(X'*W*X)*X');         
    else
      h = diag(X*inv(X'*X)*X');
    end;
  end;

  for c = 1:q                           % Regress each column of Y on X
    y = Y(:,c);                           % Current Y vector
    [b_fit,sb,beta,st,e] = linregrf(X,y,W,weighted);
    b(:,c) = b_fit;
    stderrs(:,c) = sb;
    B(:,c) = beta;
    
    if (calc_stats)                       % Save stats
      stats(:,c) = st;
    end;

    if (calc_resid & ~given_Yp)           % Save residuals
      resid(:,c) = e;
    end;

    if (calc_d)                           % Cook's D measure of obs influence
      mse = st(2);
      d(:,c) = (h./(1-h)).*((e.^2)./((p+1)*mse.*(1-h))); % Calc D values
      Fd(:,c) = fcdf(d(:,c),n,p);           % F-percentiles for D values
    end;

    if (calc_jpred | calc_jresid)             % Jackknifed predicted and residual values
      for i=1:n
        Xi = X;                               % Remove obs i from data matrices
        Xi(i,:) = [];
        yi = y;
        yi(i) = [];
        if (weighted)                         % Regression with omitted obs
          Wi = W;
          Wi(:,i) = [];
          Wi(i,:) = [];
        else
          Wi = [];
        end;
        bi = linregrf(Xi,yi,Wi,weighted);
        predi = X*bi;                         % Predicted values with omitted obs
        ei = y - predi;                       % Residuals from regr with omitted obs
        s2i(i) = ei'*ei / max([(n-p-2)],1);   % Corrected residual var for omitted obs
        jpred(i,c) = predi(i);
        jresid(i,c) = ei(i);
      end;       
      if (calc_sr & ~given_Yp)
        sr(:,c) = jresid(:,c) ./ sqrt(s2i);   % Studentized jackknifed residuals
      end;
    end;
  end;

  if (calc_pred)                          % Calc predicted values for original X
    if (given_Xp)                           % Or for substitute matrix, if provided
      pred = Xp*b;                            % Predicted values
    else
      pred = X*b;
    end;
  end;
  
  if (given_Yp)                           % If given auxiliary data,
    if (calc_resid)
      resid = Yp - pred;                    % Save residuals
    end;
    if (calc_sr)                            % Standardized predicted residuals based on
      sr = resid ./ (ones(size(resid,1),1)*std(jresid)); % jackknifed residual variances
    end;
  end;

  return;
