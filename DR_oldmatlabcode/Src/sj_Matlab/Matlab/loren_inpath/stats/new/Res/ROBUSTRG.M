% ROBUSTRG: Robust predictive multivariate linear regression of Y on X,
%           based on iterative Welsch weighting to downweight outliers.
%
%     Syntax: [b,r2,s2,pred,resid,d,sr] = robustrg(X,Y,{Xp},{noint})
%
%        X =     [n x p] matrix of independent variables.
%        Y =     [n x 1] vector of dependent variable.
%        Xp =    [m x p] optional matrix of independent variables for which
%                 predicted values are to be estimated.
%        noint = 1 (=TRUE) for regression thru the origin, 0 (=FALSE) otherwise
%                 [default 0].
%        ----------------------------------------------------------------------
%        b =     [p+1 x q] matrix of regression coefficients; the first row
%                 gives the intercepts, the second row gives the
%                 coefficients for X1, etc.  The first row is omitted for
%                 regression thru the origin.
%        r2 =    [1 x q] vector of coefficients of determination.
%        s2 =    [1 x q] vector of MSEs (residual variances).
%        pred =  [n|m x q] matrix of predicted values.  If a non-null input
%                 matrix Xp is provided, predicted values are estimated from it;
%                 otherwise predicted values are estimated from X.
%        resid = [n x q] matrix of residuals.
%        d =     [n x q] matrix of Cook's-D influence values, which measure
%                 the influence of observations on the regression.
%        sr =    [n x q] matrix of studentized residuals, distributed as t
%                 with n-3 df.  When testing the hypothesis that all n obs
%                 are not outliers, the Bonferroni-adjusted type I error
%                 is >= n*alpha.

% Staudte, R.G. and S.J. Sheather. 1990. Robust Estimation and Testing.  
%   Wiley. Interscience.
% Jobson, J.D. 1991. Applied Multivariate Data Analysis. Vol.1: Regression
%   and Experimental Design.  Springer-Verlag.

% RE Strauss, 11/5/97
%   9/20/99 - update handling of null input arguments.

function [b,r2,s2,pred,resid,d,sr] = robustrg(X,Y,Xp,noint)
  if (nargin < 3) Xp = []; end;
  if (nargin < 4) noint = []; end;

  [n,p] = size(X);
  [ny,q] = size(Y);

  given_Xp = 0;
  if (~isempty(Xp))
    given_Xp = 1;
  end;
  if (isempty(noint))
    noint = 0;
  end;

  calc_r2 = 0;                            % Output flags
  calc_s2 = 0;
  calc_pred = 0;
  calc_resid = 0;
  calc_d = 0;

  if (nargout >= 2)
    calc_r2 = 1;
    r2 = zeros(1,q);
  end;
  if (nargout >= 3)
    calc_s2 = 1;
    s2 = zeros(1,q);
  end;
  if (nargout >= 4)
    calc_pred = 1;
  end;
  if (nargout >= 5)
    calc_resid = 1;
    resid = zeros(n,q);
  end;
  if (nargout >= 6)
    calc_d = 1;
    d = resid;
  end;

  if (n==1)                       % Transpose row vectors into col vectors
    X = X';
    [n,p] = size(X);
  end;
  if (ny==1)
    Y = Y';
    [ny,q] = size(Y);
  end;
  if (q>1)
    error('ROBUSTRG: dependent variable must be a vector');
  end;
  if (n ~= ny)
    error('ROBUSTRG: X,Y must have same number of observations');
  end;

  if (~noint)
    X = [ones(n,1) X];                % Augment the X-matrix
  end;

  H = X*inv(X'*X)*X';                 % Hat matrix
  h = diag(H);                        % Residual leverage values (h-hat)
  Yhat = H*Y;                         % Predicted values
  e = Y - Yhat;                       % Residuals
  s2 = e'*e./(n-p-1);                 % MSE
  b0 = inv(X'*X)*X'*Y;                % Initial regression parameters
b0
b00 = b0;

  s2i = zeros(n,1);                   % Studentized residuals
  for i=1:n                           % Remove obs i from data matrices
    Xi = X;                        
    Yi = Y;
    Xi(i,:) = [];
    Yi(i) = [];
    bi = inv(Xi'*Xi)*Xi'*Yi;            % Regression with omitted obs
    ei = Yi - Xi*bi;                    % Residuals from regr with omitted obs
    s2i(i) = ei'*ei / (n-p-2);          % Corrected residual var for omitted obs
  end;
  sr = e ./ (sqrt(s2i).*sqrt(1-h));   % Studentized residuals
  dfits = sqrt(h./(1-h)).*sr;         % DFITS

  k = 2;                              % Welsch constants
  c = k*sqrt((p+1)/n);

  w = min(c./abs(dfits),1);           % Initial weights
  W = diag(w);                        % Diag matrix of weights

  epsilon = 1e-4*ones(size(b0));      % Accuracy of convergence
  b = inv(X'*W*X)*X'*W*Y;
b

  while (any(abs(b-b0) > epsilon))         % Iterate to convergence
    b0 = b;

    for i=1:n                           % Remove obs i from data matrices
      Xi = X;                        
      Yi = Y;
      Wi = W;
      Xi(i,:) = [];
      Yi(i) = [];
      Wi(i,:) = [];
      Wi(:,i) = [];

      bi = inv(Xi'*Wi*Xi)*Xi'*Wi*Yi;        % Weighted regr with omitted obs
      ei = Yi - Xi*bi;                    % Residuals from regr with omitted obs
      s2i(i) = ei'*ei / (n-p-2);          % Corrected residual var for omitted obs
    end;
    sr = e ./ (sqrt(s2i).*sqrt(1-h));   % Studentized residuals
    dfits = sqrt(h./(1-h)).*sr;         % DFITS
    w = min(c./abs(dfits),1);           % Initial weights
    W = diag(w);                        % Diag matrix of weights
    b = inv(X'*W*X)*X'*W*Y;             % New regr coeffs
b
  end;

  Yhat = X * b;                       % Predicted values
  e = Y - Yhat;                       % Residuals

  if (calc_r2)                        % Coeff of multiple determination
    r2(c) = 1 - (var(e)/var(Y));
  end;

  if (calc_s2)                        % MSE = residual variance
    s2(c) = e'*e / (n-p-1);
  end;

  if (calc_resid)                   % Save residuals
    resid(:,c) = e;
  end;

  if (calc_d)                       % Cook's D measure of obs influence
    d(:,c) = (h./(1-h)).*((e.^2)./((p+1)*s2(c).*(1-h))); % Calc D values
  end;

  if (calc_pred)                      % Calc predicted values for original X
    if (given_Xp)                       % Or for substitute matrix, if given one
      X = Xp;
      [nx,p] = size(X);
      if (~noint)
        X = [ones(nx,1) X];           % Augment
      end;
    end;
    pred = X*b;                        % Predicted values
  end;

close all;
Xp = [min(X);max(X)];
pred0 = Xp*b00;
predf = Xp*b;
Xp = Xp(:,2);
X = X(:,2);
plot(X,Y,'yo',Xp,pred0,'r',Xp,predf,'y');

  return;
