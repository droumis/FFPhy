% LTSREGR: Least trimmed squares (LTS) robust linear regression of Y on X.  
%          Minimizes the sum of the n-n/2 largest squared errors.
%
%     Syntax: [b,pred,resid] = ltsregr(X,Y,{Xp},{noint})
%
%        X =     [n x p] matrix of independent variables.
%        Y =     [n x q] matrix of dependent variables, regressed separately
%                 by column.
%        Xp =    [m x p] optional matrix of independent variables for which
%                 predicted values are to be estimated.
%        noint = 1 (=TRUE) for regression thru the origin, 0 (=FALSE) otherwise
%                 [default 0].
%        ----------------------------------------------------------------------
%        b =     [p+1 x q] matrix of regression coefficients; the first row
%                 gives the intercepts, the second row gives the
%                 coefficients for X1, etc.  The first row is omitted for
%                 regression thru the origin.
%        pred =  [n|m x q] matrix of predicted values.  If a non-null input
%                 matrix Xp is provided, predicted values are estimated from it;
%                 otherwise predicted values are estimated from X.
%        resid = [n x q] matrix of residuals.
%

% RE Strauss, 4/21/98
%   9/20/99 - update handling of null input arguments.

function [b,pred,resid] = ltsregr(X,Y,Xp,noint)
  if (nargin < 3) Xp = []; end;
  if (nargin < 4) noint = []; end;

  calc_pred = 0;                          % Output flags
  calc_resid = 0;
  calc_d = 0;

  [n,p] = size(X);
  [ny,q] = size(Y);

  if (nargout >= 2)                       
    calc_pred = 1;
  end;
  if (nargout >= 3)
    calc_resid = 1;
    resid = zeros(n,q);
  end;

  maxiter = 100;                          % Max iterations
  iter = 0;

  if (n==1)                               % Transpose row vectors into col vectors
    X = X';
    [n,p] = size(X);
  end;
  if (ny==1)
    Y = Y';
    [ny,q] = size(Y);
  end;

  if (n ~= ny)
    error('LTSREGR: X,Y must have same number of observations');
  end;

  given_Xp = 0;
  if (~isempty(Xp))
    given_Xp = 1;
  end;
  if (isempty(noint))
    noint = 0;
  end;

  if (~noint)
    X = [ones(n,1) X];                % Augment the X-matrix
    b = zeros(p+1,q);                 % Allocate parameter matrix
  else
    b = zeros(p,q);                   % Allocate parameter matrix
  end;

  maxobs = ceil(n/2);

  for c = 1:q                         % Regress each column of Y on X
    y = Y(:,c);                         % Current Y vector
    b_fit = [];

    b_init = inv(X'*X)*X'*y;            % Initial regression
    pr = X * b_init;                    % Predicted values
    e = y - pr;                         % Residuals

    [e,i] = sort(e.*e);                 % Find smallest squared resids
    XX = X(i(1:maxobs),:);              % Save corresponding obs
    yy = y(i(1:maxobs),:);
    b_fit = inv(XX'*XX)*XX'*yy;         % Robust regression
    pr = X * b_fit;                     % Predicted values
    e = y - pr;                         % Residuals

    while ((b_fit ~= b_init) & (iter < maxiter))
      iter = iter+1;
      b_init = b_fit;
      [e,i] = sort(e.*e);
      XX = X(i(1:maxobs),:);
      yy = y(i(1:maxobs),:);
      b_fit = inv(XX'*XX)*XX'*yy;
      pr = X * b_fit;
      e = y - pr;
    end;

    if (iter >= maxiter)
      disp(['  Warning: max iterations exceeded for dependent variable ',int2str(c),';']);
      disp('           solution might not be stable');
    end;

    b(:,c) = b_fit;                     % Accumulate parameters

    if (calc_resid)                   % Save residuals
      resid(:,c) = e;
    end;
  end;

  if (calc_pred)                      % Calc predicted values for original X
    if (given_Xp)                       % Or for substitute matrix, if given one
      X = Xp;
      [nx,p] = size(X);
      if (~noint)
        X = [ones(nx,1) X];               % Augment
      end;
    end;
    pred = X*b;                         % Predicted values
  end;

  return;
