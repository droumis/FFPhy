% SWEEPREG: Multiple regression via the sweep operator, given means and
%           covariance matrix.
%
%       Usage: [D,E,F] = sweepreg(m,C,iv)
%
%           m =  [p] vector of means
%           C =  [p x p] covariance matrix
%           iv = [q] vector specifying indices of independent variables; the
%                 indices of the r dependent variables are complementary to
%                 this list (q+r=p)
%           ----------------------------------------------------------------
%           D =  [q+1 x q+1] matrix; can be used to estimate sampling vars
%                 and covars of regression coefficients
%           E =  [q+1 x r] matrix of regression coefficients; col j gives
%                 intercept (b0) and slopes (b1..b(p-q)) of regression of
%                 dependent var j on independent vars 1...q
%           F =  [r x r] matrix of residual covariances
%

% Little,RJA & DB Rubin. 1987. Statistical Analysis with Missing Data. Wiley.
%   Section 6.5, pp. 112-115.

% RE Strauss, 7/5/95

function [D,E,F] = sweepreg(m,C,iv)
  i = size(iv,1);                     % M & iv should be row vectors
  if (i>1)
    iv = iv';
  end;
  i = size(m,1);
  if (i>1)
    m = m';
  end;

  p = length(m);                      % p = number of vars
  q = length(iv);                     % q = number of independent vars

  i = ones(1,p);                      % Find indices of dependent vars
  i(iv) = i(iv)-1;
  dv = find(i);
  r = length(dv);                     % r = number of dependent vars

  A = zeros(p+1,p+1);                 % Construct augmented covariance matrix
  A(1,1) = -1;                          % First cell is constant -1
  A(1,2:p+1) = m;                       % First rows/cols are means vectors
  A(2:p+1,1) = m';
  A(2:p+1,2:p+1) = C;                   % Remainder is covar matrix

  indx = [1 iv+1 dv+1];               % Rearrange matrix so indep vars come first
  A = A(indx,indx);

  A = sweep(A,1:q+1);                 % Sweep const + independent vars from matrix
  D = A(1:q+1,1:q+1);                 % Extract D
  E = A(1:q+1,q+2:p+1);               % Extract E (regression coeffs)
  F = A(q+2:p+1,q+2:p+1);             % Extract F (residual covar matrix)

  return;
