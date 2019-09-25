% AUTOCORR: Computes the autocorrelations along the columns of a matrix X,
%           as a function of lag.  Uses the Pearson product-moment
%           correlation coefficient.
%
%     Syntax: [R,prob,CI_low,CI_high] = autocorr(X,{maxlag},{concat},{iter},{CI_level})
%
%           X =         [n x p] data matrix.
%           maxlag =    maximum lag (defaults to n-3).
%           concat =    1 (TRUE) if vector is to be repeated, assuming continuing of
%                           pattern into a second time sequence.
%                       0 (FALSE) if vector is not to be repeated [default[.
%           iter =      number of randomization iterations for probabilities and
%                         confidence intervals [default = 0].
%           CI_level =  percent width of confidence intervals [default = 95].
%           ------------------------------------------------------------------------
%           R =         [maxlag x p] matrix of autocorrelation estimates.
%           prob =      [maxlag x p] vector of probabilities for the observed
%                           autocorrelation estimates; asymptotic if rand_iter=0,
%                           randomized if rand_iter>0.
%           CI_low,CI_high = [maxlag x p] matrices of asymmetric CI% confidence
%                         limits (low,high) of the autocorrelation coefficients;
%                         asymptotic if boot_iter=0, bootstrapped if boot_iter>0.
%

% Should be made more efficient by bootstrapping all cols at once, rather than 
% consecutively.

% RE Strauss, 12/21/95
%   6/23/00 - improve handling of input arguments.
%   1/21/01 - updated call to corr(); removed support for null limits.

function [R,prob,CI_low,CI_high] = autocorr(X,maxlag,concat,iter,CI_level)

  if (nargin < 2) maxlag = []; end;
  if (nargin < 3) concat = []; end;
  if (nargin < 4) iter = []; end;
  if (nargin < 5) CI_level = []; end;

  [N,P] = size(X);
  if (N == 1)                             % If X is row vector,
    X = X';                               %   transpose to column vector
    [N,P] = size(X);
  end;

  default_maxlag = N-3;
  default_concat = 0;
  default_iter = 0;
  default_CI_level = 95;

  if (isempty(CI_level))                  % Set arguments to default values
    CI_level = default_CI_level;          %   if not passed
  end;
  if (isempty(iter))
    iter = default_iter;
  end;
  if (isempty(concat))
    concat = default_concat;
  end;
  if (isempty(maxlag))
    maxlag = default_maxlag;
  end;

  % Calculate autocorrelation coefficients

  R = zeros(maxlag,P);                    % Allocate output matrices
  prob = zeros(maxlag,P);
  CI_low = zeros(maxlag,P);
  CI_high = zeros(maxlag,P);

  if (concat)
    X = [X;X];                            % Double the matrix
    v2 = zeros(N,maxlag);                 % Allocate lag-data matrix
    for j = 1:P;                          % Cycle thru variables
      v1 = X(1:N,j);                      % Original vector
      for lag = 1:maxlag                  % Iterate thru lag values
        v2(:,lag) = X((lag+1):(lag+N),j);   % Contrast vectors
      end;

      [r,pr,CIl,CIh] = corr(v1,v2,[],iter,CI_level);
      
      R(:,j) = r';
      prob(:,j) = pr';
      CI_low(:,j) = CIl';
      CI_high(:,j) = CIh';
    end; 
  end; 

  if (~concat)
    for j = 1:P;                          % Cycle thru variables
      for lag = 1:maxlag                  % Iterate thru lag values
        vectlen = N-lag;
        v1 = X(1:vectlen,j);                % Subset two vectors,
        v2 = X((lag+1):(vectlen+lag),j);    %   displaced by the lag

        [r,pr,CIl,CIh] = corr(v1,v2,[],iter,CI_level);

        R(lag,j) = r;                       %   Stash values in output matrices
        prob(lag,j) = pr;
        CI_low(lag,j) = CIl;
        CI_high(lag,j) = CIh;
      end;  
    end;  
  end;  

  return;
