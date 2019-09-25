% RANDCORR: Returns a randomized Pearson product-moment correlation matrix
%           by sampling a randomized data matrix.
%             If only a sample size is given, returns a matrix sampled from
%           null distribution (centered on zero).
%             If a correlation matrix is given, returns a matrix sampled from
%           corresponding sampling distributions (centered on the coefficients 
%           provided).  If a covariance matrix is given, retu rns a randomized 
%           covariance matrix.
%
%       Syntax: R = randcorr(P|C,N)
%
%           P =   order of randomized output matrix
%             OR
%           C =   correlation or covariance matrix
%           N =    scalar sample size
%           ---------------------------------------------
%           R =   randomized correlation or covariance matrix
%

% RE Strauss, 1/8/96

function R = randcorr(C,N)
  if (nargin<2) N = []; end;

  if (isempty(N))
    error('  RANDCORR: too few input arguments');
  end;

  [r,c] = size(C);
  if (r==1 & c==1)
    P = C;
    inmat = 0;
  else
    inmat = 1;

    is_sym = 1;
    if (r~=c)
      is_sym = 0;
    else
      if (sum(sum(abs(c-c')))>eps)
        is_sym = 0;
      end;
    end;
    if (~is_sym)
      error('  RANDCORR: input must be covariance or correlation matrix.');
    end;

    if (all(diag(C)==1))
      is_cov = 0;                       % Input correlation matrix
    else
      is_cov = 1;                       % Input covariance matrix
    end;
  end;

  % Sample from null distributions

  if (~inmat)
    x = randn(N,P);
    R = corrcoef(x);
  end;

  % Sample from sampling distributions

  if (inmat)
    P = r;
    x = randmvn(N,zeros(1,P),C);
    if (is_cov)
      R = cov(x);
    else
      R = corrcoef(x);
    end;
  end;

  return;

