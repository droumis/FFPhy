% RANDAUTP: For a given specified level of autocorrelation, predicts the
%           corresponding value of c to be used to generate autocorrelated
%           uniform random numbers, based on the method of Willemain & Desautels 
%           (1993).  Adjusts for sample size, if supplied.  Called by randauto().
%
%     Usage: c = randautp(r,{n})
%
%         r = correlation coefficient.
%         n = optional sample size.
%         --------------------------------
%         c = autocorrelation coefficient.
%

% Willemain, TR & PA Desautels. 1993. A method to generate autocorrelated
%   uniform random numbers. J. Statist. Comput. Simul. 45:23-31.

% RE Strauss, 2/11/99
%   6/15/00 - added optional adjustment for sample size.

function c = randautp(r,n)
  if (nargin < 2) n = []; end;

  a0 =  0.83324809;
  a1 = -1.038672;
  a2 =  0.34935697;
  a3 = -0.24098563;
  a4 =  0.20310172;
  a5 =  0.010787407;
  a6 = -0.027993053;

  if (isempty(n))
    r = log(corrz(r));
  else
    r = log(corrz(r,n));
  end;

  c = a0 + a1*r.^1 + a2*r.^2 + a3*r.^3 + a4*r.^4 + a5*r.^5 + a6*r.^6;

  [n,p] = size(c);
  if (ismatrix(c))
    [i,j] = find(c<=0);
    if (~isempty(i))
      for k = 1:length(i)
        c(i(k),j(k)) = eps;
      end;
    end;
  else
    i = find(c<=0);
    if (~isempty(i))
      if (n==1)
        c(i) = eps*ones(1,length(i));
      else
        c(i) = eps*ones(length(i),1);
      end;
    end;
  end;

  return;

