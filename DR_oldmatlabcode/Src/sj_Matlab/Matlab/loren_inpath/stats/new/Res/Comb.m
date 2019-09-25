% COMB: Calculates the number of possible combinations of N objects
%       taken R at a time.  If N and R are compatible vectors, a
%       corresponding vector is returned based on element-wise operations.
%       If N is a scalar and R a vector, N is expanded to the same size as R.
%
%     Usage: c = comb(N,R)
%

% RE Strauss, 1/13/96
%   11/4/98 - for scalar N and vector R, expand N to vector.
%   4/23/99 - change eqn to use logs to increase range of solutions.
%   1/17/01 - use samelength() to expand and check for compatalibity.

function c = comb(N,R)
  if (nargin<2)
    error('  COMBIN: Too few input arguments');
  end;

  [b,N,R] = samelength(N,R);
  if (~b)
    error('  COMBIN: Argument vectors incompatible');
  end;

  c = zeros(size(N));
  for i=1:length(N)
    if (N(i)>=R(i))
      c(i) = round(exp(sum(log(1:N(i)))-(sum(log(1:(N(i)-R(i))))+sum(log(1:R(i))))));
    else
      c(i) = 0;
    end;
  end;
  return;
