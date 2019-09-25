% RANDAUTO: Produces a vector of autocorrelated uniform random numbers,
%           using the method of Willemain & Desautels (1993).
%
%     Usage: U = randauto(n,r,is_c)
%
%           n =     number of values to be generated.
%           r =     level of correlation.
%           is_c =  boolean flag indicating that the value passed as 'r'
%                     is to be used as the value for c; if false [default],
%                     c is predicted from r.
%           --------------------------------------------------------------
%           U =     column vector of autocorrelated uniform random values.
%

% Willemain, TR & PA Desautels. 1993. A method to generate autocorrelated
%   uniform random numbers. J. Statist. Comput. Simul. 45:23-31.

% RE Strauss, 2/11/99
%   6/15/00 - pass sample size to randautp() for sample-size correction of 
%               correlation coefficient.

function U = randauto(n,r,is_c)
  if (nargin < 3) is_c = []; end;

  if (isempty(is_c))
    is_c = 0;
  end;

  if (~isscalar(n) | ~isscalar(r))
    error('  RANDAUTO: input values must be scalars');
  end;

  if (is_c)                       % Use given value as c
    c = abs(r);
  else                            % Else predict c from r
    if (r<-1 | r>1)
      error('  RANDAUTO: r out of range');
    end;

    c = randautp(abs(r),n);         
  end;

  if (c <= 0)
    error('  RANDAUTO: c out of range');
  end;

  V = rand(n+2,1);                % Uniform random variate
  U = zeros(n,1);                 % Allocate output vector

  X = V(1) + c*V(2);

  for i = 1:n
    if (c >= 1)                   % F(X) of current value
      if (X <= 1)
        FX = 0.5*X*X/c;
      elseif (X <= c)
        FX = (X-0.5)/c;
      else
        FX = 1 - 0.5*((1+c-X).^2)/c;
      end;
    else
      if (X <= c)
        FX = 0.5*X*X/c;
      elseif (X <= 1)
        FX = X-0.5*c;
      else
        FX = (1-0.5*c) + (1+1/c)*(X-1) - 0.5*(X*X-1)/c;;
      end;
    end;
    U(i) = FX;

    if (r>0)                      % Positive correlation
      X = FX  + c*V(i+2);
    else                          % Negative correlation
      X = 1 - FX + c*V(i+2);
    end;
  end;

  return;
