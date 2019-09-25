% RANDPROP: Generates a vector of random proportions, given the mean
%           proportion and size of the sampling distribution.
%           Uses arcsine transformation.
%
%     Usage: p = randprop(meanp,n,nval)
%
%           meanp = mean of sampling distribution.
%           n =     size of sampling distribution.
%           nval =  number of random proportions to be generated [default=1].
%           -----------------------------------------------------------------
%           p =     [nval x 1] column vector of random proportions.
%

% Sokal,RR & FJ Rohlf. 1981. Biometry, 2nd ed, pp. 427-428.  Freeman.

% RE Strauss, 10/6/95

function p = randprop(meanp,n,nval)
  if (nargin<3)
    nval = 1;
  end;

  meanp = asin(sqrt(meanp));          % Transform mean p
  randp = randn(nval,1);              % Vector of random-normal variates
  stderr = sqrt(1/(4*n));             % Expected standard error (radians)
  p = randp*stderr + meanp;           % Transformed random proportions

  L = find(p<=0);                     % Values less than lower bound
  G = find(p>=pi/2);                  % Values greater than upper bound
  lenL = length(L);
  lenG = length(G);

  while (lenL+lenG > 0)               % Replace vals out-of-bound with new ones
    if (lenL>0)
      randp = randn(lenL,1);
      p(L) = randp*stderr + meanp;
      L = find(p<=0);
      lenL = length(L);
    end;
    if (lenG>0)
      randp = randn(lenG,1);
      p(G) = randp*stderr + meanp;
      G = find(p>=pi/2);
      lenG = length(G);
    end;
  end;

  p = sin(p).^2;                      % Convert from arcsin sqrt(p) to proportions

  return;
