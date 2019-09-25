% VBFUNC: Objective function for fitting the von Bertalanffy function based
%         on the Schnute & Fournier parameters.
%
%     Usage: sse = vbfunc(param,A,S,W)
%
%         sse = sum of squared errors
%         param = vector of parameters [Sj,Sk,c]
%         A =     sorted vector of n ages
%         S =     corresponding vector of n sizes
%         W =     corresponding vector of n weights
%

% Ratkowsky (1986) parameterization

% RE Strauss, 6/2/95

function sse = vbfunc(param,A,S,W)
  Sj = param(1);                      % Get function parameters
  Sk = param(2);
  c =  param(3);
  n = length(A);

  if (c>=1 | c<=0)                    % Check bounds of parameters
    sse = 1e6;
    return;
  end;

  if (length(S)~=n | length(W)~=n)
    error('  Vector lengths incompatible');
  end;

  sum_wt = sum(W);
  min_age = A(1);
  range_age = A(n)-A(1);
  range_size = Sk-Sj;
  sse = 0;

  for i=1:n                           % Calculate error
    m = 1 + (n-1)*(A(i)-min_age)/range_age;
    S_hat = Sj + range_size*(1-c^(m-1))/(1-c^(n-1));
    sse = sse + (S(i)-S_hat)^2;
  end;

  return;
