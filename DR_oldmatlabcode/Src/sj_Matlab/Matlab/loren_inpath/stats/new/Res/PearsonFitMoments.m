% PearsonFitMoments: Finds the Pearson-system distribution corresponding to a given
%           set of moment statistics (mean, var, skewness, kurtosis)
%

% RE Strauss, 4/20/95
%   8/20/99 - changed plot colors for Matlab v5.
%   11/19/02 - renamed from pearcurv().

function PearsonFitMoments(M,V,S,K)
  b1 = S*S;
  b2 = K;

  d = 10*b2 - 12*b1 - 18;

  a = S * (b2+3) * sqrt(V) / d;
  c0 = (4*b2 - 3*b1) * V / d;
  c1 = a;
  c2 = (2*b2 - 3*b1 - 6) / d;

  c = [a c0 c1 c2]                   % Vector of constants
  
  t0 = -10;                            % Arbitrary range
  tfinal = +10;
  p0 = 1e-6;                          % Arbitrarily small value
  [X,P] = ode('peardot',t0,tfinal,p0,c); % Construct distribution
  P = P./trapz(X,P);                  % Standardize area --> 1

  XI = [min(X) : (max(X)-min(X))/100 : max(X)]'; % Interpolate
  PI = interp1(X,P,XI,'spline');
%[XI PI]

  p_tol = 0.002;
  indx = find(PI == max(PI));           % Find mode
  modex = indx(1);
  ub = modex;
  for i = (modex+1):length(XI)
    if (PI(i) > p_tol)
      ub = ub+1;
    else
      break;
    end;
  end;
  lb = modex;
  for i = modex:-1:1
    if (PI(i) > p_tol)
      lb = lb-1;
    else
      break;
    end;
  end;
%[1 lb modex ub length(X)]

  XI = XI(lb:ub);                      % Capture central part of distribution
  PI = PI(lb:ub);

  plot(XI,PI,'k');

  d = 4*(4*K-3*S)*(2*K-3*S-6);
  kappa = sign(d)*(S*(K+3)^2)/max(eps,abs(d))

  return;
