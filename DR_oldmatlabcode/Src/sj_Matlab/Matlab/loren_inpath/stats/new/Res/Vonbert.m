% VONBERT: Fits and optionally bootstraps a von Bertalanffy growth function to 
%          k sets of size-at-age data.  
%          Ages must be the same for all k sets of size-at-age data.
%
%     Usage: [P,pS,tp,r,curve] = vonbert(A,S,{doplot})
%
%           A =     vector (length n) of ages.
%           S =     [k x length(A)] matrix of k sets of size-at-age data.  The 
%                     matrix may include missing data (NaN) but, if so, the Sj and 
%                     Sk parameters might not be equivalent across all data sets.
%           doplot = optional boolean flag indicating whether plot is to be produced
%                     [default = 0].
%           -----------------------------------------------------------------------
%           P =     [k x 7] matrix of parameters: c, K, t0, S0, Sj, Sk, Sa.
%           pS =    [k x n] matrix of predicted sizes at each input age.
%           tp =    [k x 9] matrix of percent-saturation times: t10, t20, ..., t90.
%           r =     [k x n] matrix of instantaneous growth rates at each age.
%           curve = [k x 100] matrix of predicted sizes ages linspace(t0,A(n)), 
%                     suitable for plotting.
%

% RE Strauss 10/22/98
%   11/5/98 - allow missing data in S matrix.
%   1/4/00 -  changed fminu() to fmins().
%   3/1/00 -  allow A & S to be column vectors.
%   2/28/02 - changed fmins() to fminsearch(); make plot optional.

function [P,pS,tp,r,curve] = vonbert(A,S,doplot)
  if (nargin < 3) doplot = []; end;
  
  if (isempty(doplot))
    doplot = 0;
  end;

  if (isvector(A))
    A = [A(:)]';
  else
    error('  VONBERT: age must be vector');
  end;
  if (isvector(S))
    S = [S(:)]';
  end;
  
  nage = length(A);
  [ksets,p] = size(S);

  make_pS = 0;                        % Determine which output parameters have been
  make_tp = 0;                        % requested, set flags and allocate matrices
  make_r =  0;
  make_curve = 0;

  if (nargout > 1)
    make_pS = 1;
    pS = zeros(ksets,nage);
  end;
  if (nargout > 2)
    make_tp = 1;
    tp = zeros(ksets,9);
  end;
  if (nargout > 3)
    make_r = 1;
    r = zeros(ksets,nage);
  end;
  if (nargout > 4)
    make_curve = 1;
    curve = zeros(ksets,100);
  end;

  if (p ~= nage)
    error('  VONBERT: Age vector and size-at-age matrix incompatible');
  end;

  W = ones(1,nage);                   % Weights vector

  for ki = 1:ksets                    % Cycle thru data sets
    size = S(ki,:);

    fvals = find(finite(size));
    size = size(fvals);
    age = A(fvals);
    w = W(fvals);
    nage = length(age);

    init = [size(1),size(nage),0.5];
    p1 = fminsearch('vbfunc',init,[],age,size,w); % Fit function

    si = p1(1);
    sk = p1(2);
    c =  p1(3);

    ti = age(1);
    tk = age(nage);

    cexp = c.^(tk-ti);

    K = -log(c);
    t0 = ti - log((sk-si)/(sk-si*cexp)) / log(c);
    Sa = (sk-si*cexp)/(1-cexp);
    s0 = si + (sk-si)*(1-c.^(-ti))/(1-cexp);

    P(ki,:) = [c K t0 s0 si sk Sa];

    if (make_pS)                      % Predicted size-at-age
      pS(ki,:) = Sa .* (1-exp(-K.*(A-t0)));
    end;

    if (make_tp)
      p = 0.1:0.1:0.9;                % Percent-saturation times
      tp(ki,:) = t0 + log(1-p)/log(c);
    end;

    if (make_r)
      r(ki,:) = log(c) .* (si-sk).*(c.^A)./(c^ti-c^tk);
    end;
  end;

  if (make_curve)
    K =  P(:,2);
    t0 = P(:,3);
    Sa = P(:,7);

    minage = min(t0);
    maxage = max(A);
    t = linspace(minage,maxage);

    for ki = 1:ksets                    % Cycle thru data sets
      curve(ki,:) = Sa(ki) .* (1-exp(-K(ki).*(t-t0(ki))));
    end;

    maxcurve = max(curve(:));
    deltax = 0.05*(maxage-minage);
    deltay = 0.05*(maxcurve);

    if (doplot)
      hold on;
      plot(t,curve(1,:),'k');
      plot(A,S(1,:),'ok');
      v = putbnds([t,A]',[curve,S]');
      if (ksets>1)
        for ki = 2:ksets
          plot(t,curve(ki,:),'k');
          plot(A,S(ki,:),'ok');
          vki = putbnds([t,A]',[curve(ki,:),S(ki,:)]');
          v = max([v;vki]);
        end;
        axis(v);
      end;
      hold off;
    end;
  end;

  return;
