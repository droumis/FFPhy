% VONBERT2: Produces the same output as function 'vonbert', but based upon 
%           previously fitted von Bertalanffy parameters rather than original data.  
%
%     Usage: [P,pS,tp,r,curve] = vonbert2(param,ages)
%
%           param = [k x 3] matrix of von Bertalanffy parameters (K,L,t0) for each 
%                     of k growth curves.
%           A =     vector of n ages at which predicted sizes and growth rates
%                     are to be estimated.
%           -----------------------------------------------------------------------
%           P =  [k x 7] matrix of parameters: c, K, t0, S0, Sj, Sk, Sa.
%           pS = [k x n] matrix of predicted sizes at each input age.
%           tp = [k x 9] matrix of %-saturation times: t10, t20, ..., t90.
%           r =  [k x n] matrix of instantaneous growth rates at each age.
%           curve = [k x 100] matrix of predicted sizes ages linspace(t0,A(n)), 
%                     suitable for plotting.
%

% RE Strauss 11/6/98, modified from vonbert().

function [P,pS,tp,r,curve] = vonbert2(param,A)
  nage = length(A);
  ksets = size(param,1);

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

  for ki = 1:ksets                    % Cycle thru data sets
    K = param(ki,1);
    Sa = param(ki,2);
    t0 = param(ki,3);

    predsize = Sa .* (1-exp(-K.*(A-t0)));   % Predicted size-at-age

    c = exp(-K);
    si = predsize(1);
    sk = predsize(nage);

    ti = A(1);
    tk = A(nage);

    cexp = c.^(tk-ti);

    s0 = si + (sk-si)*(1-c.^(-ti))/(1-cexp);

    P(ki,:) = [c K t0 s0 si sk Sa];

    if (make_pS)                      % Predicted size-at-age
      pS(ki,:) = predsize;
    end;

    if (make_tp)
      p = 0.1:0.1:0.9;                % Percent-saturation times
      tp(ki,:) = t0 + log(1-p)/log(c);
    end;

    if (make_r)
      r(ki,:) = log(c) .* (si-sk).*(c.^A)./(c^ti-c^tk);
    end;

    if (make_curve)
      t = linspace(t0,A(length(A)));
      curve(ki,:) = Sa .* (1-exp(-K.*(t-t0)));
    end;
  end;

if (make_curve)
  close all;
  plot(t,curve(1,:));
  hold on;
  for ki = 2:ksets
    plot(t,curve(ki,:));
  end;
  hold off;
end;

  return;
