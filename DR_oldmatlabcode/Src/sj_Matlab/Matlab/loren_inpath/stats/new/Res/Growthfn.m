% GROWTHFN: Objective function for growth().  Returns the SSE (sum of squared 
%           residuals), given the identification and parameters of the model and 
%           vectors of age and size measures.  Also returns the predicted values 
%           corresponding to ages.  
%           See Wan, Zhong & Wang 1998 for notation.
%           

% RE Strauss, 10/21/98

function [sse,predW] = growthfn(p,m,t,W)
  p = abs(p);           % All parameters must be positive

  if (m==1)             % von Bertalanffy, =monomolecular (3 parameters: k, W0, Wf)
    k =  p(1);
    W0 = p(2);
    Wf = p(3);
    predW = Wf - (Wf-W0).*exp(-k.*t);

  elseif (m==2)         % Gompertz (3 parameters: k, W0, D)
    k =  p(1);
    W0 = p(2);
    D =  p(3);
    predW = W0.*exp(k.^(1-exp(-D.*t)));

  elseif (m==3)         % logistic (3 parameters: k, W0, Wf)
    k =  p(1);
    W0 = p(2);
    Wf = p(3);
    predW = (W0*Wf)./(W0+(Wf-W0).*exp(-k.*t));

  elseif (m==4)         % Richards (4 parameters: k, W0, Wf, n)
    k =  p(1);
    W0 = p(2);
    Wf = p(3);
    n =  p(4);
    predW = (W0*Wf)./((W0^n+(Wf^n-W0^n).*exp(-k.*t)).^(1./n));

  elseif (m==5)         % Janoschek (4 parameters: k, W0, Wf, p)
    k =  p(1);
    W0 = p(2);
    Wf = p(3);
    p =  p(4);
    predW = Wf - (Wf-W0).*exp(-k.*t.^p);

  elseif (m==6)         % Hill (4 parameters: k, W0, Wf, n)
    k =  p(1);
    W0 = p(2);
    Wf = p(3);
    n =  p(4);
    predW = (W0*k^n+Wf.*t.^n)./(k^n+t.^n);

  elseif (m==7)         % Wan (4 parameters: k, a, b, c)
    k =  p(1);
    a =  p(2);
    b =  p(3);
    c =  p(4);
    predW = a - 1./(b+c.*exp(k.*t));

  elseif (m==8)         % France (5 parameters: k, W0, Wf, c, T)
    k =  p(1);
    W0 = p(2);
    Wf = p(3);
    c =  p(4);
    T =  p(5);
    predW = Wf - (Wf-W0).*exp(-k*(t-T)+2*c*(sqrt(t)-sqrt(T)));

  else
    error('Growthfn: invalid model identifier');
  end;

  sse = sum((W-predW).^2);

  return;
