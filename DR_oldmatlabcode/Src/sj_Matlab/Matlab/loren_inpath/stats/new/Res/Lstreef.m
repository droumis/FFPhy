% LSTREEF: Objective function for lstree().
%
%     Usage: [FDr,LD,PD,r] = lstreef(D,r,delta)
%
%         D =     current estimate of optimal distance matrix (vector form).
%         r =     current PD weight.
%         delta = observed distance matrix (vector form).
%         -------------------------------------------------------------------
%         F(D,r), L(D), P(D), and r are given by De Soete (1983).
%

% RE Strauss, 6/3/99

function [FDr,r1] = lstreef(D,r0,delta)
  Dmat = trisqmat(D);
  N = size(Dmat,1);

  d = delta-D;
  LD = d'*d;

  PD = 0;
  comblist = combvals(N,4);             % List of all combs of four taxa
  for i = 1:comb(N,4)
    c = comblist(i,:);
    d = Dmat(c,c);
    d = sort(trilow(d));
    PD = PD + [0 0 -1 -1 1 1]*d;
  end;

  FDr = LD + r0*PD;

  r1 = r0;
  if (abs(PD)>eps)
    r1 = LD/PD;
  end;
  
  return;
