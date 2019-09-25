% LSTREEG: Gradients (derivatives) dF/dd for lstree().
%
%     Usage: G = lstreeg(D,r,delta)
%
%         D =     current estimate of optimal distance matrix (vector form).
%         r =     current PD weight.
%         delta = observed distance matrix (vector form).
%         -------------------------------------------------------------------
%         F(D,r), L(D), P(D), r, and partial derivatives G are given by 
%           De Soete (1983).
%

% RE Strauss, 6/3/99

function G = lstreeg(D,r,delta)
  Dmat = trisqmat(D);
  N = size(Dmat,1);

  G = -2*(D-delta);                     % Init with first term of derivatives

  comblist = combvals(N,4);             % List of all combs of four taxa
  for i = 1:comb(N,4)
    c = comblist(i,:);
    d = Dmat(c,c);
    di = [c(1)*ones(1,3) c(2)*ones(1,2) c(3)];
    dj = [c(2:4) c(3:4) c(4)];
    [d,j] = sort(trilow(d));
    di = di(j);
    dj = dj(j);
Dmat
c
j
di_dj = [di' dj']
d

    for g = 1:length(G)

    end;

pause
  end;

  return;
