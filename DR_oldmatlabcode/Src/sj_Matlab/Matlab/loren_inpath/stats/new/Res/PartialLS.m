% PartialLS:  Partial least squares via the SIMPLS algorithm for nonlinear iteration (de Jong 1993).  
%             This version zero-centers X and Y by column.  The number of extracted factors is 
%             determined by jackknifing the data and chosing the number of factors that minimizes
%             the total prediction error (PMSE).
%
%     Usage: [b] = partialls(X,Y,{cmax})
%
%         X =     [n x p] matrix of independent (predictor) variables.
%         Y =     [n x q] matrix of dependent (response) variables.
%         cmax =  maximum number of factors to be extracted [default = max([p,q])].
%         -------------------------------------------------------------------------
%         b
%
%     <<< NOT COMPLETED >>>
%

% RE Strauss, 11/20/02

function b = partialls(X,Y,cmax)
  if (nargin < 3) cmax = []; end;
  
  [n,p] = size(X);
  [m,q] = size(Y);
  
  if (isempty(cmax))
    cmax = max([p,q]);
  end;
  
  if (n~=m)
    error('  PartialLS: number of observations for X and Y must be equal.');
  end;
  
  X = zcenter(X)
  Y = zcenter(Y)

  A = X'*Y
  M = X'*X
  C = eye(size(A))
  W = [];
  P = [];
  Q = [];
  
  for h = 1:cmax
h    
    e = eigen(A'*A)
    q = e(:,1)                       
    w = A*q
    c = w'*M*w
    w = w/sqrt(c)
    W = [W w]                           % X loadings
    p = M*w
    P = [P p]
    q = A'*w
    Q = [Q q]                           % Y loadings
    v = C*p
    v = v/norm(v)
    M = M-p*p'
    A = C*A
    C = C-v*v'
    T = X*W                             % X scores
    B = W*Q'          
  end;

  return;
  