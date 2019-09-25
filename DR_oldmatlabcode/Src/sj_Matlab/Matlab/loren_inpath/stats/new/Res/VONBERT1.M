% VONBERT1: Fits and optionally bootstraps a von Bertalanffy growth function to 
%           k sets of size-at-age data, one size value per age per set.  
%           Ages must be the same for all k sets of size-at-age data.
%           Returns two vectors of parameters: a 3x1 vector of P1=[Sj,Sk,c] and 
%           a 22x1 vector of P2=[K,S0,Sa,t0,t10,...,t90,r10,...,r90].
%
%     Usage: [P1,P2] = vonbert1(A,S)
%
%           P1 = [k x 3]  matrix of parameters
%           P2 = [k x 22] matrix of parameters
%           -------------------------------------------------
%           A =  vector of ages 
%           S =  [k x length(A)] matrix of size-at-age data
%

% RE Strauss, 4/3/98
%   1/4/00 -    changed fminu() to fmins().

function [P1,P2] = vonbert1(A,S)
  na = length(A);
  [ksets,p] = size(S);

  if (p ~= na)
    error('  VONBERT1: Age vector and size-at-age matrix incompatible');
  end;

  if (nargin > 1)
    make_P2 = 1;
  else
    make_P2 = 0;
  end;

  P1 = zeros(ksets,3);                % Allocate parameter matrices
  if (make_P2)
    P2 = zeros(ksets,22);
  end;
  W = ones(1,na);                     % Weights vector

  for ki = 1:ksets                    % Cycle thru data sets
    init = [S(ki,1),S(ki,na),0.5];
    p1 = fmins('vbfunc',init,[],[],A,S(ki,:),W); % Fit function
    P1(ki,:) = p1;

    if (make_P2)
      Sj = P1(1);
      Sk = P1(2);
      c =  P1(3);
      j = 1;
      k = na;

      K = (j-k)/(A(k)-A(j))*log(c);
      S0 = Sj + (Sk-Sj)*(1-c^(-j))/(1-c^(k-j));
      Sa =(Sk-Sj*c^(k-j))/(1-c^(k-j));
      t0 = A(j) - ((A(k)-A(j))/((k-j)*log(c))) * log((Sk-Sj)/(Sk-Sj*c^(k-j)));

      P2 = [K S0 Sa t0];
    end;
  end;

  return;
