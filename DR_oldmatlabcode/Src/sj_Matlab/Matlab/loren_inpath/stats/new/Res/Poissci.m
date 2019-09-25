% POISSCI:  Finds the bounds of X delimiting a central CI of the Poisson 
%           distribution.  The bounds may be specified as conservative integers, 
%           non-conservative integers, or linearly interpolated real numbers.  
%           Note that, if the lower bound is found to be zero, the CI may be 
%           non-conservative (i.e., the realized alpha may be greater than 
%           nominal).
%
%     Usage: [lb,ub,rci] = poissci(lambda,{kind},{ci})
%
%         lambda =  matrix of expected counts (or observed mean counts) of X.
%         kind =    optional flag indicating the kind of interval desired:
%                     0 = linearly interpolated non-integer bounds [default]
%                     1 = conservative integer bounds
%                     2 = non-conservative integer bounds
%         ci =      optional level of central confidence interval
%                     [default = 0.95].
%         ------------------------------------------------------------------
%         lb =      corresponding matrix of lower bounds.
%         ub =      corresponding matrix of upper bounds.
%         rci =     corresponding matrix of realized confidence intervals.
%

% RE Strauss, 4/3/00
%   6/1/01 -  corrected problem with lower bound for low expected counts; lb must 
%               now be >= 0;
%             changed sequence of input arguments.

function [lb,ub,rci] = poissci(lambda,kind,ci)
  if (nargin < 2) kind = []; end;
  if (nargin < 3) ci = []; end;

  if (isempty(ci))
    ci = 0.95;
  end;
  if (isempty(kind))
    kind = 0;
  end;
  if (ci > 1)
    ci = ci/100;
  end;
  tailprob = (1-ci)/2;

  [r,c] = size(lambda);
  lb = zeros(r,c);
  ub = zeros(r,c);
  rci = ones(r,c);

  for i = 1:r
    for j = 1:c
      lamb = lambda(i,j);
      u = 10*ceil(lamb);
      p = poisscdf(0:u,lamb);             % Get Poisson pdf
      while (1-p(end) > tailprob)
        u = 2*u;
        p = poisscdf(0:u,lamb);
      end;

      x1 = max(find(p < tailprob))-1;     % Left tail
      if (isempty(x1))
        x1 = 0;
      end;
      x2 = x1+1;
      p1 = p(x1+1);
      p2 = p(x2+1);

      switch (kind)
        case 0,                           % Linearly interpolated
          if (p1 <= tailprob)
            lb(i,j) = x1 + (x2-x1)*(tailprob-p1)/(p2-p1);
            rci(i,j) = rci(i,j) - tailprob;
          else
            lb(i,j) = x1;
            rci(i,j) = rci(i,j) - p1;
          end;
        case 1,                           % Conservative
          lb(i,j) = x1;
          rci(i,j) = rci(i,j) - p1;

        case 2,                           % Non-conservative
          if (p1 <= tailprob)
            lb(i,j) = x2;
            rci(i,j) = rci(i,j) - p2;
          else
            lb(i,j) = x1;
            rci(i,j) = rci(i,j) - p1;
          end;

        otherwise
          error('  POISSCI: invalid kind-of-interval indicator');
      end;

      p = 1-p;                          % Right tail
      x1 = max(find(p > tailprob))-1;
      x2 = x1+1;
      p1 = p(x1+1);
      p2 = p(x2+1);

      switch (kind)
        case 0,                           % Linearly interpolated
          ub(i,j) = x1 + (x2-x1)*(tailprob-p1)/(p2-p1);
          rci(i,j) = rci(i,j) - tailprob;

        case 1,                           % Conservative
            ub(i,j) = x2;
            rci(i,j) = rci(i,j) - p2;

        case 2,                           % Non-conservative
          ub(i,j) = x1;
          rci(i,j) = rci(i,j) - p1;
      end;
    end;
  end;
  
  return;
