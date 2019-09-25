% Fratio: Randomizes the F-test for equality of two variances.  Creates a null 
%         distribution by permuting group membership.
%
%     Usage: [F,pr] = fratio(x1,x2,iter)
%
%         x1 =   [n x 1] column vector of data values for the first distribution.
%         x2 =   [m x 1] column vector of data values for the second distribution.
%         iter = optional number of randomization iterations [default = 5000].
%         ------------------------------------------------------------------------
%         F =    observed F ratio of larger variance to smaller variance.
%         pr =   significance level.
%

% RE Strauss, 6/4/98

function [F,pr] = fratio(x1,x2,iter)
  if (nargin < 3) iter = []; end;

  if (isempty(iter))
    iter = 5000;
  end;

  [r1,c1] = size(x1);
  [r2,c2] = size(x2);

  if (r1==1)
    x1 = x1';
  end;
  if (r2==1)
    x2 = x2';
  end;

  v1 = var(x1);
  v2 = var(x2);

  if (v2 > v1)
    x = x1;
    x1 = x2;
    x2 = x;
    v = v1;
    v1 = v2;
    v2 = v;
  end;

  F = v1/v2;

  if (iter)
    pr = 0;

    x = [x1;x2];
    g = makegrps([1 2],[length(x1) length(x2)]);
    ng = length(g);

    for it = 1:iter
      g = g(randperm(ng));
      v1 = var(x(g==1));
      v2 = var(x(g==2));
      Fit = v1/v2;
      if (Fit >= F)
        pr = pr + 1/iter;
      end;
    end;
  end;

  return;


