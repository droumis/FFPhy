% RANDNT: Generate random observations from a symmetrically truncated normal 
%         distribution.
%
%     Usage: randdata = randnt(r,c,l,u,{alpha})
%
%         r,c =   numbers of rows and columns of resulting random matrix.
%         l,u =   lower and upper limits of range.
%         alpha = central area of normal distribution to be represented by range 
%                   [default = 0.95].
%         -----------------------------------------------------------------------
%         randdata = [r x c] matrix of random numbers.
%

function randdata = randnt(r,c,l,u,alpha)
  if (nargin < 5)
    alpha = [];
  end;

  if (isempty(alpha))
    alpha = 0.05;
  end;
  if (alpha > 1)
    alpha = 1-(0.01 * alpha);
  end;

  if (l > u)
    [l,u] = exchange(l,u);
  end;

  lim =   norminv(1-alpha/2);
  w =     2*lim;
  rng =   u-l;
  cells = r*c;

  randdata = randn(cells,1);
  i = abs(randdata) > lim;

  while (any(i))
    rd = randn(sum(i),1);
    randdata(i) = rd;
    i = abs(randdata) > lim;
  end;

  randdata = randdata*rng/w + lim*rng/w + l;
  randdata = reshape(randdata,r,c);

  return;

