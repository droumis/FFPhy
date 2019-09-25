% MANTEL: Mantel's test (Daniels 1944) for the positive association between two square 
%         symmetric distance matrices of identical size.  Diagonal values are 
%         assumed to be zeros; non-zero diagonal elements are ignored.
%         Estimates post-hoc power via a parametric bootstrap if standard errors 
%         are provided for the input distance matrices, and via nonparametric 
%         bootstrap if standard errors are not provided.
%
%     Syntax: [r,Z,pr,power] = mantel(X,Y,{iter},{Xstderr},{Ystderr})
%
%         X,Y =     square symmetric matrices with zeros on the diagonal.
%         iter =    optional number of permutation iterations [default=0].
%         Xstderr = optional square symmetric matrix of standard errors for the 
%                     distances in X, for power estimation.
%         Ystderr = optional standard-error matrix for Y.
%         ----------------------------------------------------------------------
%         r =       correlation between matrix cells.
%         Z =       sum of cross products of matrix cells.
%         pr =      right-tailed significance level (p-value), either 
%                     asymptotic (if iter=0) or randomized.
%         power =   randomized estimate of post-hoc power.
%

% Daniels, HE. 1944. Biometrika 33:129-135.

% RE Strauss, 10/27/95
%   11/19/99 - error messages changed; ignore non-zero diagonal values rather 
%               than treat as error condition.
%   10/19/00 - added power estimation.
%   6/18/02 -  added reference to original work by Daniels 1944.

function [r,Z,pr,power] = mantel(X,Y,iter,Xstderr,Ystderr)
  if (nargin < 3) iter = []; end;
  if (nargin < 4) Xstderr = []; end;
  if (nargin < 5) Ystderr = []; end;

  get_power = 0;
  param_boot = 0;
  power = [];

  if (nargout > 3)
    get_power = 1;
  end;
  if (xor(~isempty(Xstderr), ~isempty(Ystderr)))
    disp( '  MANTEL: if one matrix of standard errors is provided,');
    error('          both must be provided.');
  end;
  if (~isempty(Xstderr))
    param_boot = 1;
  end;

  if (isempty(iter))
    iter = 0;
  end;
  itersave = iter;

  [n,c] = size(X);
  if (~issqsym(X))
    error('  MANTEL: matrices must be square, symmetric, and of same size.');
  end;
  if (get_power)
    if (param_boot)
      if (~issqsym(Xstderr) | ~issqsym(Ystderr))
        error('  Mantel: matrices must be square, symmetric, and of same size.');
      end;
    end;
    if (~iter)
      error('  Mantel: number of iterations must be provided for power estimation.');
    end;
  end;
  if (sum(diag(X)) | sum(diag(Y)))
    disp('  MANTEL warning: diagonal values ignored in test');
  end;

  % Observed test statistics

  r = corr(trilow(X),trilow(Y));
  Z = sum(sum(X.*Y));

  if (~iter)                          % Asymptotic probability
    sumX = sum(X);
    Ax = sum(sumX);
    Bx = sum(sum(X.*X));
    Dx = sum(sumX.*sumX);
    Gx = Ax.*Ax;
    Hx = Dx-Bx;
    Kx = Gx + 2*Bx - 4*Dx;

    sumY = sum(Y);
    Ay = sum(sumY);
    By = sum(sum(Y.*Y));
    Dy = sum(sumY.*sumY);
    Gy = Ay.*Ay;
    Hy = Dy-By;
    Ky = Gy + 2*By - 4*Dy;

    L = 2*Bx*By;
    O = 4*Hx*Hy/(n-2);
    P = Kx*Ky/((n-2)*(n-3));
    Q = Gx*Gy/(n*(n-1));
    R = L+O+P-Q;

    expZ = Ax*Ay/(n*(n-1));
    seZ = sqrt(R/(n*(n-1)));
    z = (Z-expZ)/seZ;
    pr = 1-normcdf(z);
  end;

  if (iter)                           % Randomized probability
    nperm = prod(1:n);
    if (nperm <= iter)                % Use all possible permutations
      disp(sprintf('  Mantel: using all %d possible permutations',nperm));
      iter = nperm;
      Zp = zeros(nperm,1);
      Zp(1) = Z;
      i = 1:n;

      for it = 2:nperm
        i = permnext(i);
        Xp = X(i,i);
        Zp(it) = sum(sum(Xp.*Y));
      end;

    else                                % Use random permutations
      disp(sprintf('  Note: there are %d possible permutations',nperm));
      Zp = zeros(iter,1);               % Allocate matrix for permutation 
      Zp(1) = Z;                        %   test statistic values

      for it = 2:iter
        i = randperm(n);
        Xp = X(i,i);
        Zp(it) = sum(sum(Xp.*Y));
      end;
    end;

    Zp = sort(Zp);
    pr = randprob(Z,Zp);

    if (get_power)                      % Power estimation
      iter = itersave;
      Zcrit = prctiles(Zp,95);            % Critical value from null distribution
      Zp = zeros(iter,1);                 % Reallocate matrix

      for it = 1:iter                     % Randomize
        if (param_boot)                     % Parametric bootstrap
          Xb = mantelb(X,Xstderr);            
          Yb = mantelb(Y,Ystderr);
        else                                % Nonparametric bootstrap
          b = ceil(n*rand(n,1));
          Xb = X(b,b);
          Yb = Y(b,b);
        end;
        Zp(it) = sum(sum(Xb.*Yb));
      end;
      power = sum(Zp >= Zcrit)/iter;
    end;
  end;

  return;


