% WRIGHT: Wright-style factor analysis, based on Wright (1954).
%         Estimate a general/primary factor plus one or more orthogonal
%         secondary factors from a covariance (or correlation) matrix,
%         using a least-squares minimization of residual covariances.
%         Alternately, a general/primary factor plus at most one oblique
%         secondary factor can be estimated.
%
%     Syntax: [f,r2a,rc_full,rc_primary]
%                 = wright(c,{secnd},{orthog},{primary},{suppress})
%
%         c =       [p x p] covariance or correlation matrix.
%         secnd =     [n x s] matrix of boolean values (T/F = 1/0)
%                     indicating by 1's the submatrices of variables
%                     to be included in the secondary factors (s of them).
%         orthog =  optional flag indicating whether secondary factors are to be:
%                     0 - oblique to primary factor and one another;
%                     1 - orthogonal to primary factor but oblique to other
%                           secondary factors;
%                     2 - orthogonal to primary factor and all other secondary
%                           factors;
%                     [default = 0].
%         primary = optional boolean flag indicating whether first factor is
%                     to be a primary factor (ignoring covariances sequestered
%                     into secondary factor) or a general factor (fitted using
%                     all covariances).
%                     [default = primary if secondary factors are specified,
%                     general otherwise].
%         suppress = optional boolean flag indicating whether warning messages
%                     are to be suppressed (invoked for bootstrapping and
%                     randomization).
%         ----------------------------------------------------------------------
%         f =          [p x s+1] matrix of factor loadings (primary + secondary).
%         r2a =        adjusted proportion of total sums-of-squares accounted for by
%                        model; if secondary factor(s) specified, a 2-element row
%                        vector is returned, where the first element is the R2_a for
%                        the primary factor and the second element is the R2_a for
%                        the primary + secondary factors (full model).
%         rc_full =    [p x p] matrix of residual covariances for primary + secondary
%                        factors (full model).
%         rc_primary = [p x p] matrix of residual covariances for primary
%                        factor only.
%

% RE Strauss
%   9/1/98 -  changed R2 to adjusted R2 (Zar 1996:417), accounting for number of parameters
%               in the model.
%          -  if secondary factor passed as a row vector, transpose automatically to col vector.

function [f,r2a,rc_full,rc_primary] = wright(c,secnd,orthog,primary,suppress)
  if (nargin < 2) secnd = []; end;
  if (nargin < 3) orthog = []; end;
  if (nargin < 4) primary = []; end;
  if (nargin < 5) suppress = []; end;

  [p,nvar] = size(c);                 % Number of variables

  % Check for valid input matrix

  is_cov = 1;
  if (p ~= nvar)
    is_cov = 0;
  else
    if (sum(sum(abs(c-c'))))
      is_cov = 0;
    end;
  end;
  if (~is_cov)
    error('    Error: input must be covariance or correlation matrix.');
  end;

  % Default input arguments

  if (isempty(suppress))
    suppress = 0;
  end;
  if (isempty(primary))
    primary = 1;
  end;
  if (isempty(orthog))
    orthog = 0;
  end;

  % Check for valid secondary-factor specifications

  if (isempty(secnd))
    nsec = 0;
  else
    [p,nsec] = size(secnd);
    if (p==1 & nsec>1)                % Transpose row vector to column vector
      secnd = secnd';
      [p,nsec] = size(secnd);
    end;

    if (p>0 & p~=nvar)
      error('    Error: secondary-factor matrix must have P rows.');
    end;

    ssec = sum(secnd);
    if (any(ssec < 2))
      error('    Error: secondary factors must have at least 2 variables.');
    end;
  end;

  % Estimate factor loadings and residual matrices

%  if (~suppress)
%    disp('  Fitting primary factor...');
%  end;
  [f,rc_primary] = factorp(c,secnd,primary);  % Estimate general or primary factor

  if (nsec > 0)
%    if (~suppress)
%      disp('  Fitting secondary factors...');
%    end;
    [f,rc_full] = factors(rc_primary,f,secnd,orthog,suppress); % Est secondary factors

%    if (~suppress)
%      disp('  Optimizing full model...');
%    end;
%    [f,rc_full] = factorw(c,f,secnd,orthog);      % Final global optimization
  else
    rc_full = rc_primary;
  end;


  % Estimate adjusted R-squared values for primary and full models (Zar 1996:417)

  n = p*(p-1)/2;
  m = p;

  sse =  sum(trilow(rc_primary).^2);    % SSE of covars, primary factor only
  ssto = sum(trilow(c).^2);             % SSTO of covars
  r2 = 1-sse/ssto;                      % Proportion SS accounted for
  r2a = 1 - ((1-r2)*(n-1)/(n-m-1));     % Adjust for number of parameters

  if (nsec > 0)
    m = m + sum(sum(secnd));

    sse =  sum(trilow(rc_full).^2);       % SSE of covars, full model
    r2 = 1-sse/ssto;                      % Proportion SS accounted for, primary
    r2 = 1 - ((1-r2)*(n-1)/(n-m-1));      % Adjust for number of parameters

    r2a = [r2a r2];
  end;

  % Check orthogonality of secondary factors

  if (nsec>0 & ~suppress)
    p = f(:,1);
    s = f(:,2:(nsec+1));

    if (orthog>0)
      vc = vectcorr(p,s);
      if (sum(vc)>eps)
        disp('    Warning: not all secondary factors are orthogonal');
        disp('             to primary factor.');
        disp(' ');
      end;
    end;

    if (nsec>1)
      vc = vectcorr(s);
      if (sum(trilow(vc))>eps)
        if (orthog==2)
          disp('    Warning: not all secondary factors are orthogonal');
          disp('             to one another.');
          disp(' ');
        end;
        if (max(sum(secnd'))>1);
          disp('    Warning: secondary residual covariances may be underestimated');
          disp('             when secondary factors are overlapping and oblique.');
          disp(' ');
        end;
      end;
    end;
  end;

  return;
