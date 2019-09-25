% WRIGHTCI: Determines confidence intervals about factor loadings, R-squared 
%           values, and residual correlations for Wright's path analysis, either 
%           by bootstrapping the original character matrix or randomizing a 
%           correlation or covariance matrix.  See documentation for function WRIGHT.
%
%     Syntax: [f,fl,fh,r2,r2l,r2h,rcf,rcfl,rcfh,rcp,rcpl,rcph]
%             = wrightci(X|C,{secnd},{covar}|N,{iter},{orthog},{primary},{CI_level})
%
%        |X =       [n x p] data matrix
%        |  OR
%        |C =       [p x p] covariance or correlation matrix.
%         secnd =   [n x s] matrix of boolean values (T/F = 1/0)
%                     indicating by 1's the submatrices of variables
%                     to be included in the secondary factors (s of them).
%                     Can be passed as a null matrix.
%        |covar =   boolean flag indicating whether the covariance
%        |            matrix (T=1) or correlation matrix (F=0) is to be factored
%        |            from the data matrix [default = covariance].
%        |            Can be passed as a null matrix if cov/corr matrix is supplied.
%           OR
%        |N =       sample size, if cov/corr matrix is supplied.
%         iter =    number of bootstrap iterations (default = 1000).
%         orthog =  optional flag indicating whether secondary factors are to be:
%                     0 - oblique to primary factor and one another;
%                     1 - orthogonal to primary factor but oblique to other
%                           secondary factors;
%                     2 - orthogonal to primary factor and all other secondary
%                           factors;
%                     [default = 2].
%         primary =  boolean flag indicating whether, in the presence of
%                       secondary factors, the first factor is to be primary (T=1)
%                       or general (F=0) [default = primary].
%         CI_level = percent width of confidence intervals (default = 95)
%         ----------------------------------------------------------------------
%         f =   [n x s+1] matrix of factor loadings (primary + secondary).
%         fl =  corresponding low confidence limits.
%         fh =  corresponding high confidence limits.
%         r2 =  proportion of total sums-of-squares accounted for by model(s).
%         r2l = corresponding low confidence limits.
%         r2h = corresponding high confidence limits.
%         rcf =  [n x n] matrix of residual covariances, full model.
%         rcfl = corresponding low confidence limits.
%         rcfh = corresponding high confidence limits.
%         rcp =  [n x n] matrix of residual covariances, primary factor only.
%         rcpl = corresponding low confidence limits.
%         rcph = corresponding high confidence limits.
%

function [f,fl,fh,r2,r2l,r2h,rcf,rcfl,rcfh,rcp,rcpl,rcph] ...
            = wrightci(X,secnd,covar,iter,orthog,primary,CI_level)

  TRUE = 1; FALSE = 0;

  % Determine whether data matrix or cov/corr matrix has been passed

  [N,P] = size(X);
  is_data = TRUE;
  if (N==P)
    if (sum(abs(X-X'))==0)
      is_data = FALSE;
      C = X;
      if (nargin < 4)
        error('  Error: sample size required if cov/corr matrix is to be randomized')
      else
        N = covar;
      end;
    end;
  end;
  
  % Initialize arguments and flags

  default_secnd = [];
  default_iter = 1000;
  default_iter = 0;
  default_covar = TRUE;
  default_orthog = 2;
  default_primary = TRUE;
  default_CI_level = 95;

  if (nargin < 2)
    secnd = default_secnd;
  end;
  if (nargin < 3)
    iter = default_iter;
  end;
  if (nargin < 4)
    covar = default_covar;
  end;
  if (nargin < 5)
    orthog = default_orthog;
  end;
  if (nargin < 6)
    primary = default_primary;
  end;
  if (nargin < 7)
    CI_level = default_CI_level;
  end;

  if (nargout > 3)
    boot_r2 = TRUE;
  else
    boot_r2 = FALSE;
  end;
  if (nargout > 6)
    boot_rcf = TRUE;
  else
    boot_rcf = FALSE;
  end;
  if (nargout > 9)
    boot_rcp = TRUE;
  else
    boot_rcp = FALSE;
  end;

  fl = [];                            % Allocate output arguments
  fh = [];
  r2l = [];
  r2h = [];
  rcpl = [];
  rcph = [];
  rcfl = [];
  rcfh = [];

  % Covariances or correlations from data matrix

  if (is_data)
    if (covar)
      C = cov(X);
    else
      C = corrcoef(X);
    end;
  end;

  % Get point estimates for factor-analysis parameters

  [f,r2,rcf,rcp] = wright(C,secnd,orthog,primary);
f
r2
rcf
rcp

  % Bootstrap or randomize the factor analysis

  if (iter > 0)
    [fr,fc] = size(f);                % Dimensions of factor matrix
    Bf = zeros(iter,fr*fc);           % Allocate working bootstrap matrices
    if (boot_r2)
      Br2 = zeros(iter,length(r2));
    end;
    if (boot_rcp)
      Brcp = zeros(iter,P*P);
    end;
    if (boot_rcf)
      Brcf = zeros(iter,P*P);
    end;

    clk = clock;                      % Seed random-number generator
    rand('seed',clk(6));              % from system clock

    Bf(1,:) = f(:)';                  % Begin with initial solution, transposed
    if (boot_r2)
      Br2(1,:) = r2;
    end;
    if (boot_rcf)
      Brcf(1,:) = rcf(:)';            % Store corr matrix as transposed column vector
    end;
    if (boot_rcp & fc>1)
      Brcp(1,:) = rcp(:)';            % Store corr matrix as transposed column vector
    end;

    for bit = 2:iter                  % Bootstrap iterations
%if (rem(bit,50)==0)
  bit
%end;
      if (is_data)                      % If data matrix supplied,
        BX = bootsamp(X);                 % Sample rows of X
        if (covar)                        % Covar/corr matrix
          BC = cov(BX);
        else
          BC = corrcoef(BX);
        end;
      else                              % If cov/corr matrix supplied,
        BC = randcorr(C,N);
      end;

      [bf,br2,brcf,brcp] = wright(BC,secnd,orthog,primary,1); % Factor analysis

      vc = diag(vectcorr(f,bf));        % Check directions of factors
      if (any(vc<0))                      % Reverse if negative of original
        i = find(vc<0);
        bf(:,i) = -bf(:,i);
      end;

      Bf(bit,:) = bf(:)';               % Stash results
      if (boot_r2)
        Br2(bit,:) = br2;
      end;
      if (boot_rcf)
        Brcf(bit,:) = brcf(:)';
      end;
      if (boot_rcp & fc>1)
        Brcp(bit,:) = brcp(:)';
      end;
    end;  % End bootstrap iterations

    high = ceil(iter*CI_level/100);   % CI indices
    low = iter - high + 1;

    Bf = sort(Bf);                    % Conf limits for factor loadings
    fl = reshape(Bf(low,:)',fr,fc);
    fh = reshape(Bf(high,:)',fr,fc);

    if (boot_r2)                      % Conf limits for proportion fit
      Br2 = sort(Br2);
      r2l = Br2(low,:);
      r2h = Br2(high,:);
    end;

    if (boot_rcp)                     % Conf limits for resid covars, primary
      Brcp = sort(Brcp);
      rcpl = reshape(Brcp(low,:)',P,P);
      rcph = reshape(Brcp(high,:)',P,P);
    end;

    if (boot_rcf)                     % Conf limits for resid covars, full
      Brcf = sort(Brcf);
      rcfl = reshape(Brcf(low,:)',P,P);
      rcfh = reshape(Brcf(high,:)',P,P);
    end;

    disp(sprintf(' Bootstrap time: %5.2f minutes\n', ...
          etime(clock,clk)/60));
  end;  % End bootstrap

  return;
