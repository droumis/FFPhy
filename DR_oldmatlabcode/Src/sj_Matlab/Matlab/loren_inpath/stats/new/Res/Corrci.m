% CORRCI: Estimates the simultaneous (Bonferroni) upper and lower confidence
%         limits for a correlation or covariance matrix, either asymptotically,
%         using Fisher's z-transform. or via randomization.  For a single 
%         correlation coefficient, returns the asymptotic confidence limits.
%
%         If a covariance matrix is supplied, confidence limits are estimated
%         directly by randomization, but indirectly (via the correlation matrix)
%         for asymptotic estimates.
%
%     Syntax: [low,high] = corrci(R,N,{iter},{CI_level})
%
%         R =     [P x P] matrix of correlation coefficients.
%         N =     sample size (scalar).
%         iter =  number of iterations for randomized intervals,
%                   or zero for asymptotic estimates [default = 0].
%         CI_level = optional significance level to determine critical value
%                   [default = 0.95].
%         ------------------------------------------------------------------
%         low =   corresponding lower confidence limits
%         high =  corresponding upper confidence limits
%

% RE Strauss, 1/23/96
%   12/27/99 - changed handling of input arguments;
%              added option of single correlation coefficient.

function [low,high] = corrci(R,N,iter,CI_level)
  if (nargin < 2) N = []; end;
  if (nargin < 3) iter = []; end;
  if (nargin < 4) CI_level = []; end;

  if (isempty(N))
    error('  CORRCI: sample size required');
  end;
  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(CI_level))
    CI_level = 0.95;
  end;

  if (CI_level > 1)
    CI_level = CI_level/100;
  end;
  alpha = 1-CI_level;

  iter_display = 500;                   % Iteration display interval

  [r,P] = size(R);
  if (r==1 & P==1)                      % Single correlation coefficient
    if ((R<-1) | (R>1))  
      error('  CORRCI: correlation coefficient out of range');
    end;    
    
    [z,stderr] = corrz(R,N);            % Fisher's z
    low = norminv(alpha/2,z,stderr);    % Interval on z
    high = norminv(1-alpha/2,z,stderr);
    low =  tanh(low);                   % Retransform to correlation
    high = tanh(high);

    return;
  end;

  is_corr = 1;
  is_cov =  1;
  if (r~=P)                             % Check for square matrix
    is_corr = 0;
    is_cov =  0;
  else
    if (sum(sum(abs(R-R'))) > eps)      % Check for symmetric matrix
      is_corr = 0;
      is_cov =  0;
    else
      if (any(abs(diag(R)-ones(P,1))) > eps) % Check for diagonal of 1s
        is_corr = 0;
      else
        is_cov =  0;
      end;
    end;
  end;
  if (~is_corr & ~is_cov)
    error('  Error: invalid correlation/covariance matrix');
  end;

  ncorr = P*(P-1)/2;                    % Number of pairwise correlations
  alpha = alpha / ncorr;                % Bonferroni adjustment of signif level

  % Asymptotic estimates

  if (iter == 0)
    if (is_cov)                         % Convert covariances to correlations
      [R,S] = covcorr(R);
    end;

    [z,stderr] = corrz(R,N);            % Fisher's z
    low = norminv(alpha/2,z,stderr);    % Interval on z
    high = norminv(1-alpha/2,z,stderr);
    low =  tanh(low);                   % Retransform to correlations
    high = tanh(high);

    if (is_cov)
      low = covcorr(low,S);             % Rescale correlations to covariances
      high = covcorr(high,S);
    end;
  end;

  % Randomized estimates

  if (iter > 0)
    Br = zeros(iter,ncorr);             % Allocate randomization matrix
    Br(1,:) = trilow(R)';               % Stash input matrix
    if (is_cov)                         % For covariances,
      Bd = zeros(iter,P);               %   also allocate diagonal rand matrix
      Bd(1,:) = diag(R)';               %   and stash input
    end;

    mu = zeros(1,P);                    % Params of multinormal sampling distrib
    if (is_corr)
      S =  ones(1,P);
    end;

    for it = 2:iter                     % Randomization iterations
      if (rem(it,iter_display)==0)
        disp(sprintf('    Iterations = %d',it));
      end;

      if (is_corr)                        % For correlations,
        x = randmvn(N,mu,S,R);              % standardized multinormal sample
        r = corrcoef(x);                    % random corr matrix
      else                                % For covariances,
        x = randmvn(N,mu,R);                % multinormal sample
        r = cov(x);                         % random covariance matrix
        Bd(it,:) = diag(r)';                % stash diagonal
      end;
      Br(it,:) = trilow(r)';              % Stash correlations
    end;

    Br = sort(Br);                      % Sort all accumulated corrs
    limits =  bootci(Br,alpha);         % Retrieve confidence limits

    if (is_corr)                        % For correlations,
      low =  trisqmat(limits(1,:)',1);    % stash into square matrices
      high = trisqmat(limits(3,:)',1);
    else                                % For covariances,
      Bd = sort(Bd);                      % sort diagonal elements
      ld = bootci(Bd,alpha);              % retrieve confidence limits
      low =  trisqmat(limits(1,:)',ld(1,:)'); % stash into square matrices
      high = trisqmat(limits(3,:)',ld(3,:)');
    end;
  end;

  return;
