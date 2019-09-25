% CORRPROB: Estimates the probability that the entire correlation matrix
%           does not differ from a diagonal matrix; equivalent to testing
%           H0:variables are independent.  The test is done either asymptotically
%           using Lawley's test (Morrison 1976:116-119) or via randomization,
%           for which the test statistic is the sum of squared correlations.
%
%           Estimates the 2-tailed significance of a Pearson product-moment
%           correlation coefficient (against H0:corr=0), either asymptotically,
%           using Fisher's z-transform (Sokal & Rohlf 1981:584-587), or
%           via randomization of the correlation matrix.
%             In either case, a critical absolute correlation is returned for
%           Bonferroni-adjusted decisions as to which individual correlation
%           coefficients are significant, and the correlations greater in
%           magnitude than this critical value are returned in decreasing order
%           of magnitude.
%
%           If a covariance matrix (square symmetric, diagonals not all unity)
%           is given, it is converted to a correlation matrix for any
%           asymptotic tests or for finding probabilities of correlation
%           coefficients via randomization.
%             However, if only the diagonal test is requested via randomization,
%           the covariance matrix is used directly.  This allows the testing
%           of matrices having zeros on the diagonal (eg, residual covariances
%           from WRIGHT() ).
%
%     Syntax: [matprob,prob,critcorr,sigcorr] = corrprob(R,N,{iter},{alpha})
%
%           R =    [P x P] matrix of correlation coefficients.
%           N =    sample size (scalar).
%           iter = number of iterations for randomized probabilities, or zero
%                    for asymptotic estimates [default = 0].
%           alpha = optional significance level to determine critical value
%                    [default = 0.05].
%           ------------------------------------------------------------------
%           matprob =  1-tailed probability level against H0:variables are
%                        independent.
%           prob =     [p x p] matrix of simultaneous 2-tailed probability levels
%                        against H0:corr=0.
%           critcorr = critical multiple-comparisons magnitude for testing
%                        individual correlation coefficients.
%           sigcorr =  3-column matrix of significant correlations (based on
%                        'critcorr', sorted by decreasing magnitude.
%                           Col 1 = correlation coefficient
%                               2 = row subscript
%                               3 = column subscript
%

% RE Strauss, 2/16/97

function [matprob,prob,critcorr,sigcorr] = corrprob(R,N,iter,alpha)
  if (nargin < 3) iter = []; end;
  if (nargin < 4) iter = []; end;

  iter_display = 500;               % Iteration display interval

  if (isempty(iter))                % Default arguments
    iter = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  else
    if (alpha > 1)
      alpha = alpha / 100;
    end;
  end;

  get_prob = 0;
  get_critcorr = 0;
  get_sigcorr = 0;

  if (nargout > 1)
    get_prob = 1;
  end;
  if (nargout > 2)
    get_critcorr = 1;
  end;
  if (nargout > 3)
    get_sigcorr = 1;
  end;

  [r,P] = size(R);

  if (r==1 & P==1)                    % Convert scalar to 2x2 matrix
    R = [1 R;R 1];
    [r,P] = size(R);
  end;

  if (r ~= P)                         % Check for square matrix
    error('  CORRPROB: matrix not square');
  end;
  if (sum(sum(abs(R-R'))) > eps)      % Check for symmetric matrix
    error('  CORRPROB: matrix not symmetric');
  end;
  if (~all(diag(R)>0))                % Check for positive diagonal elements
    error('  CORRPROB: matrix has zero or negative diagonal elements');
  end;

  if (any(abs(diag(R)-ones(P,1))) > eps) % If diagonal elements not all 1s,
    if (iter==0 | get_prob)              %   and corr coeff probs are requested,
      R = covcorr(R);                    %   convert covar matrix to corr matrix
    end;
  end;

  ncorr = P*(P-1)/2;                  % Number of pairwise correlations
  Ra = abs(R);                        % Convert to positive for rt-tailed probs

  % Asymptotic estimates

  if (iter == 0)
    t = trilow(R);
    ssr = t'*t;                       % Sum of squared off-diag correlations
    adj = N-1 - (2*P+5)/6;            % Adjustment for numbers of obs & vars
    L = adj * ssr;                    % Lawley's chi-square statistic
    df = ncorr;
    matprob = 1 - chi2cdf(L,df);      % Diagonal-test probability

    if (get_prob)                     % Correlation-coefficient probabilities
      [z,stderr,prob] = corrz(Ra,N);    % Fisher's z
    end;

    if (get_critcorr)                 % Critical correlation magnitude
      alpha = alpha / ncorr;          % Bonferroni-adjusted signif level
      cz = stderr * (norminv(1-alpha/2));  % Critical Fisher-z (2-tailed)
      critcorr = tanh(cz);                 % Critical absolute correlation
    end;
  end;

  % Randomized estimates

  if (iter > 0)
    Bssr = zeros(iter,1);             % Allocate randomization matrices
    if (get_prob)
      Br = zeros(iter,ncorr);
    end;

    for it = 1:iter                   % Randomization iterations
      if (rem(it,iter_display)==0)
        disp(sprintf('    Iterations = %d',it));
      end;
      r = corrcoef(randn(N,P));         % Random corr matrix
      t = trilow(abs(r));               % Lower triangular matrix
      Bssr(it) = t'*t;                  % Stash sum of squared correlations
      if (get_prob)
        Br(it,:) = t';                  % Stash absolute correlations
      end;
    end;

    t = trilow(R);                    % Prob of sum of squared correlations
    Hssr = t'*t;                        % Observed value
    Bssr = sort(Bssr);                  % Sort null distribution
    matprob = randprob(Hssr,Bssr)';     % Right-tailed probability

    if (get_prob)                     % Probs of corr coefficients
      t = trilow(Ra);                   % Use abs(R) for rt-tailed probs
      Hr = t;                           % Vector of observed values
      Br = sort(Br(:));                 % Sort all accumulated abs null corrs

      prob = zeros(P,P);                % Allocate probilities matrix
      k = 0;
      for i = 1:(P-1)
        for j = (i+1):P
          k = k+1;
          prob(i,j) = randprob(Hr(k),Br)'; % Stash probs in upper diagonal
        end;
      end;
      prob = max(prob,prob');           % Make symmetric
    end;

    if (get_critcorr)                 % Critical absolute correlation value
      alpha = alpha / ncorr;            % Bonferroni-adjusted signif level
      indx = iter*ncorr * alpha;        % Indices into sorted null distrib
      crit_high = length(Br) - floor(indx); % Interpolate critical abs corr
      crit_low  = length(Br) - ceil(indx);
      if (crit_low ~= crit_high)
        delta = (ceil(indx)-indx)/(ceil(indx)-floor(indx));
      else
        delta = 0;
      end;
      critcorr_low  = Br(crit_low);
      critcorr_high = Br(crit_high);
      critcorr = critcorr_low + delta*(critcorr_high - critcorr_low);
    end;
  end;

  % Return sorted significant correlation coefficients

  if (get_sigcorr)
    sigcorr = [];
    b = Ra > critcorr;                % Find significant absolute coefficients
    [t,j,i] = trilow(b);              % Lower triangular flags
    if (sum(t))                       % If any coefficients are significant,
      c = trilow(Ra);
      indx = find(t);                   % Isolate significant coefficients
      c = c(indx);
      i = i(indx);
      j = j(indx);
      [csort,indx] = sort(-c);          % Sort signif corrs by decreasing magnitude
      sigcorr = zeros(length(i),3);     % Allocate output matrix
      for k = 1:length(indx)            % Stash original signed correlations
        ii = i(indx(k));                %    by decreasing magnitude
        jj = j(indx(k));
        sigcorr(k,1) = R(ii,jj);
        sigcorr(k,2:3) = [ii jj];
      end;
    end;
  end;

  return;
