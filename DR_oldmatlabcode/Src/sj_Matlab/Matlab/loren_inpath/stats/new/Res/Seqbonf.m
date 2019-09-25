% SEQBONF: Sequential Bonferroni test, given the probabilities from a set 
%          of k previously executed tests.
%
%     Syntax: [signif,critlim] = seqbonf(pvalues,{alpha},{independent})
%
%           pvalues =     [k x p] matrix of probabilities (p-values) from 
%                           k tests, evaluated separately by column, OR
%                         [k x k] square symmetric matrix of probabilities, 
%                           of which the lower triangular portion is evaluated.
%           alpha =       desired family-wide critical value [default = 0.05].
%           independent = boolean variables indicating whether (=1) or 
%                           not (=0) the set of tests can be considered to be 
%                           independent (producing a slight increase in power) 
%                           [default = 0].
%           -------------------------------------------------------------------
%           signif =      boolean vector (corresponding to 'pvalues') 
%                           indicating tests to be considered significant.
%           critlim =     critical limits corresponding to each observed 
%                           probability level.
%

% References:
%   Holm, S.  1979.  A simple sequentially rejective multiple test 
%     procedure.  Scand. J. Stat. 6:65-70.
%   Miller, R.G. Jr.  1981.  Simultaneous Statistical Inference.  McGraw Hill.
%   Rice, W.R. 1989. Analyzing tables of statistical tests. Evolution 43:223-225.

% RE Strauss, 12/16/98
%    9/19/99 - misc changes for Matlab v5.
%   12/25/99 - output the critical limits.

function [signif,critlim] = seqbonf(pvalues,alpha,independent)
  if (nargin < 2) alpha = []; end;
  if (nargin < 3) independent = []; end;

  [k,p] = size(pvalues);

  row_vector = 0;
  sqsymmat = 0;

  if (isempty(alpha))                 % Default input arguments
    alpha = 0.05;
  end;
  if (isempty(independent))
    independent = 0;
  end;

  if (alpha > 1)                      % Convert percentage to proportion
    alpha = alpha / 100;
  end;

  if (k==1 & p>1)                     % Transpose row vector
    pvalues = pvalues';
    row_vector = 1;
    k = p;
    p = 1;
  elseif (k==p)                       % Extract lower triangular portion
    trilowp = trilow(pvalues);        %   of square symmetric matrix
    if (trilowp == trilow(pvalues'))
      pvalues = trilowp;
      sqsymmat = 1;
      k = length(pvalues);
      p = 1;
    end;
  end;

  if (independent)                    % Critical limits
    critlim = 1 - (1-alpha).^(1./(1+k-([1:k]')));
  else
    critlim = alpha ./ (1+k-([1:k]'));  
  end;

  [pvalues,indx] = sort(pvalues);     % Sort p-values, keeping track of indices
  signif = zeros(size(pvalues));

  for c = 1:p                         % Cycle through columns
    signif(:,c) = (pvalues(:,c) <= critlim); % Boolean significance decision
    s = signif(:,c);                  % All p-values after the first
    found_insignif = 0;               %   non-significant one are non-significant
    for i = 1:k
      if (~signif(i,c))
        found_insignif = 1;
      elseif (found_insignif)
        signif(i,c) = 0;
      end;
    end;

    for i = 1:k                         % Resort into original sequence
      s(indx(i,c)) = signif(i,c);
    end;
    signif(:,c) = s;
  end;

  pvalues(indx) = pvalues;
  critlim(indx) = critlim;

  if (row_vector)                     % Transpose for row vector
    signif = signif';
  end;

  if (sqsymmat)                       % Reconstruct square symmetric matrix
    signif = trisqmat(signif);
    critlim = trisqmat(critlim);
  end;

  return;

