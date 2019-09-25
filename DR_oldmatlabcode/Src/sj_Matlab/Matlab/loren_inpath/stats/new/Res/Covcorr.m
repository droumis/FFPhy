% COVCORR: Rescales a covariance matrix to a correlation matrix + vector of
%          standard deviations.  Or, rescales correlation matrix + vector of
%          standard deviations to a covariance matrix.
%
%     Syntax: [R,S] = covcorr(C)
%               OR
%             C = covcorr(R,S)
%
%          C = covariance matrix
%          R = correlation matrix
%          S = column vector of standard deviations
%

% RE Strauss, 1/13/96
%   5/17/00 - use iscorr() to check for valid matrix form;
%             check and disallow missing data.

function [Mo,So] = covcorr(Mi,Si)
  if (~isfinite(sum(sum(Mi))))
    error('  COVCORR: matrix contains missing data.');
  end;

  % Covariances to correlations

  if (nargin == 1)
    [r,c] = size(Mi);
    dc = diag(Mi);

    if (~iscov(Mi))
      error('  COVCORR: invalid covariance matrix');
    end;

    So = sqrt(dc);
    ds = diag(1./So);
    Mo = ds * Mi * ds;
    Mo = Mo - 0.5*(Mo-Mo');           % Correct any round-off error
    Mo = Mo - diag(diag(Mo)-ones(r,1));
  end;

  % Correlations to covariances

  if (nargin == 2)
    [r,c] = size(Mi);
    dc = diag(Mi);
    Si = diag(diag(Si));              % Convert to column vector

    if (~iscorr(Mi))
      error('  COVCORR: invalid correlation matrix');
    end;

    ds = diag(Si);
    Mo = ds * Mi * ds;
    Mo = Mo - 0.5*(Mo-Mo');           % Correct any round-off error
    Mo = Mo - diag(diag(Mo)-Si.*Si);
  end;

  return;
