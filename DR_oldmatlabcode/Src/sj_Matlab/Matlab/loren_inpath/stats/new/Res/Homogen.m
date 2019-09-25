% HOMOGEN:  Measures the degree of homogeneity of a data matrix as the 
%           complement of the mean clustering level from a UPGMA of 
%           1-abs(R), where R is the correlation matrix.
%
%     Usage:[h,m,e] = homogen(X,{noplot})
%
%           X =       [n x p] data matrix.
%           noplot =  boolean variable indicating that all plots 
%                     should be suppressed [default = 0].
%           ----------------------------------------------------
%           h =       mean UPGMA clustering level.
%           m =       mean absolute correlation.
%           e =       proportion variance of first eigenvalue.
%

% RE Strauss, 10/4/99
%   5/12/00 - make corr plot conditional on variation in corr(X).

function [h,m,e] = homogen(X,noplot)
  if (nargin < 2) noplot = []; end;

  get_evals = 0;
  if (nargout > 2)
    get_evals = 1;
  end;

  if (isempty(noplot))
    noplot = 0;
  end;

  R = corr(X);
  dR = abs(R);
  t = trilow(dR);
  m = mean(t);

  if (get_evals)
    evals = eig(R);
    e = max(evals)/sum(evals);
  end;

  d = 1-dR;
  d = trisqmat(trilow(d));
  topo = upgma(d,[],noplot);
  h = 1 - mean(topo(:,4));

  if (~noplot)
    putxbnd(0,1);
    puttitle('UPGMA of variables');

    if (var(t)>1e-6)
      figure;
      histgram(t);
      putxbnd(0,1);
      puttitle('Mean Absolute Correlations');
    end;
  end;

  return;
