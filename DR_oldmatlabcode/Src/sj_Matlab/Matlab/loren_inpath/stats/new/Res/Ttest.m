% TTEST:  Pairwise t-tests for all possible pairs of groups, with sequential 
%         Bonferroni corrections for significance.
%
%     Usage: [t,pr,df,signif,pairs] = ttest(x,g,{tail},{alpha})
%
%         x =       [n x 1] observations for a single variable.
%         g =       [n x 1] classification variable for k groups.
%         tail =    flag indicating the "direction" of the test for groups:
%                     -1: left-tailed  (H0: g1 >= g2)
%                      0: 2-tailed     (H0: g1 =  g2)  [default]
%                     +1: right-tailed (H0: g1 <= g2)
%                   where g1 has the smaller group id and g2 the larger.
%         alpha =   optional criterion for statistical significance 
%                     [default=0.05].
%         -------------------------------------------------------------------
%         t =       vector of t values for k(k-1)/2 pairwise tests.
%         pr =      corresponding vector of probabilities.
%         df =      corresponding vector of degrees of freedom.
%         signif =  corresponding boolean vector indicating statistical 
%                     significance, with sequential Bonferroni adjustment.
%         pairs =   [k(k-1)/2 x 2] matrix of group identifiers for tests.
%

% RE Strauss, 12/22/98
%   6/11/99 - allow for tailed tests.
%   1/26/00 - correct the incorrect two-tailed test.

function [t,pr,df,signif,pairs] = ttest(x,g,tail,alpha)
  if (nargin < 3) tail = []; end;
  if (nargin < 4) alpha = []; end;

  [nx,px] = size(x);
  [ng,pg] = size(g);
  err = 0;
  if (min([nx,px])>1 | min([ng,pg])>1)
    err = 1;
  end;
  if (length(x) ~= length(g))
    err = 1;
  end;
  if (err)
    error('TTEST: input data and group matrix must be compatible vectors');
  end;

  if (isempty(tail))
    tail = 0;
  end;
  if (isempty(alpha))
    alpha = 0.05;
  end;

  grpid = uniquef(g,1);               % Unique group identifiers, in sequence
  k = length(grpid);                  % Number of groups
  npairs = k*(k-1)/2;                 % Number of pairs of groups

  if (k<2)
    error('  TTEST: need at least two groups');
  end;

  t =     zeros(npairs,1);            % Allocate output matrices
  df =    zeros(npairs,1);
  pairs = zeros(npairs,2);

  test = 0;
  for i1 = 1:(k-1)                    % Cycle thru pairs of groups
    g1 = grpid(i1);
    for i2 = (i1+1):k
      test = test+1;
      g2 = grpid(i2);
      pairs(test,:) = [g1 g2];
      i = find(g==g1 | g==g2);
      xx = x(i);
      gg = g(i);
      t(test) = tval(xx,gg);
      df(test) = length(i)-2;
    end;
  end;

  if (tail==0)
    pr = 2*tcdf(-abs(t),df);
  elseif (tail>0)
    pr = 1-tcdf(t,df);
  else
    pr = tcdf(t,df);
  end;

  if (k==2)
    signif = 0;
    if (pr <= alpha)
      signif = 1;
    end;
  else
    signif = seqbonf(pr,alpha);
  end;

  return;




