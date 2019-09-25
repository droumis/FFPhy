% MANNWHIT: Mann-Whitney-Wilcoxon 2-sample rank-sum test for group differences.  
%           Groups are identified by the collating sequence of the 
%           classification variable.
%
%   Syntax: [pr,U,N,R] = mannwhit(x,grps,{tail})
%
%         x =     [n x 1] observations for a single variable.
%         grps =  [n x 1] classification variable.
%         tail =  flag indicating the "direction" of the test for groups g1 & g2:
%                   -1: left-tailed  (H0: g1 <= g2)
%                    0: 2-tailed     (H0: g1 =  g2)  [default]
%                   +1: right-tailed (H0: g1 >= g2)
%         --------------------------------------------------------------------------
%         pr =    significance level of the test.
%         U =     observed test-statistic values for the two groups.
%         N =     sample sizes for the two groups.
%         R =     sums of ranks for the two groups.
%

% Exact probabilities based on Applied Statistics algorithm
%   AS 62, Appl. Statist. 22(2), 1973.

% RE Strauss, 11/3/97
%   4/26/99 - exact probabilities implemented.
%   4/10/01 - miscellaneous improvements.
%   5/5/03 - corrected error in calculation of test-statistic value.

function [pr,U,N,R] = mannwhit(x,grps,tail)
  if (~nargin) help mannwhit; return; end;

  if (nargin < 3) tail = []; end;

  [n,p] = size(x);
  if (n==1)                                 % Input must be column vector
    x = x';
    [n,p] = size(x);
  end;
  if (p>1)
    error(' MANNWHIT: dependent variable must be a vector');
  end;

  if (isempty(tail))                        % Default input arguments
    tail = 0;
  end;

  [grpid,N] = uniquef(grps,1);              % Get grp id's and sample sizes
  if (length(grpid)~=2)
    error('  MANNWHIT: 2-group test only.');
  end;

  R = ranks(x);                             % Ranks of pooled data
  R1 = sum(R(grps==grpid(1)));              % Sums of ranks
  R2 = sum(R(grps==grpid(2)));
  R = [R1;R2];

  [Nmax,imax] = max(N);
  [Nmin,imin] = min(N);
  
  U1 = prod(N) + Nmax*(Nmax+1)/2 - R(imax); % Mann-Whitney U statistics
  U2 = prod(N) - U1;
  U = [U1;U2];

  pr = mannprob(U,N,tail);                % Probability level
      
  return;
