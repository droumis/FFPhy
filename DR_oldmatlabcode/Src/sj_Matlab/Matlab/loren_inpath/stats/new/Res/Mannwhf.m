% MANNWHF: Objective function for mannwhit(), called by bootstrp().
%
%     Syntax: W = mannwhf(R,grps,[],[],tail)
%
%         R =     [n x 1] vector of ranks of data.
%         grps =  [n x 1] vector of group memberships.
%         grpid = vector (length 2) of unique group identifiers.
%         tail =  tail-direction of test.
%         ---------------------------------------------------------
%         W =     Wilcoxon statistic.
%

% RE Strauss, 10/30/97
%   1/28/99 - modified to calc Wilcoxon W rather than Mann-Whitney U,
%               to allow for negative z-scores

function W = mannwhf(R,grps,not_used1,not_used2,tail)
  [grpid,n] = uniquef(grps);

  if (tail>1)                         % Right-tailed test (g1 < g2)
    W = sum(R(grps==grpid(2)));         % Rank sum for group 2

  elseif (tail<1)                     % Left-tailed test  (g2 < g1)
    W = sum(R) - sum(R(grps==grpid(1)));  % Complement of rank sum for group 1

  else                                % Two-tailed test   (g1 < g2 or g2 > g1)
    w1 = sum(R(grps==grpid(2)));
    w2 = sum(R) - sum(R(grps==grpid(1)));
    W = max([w1 w2]);
  end;

  return;


