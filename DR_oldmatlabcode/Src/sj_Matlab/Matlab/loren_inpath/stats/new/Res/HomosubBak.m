% HOMOSUB:  Composes a matrix of homogeneous subsets (rows), coded in lower-case 
%           letters, given either (1) an NxN symmetric binary matrix indicating 
%           which of N groups are significantly different from one another, or 
%           (2) the original data, from which means and standard errors are 
%           estimated.  The first method is qualitative, based on presence/absence 
%           of significant difference, and the measure of consistency is the 
%           proportion of presence/absence significance states that is in 
%           agreement with the assignment of groups to homogeneous subsets.
%           The second method is quantitative, and the measure of consistency is 
%           the Spearman rank correlation between intergroups differences (as 
%           measured by t-statistics) and the assignment of groups to homogeneous 
%           subsets.  Subsets often comprise ambiguous, overlapping sets of groups.
%           Given are the subsets most consistent with the input matrix; all 
%           subsets returned are equally consistent.
%
%     Usage: [hs,consistency] = homosub(S)
%                   OR
%            [hs,consistency] = homosub(x,grps)
%
%           S =     [n x n] symmetric boolean matrix indicating which pairwise 
%                     comparisons are significantly different.
%                       OR
%           x =     vector (length n) of data.
%           grps =  corresponding group-membership vector.
%           ---------------------------------------------------------------------
%           hs =    [ns x n] character matrix with possible homogeneous subsets.
%           consistency = measure of total pairwise consistencies [0-1] 
%                           of subsets with S.
%

% RE Strauss, 12/25/99
%   7/24/00 - for binary comparisons, give sig difference a slight bias over 
%             non-sig difference in consistency measure.

function [hs,consistency] = homosub(x,grps)
  if (nargin < 2) grps = []; end;

  if (isempty(grps))
    use_data = 0;
    S = x;
    [n,m] = size(S);
  else
    use_data = 1;
    x = x(:);
    grps = grps(:);
    ugrps = uniquef(grps,1);
    n = length(ugrps);
  end;

  hs = [];                                % Find all potential subsets
  hs = [];
  for i = 1:n
    hs = [hs; allgrps(n,i,1)];
  end;
  [rp,cp] = size(hs);
%hs
  
  consistency = zeros(rp,1);              % Allocate consistency vector

  if (use_data)                           % Use original data
    if (length(x) ~= length(grps))
      error('  HOMOSUB: input matrices not compatible');
    end;

    t = zeros(n,n);                         % Use t-stats as intergrp distances
    for i = 1:(n-1)
      ii = find(grps == ugrps(i));
      for j = (i+1):n
        jj = find(grps == ugrps(j));
        t(j,i) = tval(x([ii;jj]),grps([ii;jj]));
      end;
    end;
    tvect = trilow(t);                      % Convert to single vector

    g = zeros(n,n);
    for ip = 1:rp                           % Cycle thru potential subsets,
      for i = 1:(n-1)                       %   scoring consistencies
        for j = (i+1):n
          g(j,i) = (hs(ip,i)~=hs(ip,j));    % 0 = same subset, 1 = diff subset
        end;
      end;
      gvect = trilow(g);
%      consistency(ip) = abs(spearman(tvect,gvect));

      tv = tvect;
      iv = find(gvect==1);
      tv(iv) = tv(iv) - mean(tv(iv));
      consistency(ip) = 1./std(tv);
    end;

  else                                    % Use presence/absence of signif
    err = 0;                                % Check for valid input matrix 
    if (n~=m)
      err = 1;
    end;
    if (~err)
      u = uniquef(S,1);
      if (u~=[0;1])
        err = 1;
      end;
    end;
    if (~err & sum(sum(S-S'))>eps)
      err = 1;
    end;
    if (err)
      error('  HOMOSUB: Input matrix must be square, symmetric, binary');
    end;

    incr = 1./(n*(n-1)/2);                  % Increment for consistency index
    delta = 0.01*incr;                      % Bias for difference over similarity
  
    for ip = 1:rp                           % Cycle thru potential subsets,
      c = 0;                                %   scoring consistencies
      for i = 1:(n-1)
        for j = (i+1):n
          if (S(i,j) & (hs(ip,i)~=hs(ip,j)))
            c = c + incr + delta;
          elseif (~S(i,j) & (hs(ip,i)==hs(ip,j)))
            c = c + incr;
          end;
        end;
      end;       
      consistency(ip) = c;
    end;
  end;
[consistency hs]

  cmax = max(consistency);                % Retain subset of max consistency
  i = find(consistency==cmax);
  hs = hs(i,:);
  consistency = cmax;

  for i = 1:size(hs,1)                    % Sequence values on each line
    if (any(hs(i,:)>1))
      j = find(hs(i,:)>0);
      h = hs(i,j);
      u = uniquef(h);
      hs(i,j) = replace(h,u,[1:length(u)]);
    end;
  end;
  
  rv = rowtoval(hs);                      % Delete any identical subsets
  [u,f] = uniquef(rv);
  if (any(f>1))
    hsave = hs;
    hs = [];
    for iu = 1:length(u)
      i = find(rv==u(iu));
      hs = [hs; hsave(i(1),:)];
    end;
  end;

  na = double('a');                       % Translate numbers into letters
  maxhs = max(max(hs));
  hs = char(replace(hs,[1:maxhs],[na:(na+maxhs-1)]));
  hs = sortrows(hs);                      % Sort into lexicological sequence
  
  return;