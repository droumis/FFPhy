% HOMOSUB:  Composes a matrix of homogeneous subsets (rows), coded in lower-case 
%           letters, given either (1) an NxN symmetric binary matrix indicating 
%           which of N groups are significantly different from one another, or 
%           (2) the original data, from which means and standard errors are 
%           estimated.  The first method is qualitative, based on presence/absence 
%           of significant difference, and the measure of consistency is the 
%           proportion of presence/absence significance states that is in 
%           agreement with the assignment of groups to homogeneous subsets.
%           The second method is quantitative, based on pairwise t-tests with a 
%           sequential Bonferroni correction.  If subgroup identifiers are 
%           provided (for a nested design), nested anovas are done rather than 
%           t-tests, again with a sequential Bonferroni correction.
%             Subsets often comprise ambiguous, overlapping (intransitive) sets 
%           of groups.  Given are the subsets most consistent with the input matrix;  
%           all subsets returned are equally consistent.
%
%     Usage: [hs,consistency] = homosub(S)
%                   OR
%            [hs,consistency] = homosub(x,grps,{subgrps},{alpha})
%
%           S =       [n x n] symmetric boolean matrix indicating which pairwise 
%                       comparisons are significantly different.
%                         OR
%           x =       vector (length n) of data.
%           grps =    corresponding group-membership vector.
%           subgrps = optional corresponding vector of subgroup membership.
%           alpha =   optional overall criterion for statistical significance 
%                       [default = 0.05].
%           ---------------------------------------------------------------------
%           hs =      [ns x n] character matrix with possible homogeneous subsets.
%           consistency = measure of total pairwise consistencies [0-1] 
%                           of subsets with S.
%

% RE Strauss, 12/25/99
%   7/24/00 -  for binary comparisons, give sig difference a slight bias over 
%                non-sig difference in consistency measure.
%   10/17/00 - for data input, change from original measure of consistency (rank 
%                correlation between pairwise t-statistic values and binary 
%                vector indicating membership in same or different groups) to 
%                one based solely on pairwise t-tests with sequential 
%                Bonferroni correction.
%   10/19/00 - added nested anova in presence of subgroups within groups.

function [hs,consistency] = homosub(x,grps,subgrps,alpha)
  if (nargin < 2) grps = []; end;
  if (nargin < 3) subgrps = []; end;
  if (nargin < 4) alpha = []; end;

  if (isempty(alpha))
    alpha = 0.05;
  end;
  if (alpha > 1)
    alpha = alpha/100;
  end;

  if (isempty(grps))                      % If input is significance matrix,
    S = x;                                  % Switch assignment of variables
    if (~issqsym(S))                        % Check for valid input matrix 
      error('  HOMOSUB: Input matrix must be square, symmetric, binary');
    end;
  else                                    % Else if input is data,
    x = x(:);                               % Convert input to col vectors
    grps = grps(:);

    if (~isempty(subgrps))                  % Do pairwise tests
      ugrps = uniquef(grps,1);                % Group identifiers
      ngrps = length(ugrps);
      subgrps = subgrps(:);                   
      pr = zeros(ngrps*(ngrps-1)/2,1);        % Allocate probability matrix
      p = 0;
      for i = 1:(ngrps-1)                     % Pairwise nested anovas if have subgroups
        for j = (i+1):ngrps
          k = find(grps==ugrps(i) | grps==ugrps(j));
          p = p+1;
          if (length(uniquef(subgrps(k)))>1)
            [F,prp] = anovanst(x(k),grps(k),subgrps(k));
            pr(p) = prp(1);
          else
            [F,pr(p)] = anova(x(k),grps(k));
          end;
        end;
      end;
      signif = (pr <= alpha);
    else
      [t,pr,df,signif] = ttest(x,grps,[],alpha);  % Or t-tests if not
    end;
    S = trisqmat(signif);                       % Form significance matrix
  end;
  [n,m] = size(S);

  hs = [];                                % Find all potential subsets
%   hs = [];
%   for i = 1:n
%     hs = [hs; allgrps(n,i,1)];
%   end;
  hs = allgrps(n,1:n,1);
  [rp,cp] = size(hs);
%hs
  
  consistency = zeros(rp,1);              % Allocate consistency vector
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
%[consistency hs]

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