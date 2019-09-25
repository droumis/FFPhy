% Allgrps:  Returns a list of all possible groupings of N objects into k groups, 
%           labeled or unlabeled.  If group sample sizes are not specified, then 
%           all possible arrangements are considered.
%
%     Usage: grps = allgrps(N,k,{labeled},{sortmat})
%           
%             N =         number of objects.
%             k =         number of groups, or vector of group sample sizes.
%             labeled = boolean flag indicating that groups are labeled 
%                           (ie, groups of equal sample size are distinguishable and 
%                           equivalent) [default=0=unlabeled].
%             sortmat =   boolean flag indicating that output matrix is to be 
%                           sorted into lexicological sequence [default=0].
%             -----------------------------------------------------------------
%             grpid =     [m x N] matrix of group identifiers.
%

% RE Strauss, 8/29/99
%   10/3/99 -  added option of labeled groups.
%   10/27/99 - added option to sort output matrix.

function grpid = allgrps(N,k,labeled,sortmat)
  if (nargin < 3) labeled = []; end;
  if (nargin < 4) sortmat = []; end;

  if (isempty(labeled))
    labeled = 0;
  end;
  if (isempty(sortmat))
    sortmat = 0;
  end;

  k = k(:);
  kn = length(k);

  if (kn==1 & k==1)
    grpid = ones(1,N);
    return;
  end;
  if (kn==1 & k==N)
    grpid = [1:N];
    return;
  end;
  if (kn>1 & all(k==1))
    grpid = permlist(kn);
    return;
  end;

  if (kn==1)
    sizes = allsizes(N,k);                % Get list of all possible sample sizes
  else
    sizes = k;
  end;
  nsize = size(sizes,1);                  % Number of sample-size combinations

  grpid = [];
  N = 1:N;
  for i = 1:nsize                         % Cycle thru sample-size combinations
    k = sizes(i,:);                       % For current sample-size combination,
    g = allgrpf(N,k,~labeled);            %   get list of object permutations
      
    gid = zeros(size(g));                 % Convert objects to group identifiers
    leng = size(g,1);
    lenk = length(k);

    for ig = 1:leng
      kmax = 0;
      for ik = 1:lenk
        kmin = kmax+1;
        kmax = kmin+k(ik)-1;
        gid(ig,g(ig,kmin:kmax)) = ik*ones(1,k(ik));
      end;
    end;

    grpid = [grpid; gid];                 % Concatenate to full list
  end;

  if (sortmat)                          % Sort into lexicological sequence
    grpid = sortrows(grpid);
  end;

  return;
