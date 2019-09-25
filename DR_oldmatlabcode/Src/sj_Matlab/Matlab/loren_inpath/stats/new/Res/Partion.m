% Partion:  Returns a list of all possible partitions of N linearly arranged 
%           objects into k groups.  Practical for N<=20 for a single k.
%
%     Usage: grps = partion(N,{k},{sortmat})
%           
%             N =       number of objects.
%             k =       optional number of groups [default = 2:N].
%             sortmat = optional flag indicting that the output list is to be 
%                         sorted into lexicological sequence:
%                           0 = no sort [default].
%                           1 = sort within k.
%                           2 = sort across k.
%             ---------------------------------------------------------------
%             grps =    [npart x N] matrix of partitions.
%

% RE Strauss, 8/27/99
%   9/17/99 - sort output matrix.
%   9/21/99 - allow k to be vector; make sort optional.

function grps = partion(N,k,sortmat)
  if (nargin < 2) k = []; end;
  if (nargin < 3) sortmat = []; end;

  if (isempty(k))
    k = [2:N];
  end;
  if (isempty(sortmat))
    sortmat = 0;
  end;

  sort_within = 0;
  sort_among = 0;
  switch sortmat
    case 0
    case 1, 
      sort_within = 1;
    case 2, 
      sort_among = 1;
    otherwise
      error('  PARTION: invalid sort flag');
  end;

  grps = [];

  for ik = 1:length(k)                   % Cycle thru number of partitions
    if (k(ik)==1)
      gki = ones(1,N);                    % One group
    elseif (k(ik)==N)
      gki = [1:N];                        % N groups
    else
      g = partf(N,k(ik));                 % Get list of partitions
      gki = fliplr(g);                    % Resequence
    end;

    if (sort_within)                      % Sort list for current k
      gki = sortrows(gki);
    end;

    grps = [grps; gki];
  end;

  if (sort_among)                         % Sort across all k
    grps = sortrows(grps);
  end;


  return;
