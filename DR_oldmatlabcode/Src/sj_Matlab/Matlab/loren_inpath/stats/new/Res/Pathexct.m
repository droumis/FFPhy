% PATHEXCT: Determines the shortest (and optionally longest) paths among a small set 
%           of points by brute force, by determining the path length for all possible 
%           permutations.  A starting point must be specified, and an optional end point 
%           may be specified; if they are equal, the resulting path is a Hamiltonian circuit.  
%           Interpoint distances are given by a square distance matrix (having zeros on 
%           the diagonal), which need be neither symmetric nor complete.  Incomplete links 
%           can be specified by arbitrarily large values.
%
%           The maximum number of points feasibly examined by brute force is 8-10.
%
%     Usage: [minpath,minlen,maxpath,maxlen] = pathexct(dist,sp,{ep})
%
%           dist =    [n x n] distance matrix.
%           sp =      index of the starting point.
%           ep =      optional index of the end point.
%           -------------------------------------------------------------------------------
%           minpath = [n x 1] vector of indices of points in the sequence connected; if the 
%                       starting point is also the end point, its value is not repeated.
%           minlen =  path length of the shortest path (including the final link for a 
%                       Hamiltonian circuit).
%           maxpath = comparable path of maximum length.
%           maxlen =  path length of the longest path.
%

% RE Strauss, 7/6/98

function [minpath,minlen,maxpath,maxlen] = pathexct(dist,sp,ep)
  if (nargin < 3) ep = []; end;

  [n,m] = size(dist);

  if (isempty(ep))
    ep = 0;
  end;

  m = n-1;                                % Number of values to be permuted
  if (ep & ep~=sp)
    m = m-1;
  end;
  nperm = prod(2:m);                      % Number of possible permutations

  minlen = 1/eps;
  maxlen = eps;

  perm = 1:m;
  minpath = perm;
  maxpath = perm;

  for it = 1:nperm
    p = perm;                               % Current permuation

    i = (p >= sp);                            % Adjust indices to allow for start pt
    p(i) = p(i)+1;

    if (ep>0 & ep~=sp)                        % Adjust to allow for end pt
      i = (p >= ep);
      p(i) = p(i)+1;
      p = [sp p ep];                          % Put start & end pts at beginning and end
    else                                      %   of list
      p = [sp p];
    end;

    len = 0;  
    for i = 1:length(p)-1                     % Get path length
      len = len + dist(p(i),p(i+1));
    end;

    if (len < minlen)
      minlen = len;
      minpath = p;
    end;

    if (len > maxlen)
      maxlen = len;
      maxpath = p;
    end;

    perm = permnext(perm);                  % Next permutation
  end;

  return;

