% BINMATRN: Creates a randomized binary matrix (containing 0's and 1's), subject to 
%           one of three levels of marginal constraint.
%
%     Usage: matrix = binmatrn(matsize,total,{constr})
%                         OR
%            matrix = binmatrn(rowtot,coltot,constr)
%
%             if constr = 0:
%               matsize = 2-element vector [r,c] indicating the numbers of rows and columns
%                           for the random matrix.
%               total =   total number of 1's to be randomly placed in matrix.
%
%             if constr = 1 or 2:
%               rowtot =  vector of row totals.
%               coltot =  vector of column totals.
%
%             constr =          flag indicating the level of marginal constraints:
%                                 0 = unconstrained  [default].
%                                 1 = marginal totals proportional to probabilities.
%                                 2 = marginal totals fixed.
%             --------------------------------------------------------------------------
%             matrix =          randomized binary matrix.
%

function matrix = binmatrn(rowtot,coltot,constr)
  if (nargin < 3)
    constr = [];
  end;
  if (isempty(constr))
    constr = 0;
  end;

  if (constr == 0)
    if (length(rowtot) ~= 2)
      error('  BINMATRN error: "matsize" must be 2-element vector');
    end;
    nrows = rowtot(1);
    ncols = rowtot(2);
    total = coltot;
  else
    total = sum(rowtot);
    if (sum(coltot) ~= total)
      error('  BINMATRN error: row totals and column totals must agree');
    end;
    nrows = length(rowtot);
    ncols = length(coltot);
  end;

  matrix = zeros(nrows,ncols);      % Allocate output matrix

  if (constr == 0)                  % Random matrix unconstrained by marginal totals
    v = randperm(nrows*ncols);        % Randomly permuted vector of cell addresses
    celladdr = v(1:total)';           % Address to be filled with 1's

    c = rem(celladdr,ncols);          % Convert cell addrs to rows & cols
    i = find(c==0);
    c(i) = ncols*ones(length(i),1);
    r = (celladdr - c)/ncols + 1;
    
    for i = 1:total                   % Fill matrix with 1's
      matrix(r(i),c(i)) = 1;
    end;

  elseif (constr == 1)              % Marginal totals proportional to probabilities



  elseif (constr == 2)              % Fixed marginal totals
    obs = ones(total,1);         % Col vector of ones
    cumcoltot = cumsum(coltot);       % Cumulative sums of col totals
%    cumrowtot = cumsum(rowtot);

    b = 1;
    for c = 1:ncols                   % Label obs with col ids
      e = cumcoltot(c);
      obs([b:e]) = obs([b:e])*c;
      b = e+1;
    end;
%    p = randperm(total);              % Random permutation of col-id vector
%    obs = obs(p);
rowtot
coltot

    for r = 1:nrows
      rt = rowtot(r);
      o = obs;
      for s = 1:rt
        i = find(o>0);
if (isempty(i))
  matrix
  r
  s
  obs = obs'
  o = o'
end;
        i = i(ceil(rand*length(i)));
        matrix(r,o(i)) = 1;
        obs(i) = 0;
        if (s < rt)
          j = find(o == o(i));
          o(j) = zeros(length(j),1);
        end;
      end;
    end;

  else
    error('  BINMATRN error: invalid constraint value');
  end;

  return;
