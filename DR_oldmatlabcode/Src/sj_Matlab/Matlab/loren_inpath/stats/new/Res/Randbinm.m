% RANDBINM: Generates a pseudorandom [r x c] binary matrix (containing 0's and 1's)
%           with either fixed, probabilistic, or unconstrained (random binomial)
%           marginal (row and column) totals.
%
%     Syntax: matrix = randbinm(rows,cols,N)
%
%            rows = vector of row totals or probabilities, or scalar
%                     indicating number of rows for unconstrained row totals.
%            cols = vector of column totals or probabilities, or scalar
%                     indicating number of columns for unconstrained column totals.
%            N =    total number of 1's in matrix (optional for fixed marginal
%                     totals, required for probabilistic or unconstrained totals).
%            ----------------------------------------------------------------------
%            matrix = random binary table.
%

% RE Strauss, 5/4/98
%   1/16/99 - allow for row or col totals, but not both, to be constrained.

function matrix = randbinm(rows,cols,N)
  if (nargin < 3) N = []; end;

  lr = length(rows);
  lc = length(cols);

  constr_rows = 0;                      % Flags
  constr_cols = 0;
  needN = 0;

  if (lr==1 & lc==1)                    % Unconstrained
    nrows = rows;
    ncols = cols;
    if (isempty(N))
      needN = 1;
    end;

  elseif (lr>1 & lc==1)                 % Constrained rows
    nrows = lr;
    ncols = cols;
    constr_rows = 1;

    if (all(rows<1))                      % Convert probs to counts
      if (isempty(N))
        needN = 1;
      else
        rows = prbcount(rows,N,ncols);
      end;
    else
      if (any(rows>ncols))
        error('  RANDBINM: At least one row total exceeds number of columns');
      end;
      if (isempty(N))
        N = sum(rows);
      else
        if (sum(rows)~=N)
          error('  RANDBINM: Sum of row totals not equal to N');
        end;
      end;
    end;

  elseif (lr==1 & lc>1)                 % Constrained cols
    nrows = rows;
    ncols = lc;
    constr_cols = 1;

    if (all(cols<1))                      % Convert probs to counts
      if (isempty(N))
        needN = 1;
      else
        cols = prbcount(cols,N,nrows);
      end;
    else
      if (any(cols>nrows))
        error('  RANDBINM: At least one column total exceeds number of rows');
      end;
      if (isempty(N))
        N = sum(cols);
      else
        if (sum(cols)~=N)
          error('  RANDBINM: Sum of column totals not equal to N');
        end;
      end;
    end;

  else % if (lr>1 & lc>1)               % Constrained rows & cols
    nrows = lr;
    ncols = lc;
    constr_rows = 1;
    constr_cols = 1;

    if (any(rows>ncols))
      error('  RANDBINM: At least one row total exceeds number of columns');
    end;
    if (any(cols>nrows))
      error('RANDBINM: At least one column total exceeds number of rows');
    end;

    if (all(rows>=1) & all(cols>=1))      % Counts
      if (sum(rows) ~= sum(cols))
        error('  RANDBINM: Sums of row and column marginals not equal');
      end;
      if (~isempty(N))
        if (sum(rows) ~= N)
          error('  RANDBINM: Sums of row and column marginals not equal to N');
        end;
      else
        N = sum(rows);
      end;
    else                                  % Probabilities
      rows = prbcount(rows,N,ncols);
      cols = prbcount(cols,N,nrows);
    end;
  end;

  if (needN)
    error('  RANDBINM: N = number of 1s is required');
  end;

  ncells = nrows*ncols;
  if (N > ncells)
    error('  RANDBINM: N exceeds matrix size');
  end;


  % Unconstrained marginal totals

  if (~constr_rows & ~constr_cols)      % Unconstrained marginals
    c = [ones(N,1); zeros(ncells-N,1)];   % Vector of ones and zeros
    c = c(randperm(ncells));              % Randomly permute vector
    matrix = reshape(c,nrows,ncols);      % Reshape into matrix
    return;
  end;


  % Constrained row totals only

  if (constr_rows & ~constr_cols)
    matrix = zeros(nrows,ncols);        % Allocate matrix
    for i = 1:nrows                     % Construct randomized matrix by row
      nr = rows(i);
      c = [ones(1,nr) zeros(1,ncols-nr)]; % Vector of ones and zeros
      matrix(i,:) = c(randperm(ncols));   % Randomly permute vector
    end;
    return;
  end;


  % Constrained column totals only

  if (~constr_rows & constr_cols)
    matrix = zeros(nrows,ncols);        % Allocate matrix
    for j = 1:ncols                     % Construct randomized matrix by col
      nc = cols(j);
      c = [ones(nc,1); zeros(nrows-nc,1)];  % Vector of ones and zeros
      matrix(:,j) = c(randperm(nrows));     % Randomly permute vector
    end;
    return;
  end;


  % Constrained row and col totals

  matrix = zeros(nrows,ncols);          % Allocate matrix
  restart = 1;
  
count = 0;
  while (restart)
    restart = 0;

    [cols,indx] = sort(-cols);          % Sort cols by descending magnitude
    cols = -cols;

    for i = 1:nrows                     % Fill rows with required 0's & 1's
      matrix(i,:) = [ones(1,rows(i)) zeros(1,ncols-rows(i))];
    end;
count = count+1;

    for j = 1:(ncols-1)                   % Shift zeros left randomly within rows
      while ((sum(matrix(:,j))>cols(j)) & ~restart)  %   to reduce column sums
        withones = find(matrix(:,j)==1);    % Rows containing ones

        subm = matrix(withones,(j+1):ncols);
        if (j+1 < ncols)
          r0 = sum((subm==0)');
        else
          r0 = (subm==0)';
        end;

if (count > 1000)
  matrix
  nrows
  ncols
  j
  withones
  subm
  r0
  error('Breaking randbinm');
end;
        if (sum(r0)>0)
          cumrow = [0, cumsum(r0)/sum(r0)];
          k = [];
          while (isempty(k))
            ii = max(find(rand>cumrow));        % Select row in proportion to number of 0's
            k = min(find(subm(ii,:)==0))+j;     % Find first zero in row after current col
          end;
          i = withones(ii);
          matrix(i,k) = 1;                    % Exchange
          matrix(i,j) = 0;
        else
          restart = 1;
        end;
      end;  % while (sum(matrix(:,j)) > cols(j) & ~restart)
    end;  % for j = 1:(ncols-1)

    if (~restart)
      t = matrix;                         % Re-sort columns
      for j = 1:ncols
        matrix(:,indx(j)) = t(:,j);
      end;
    end;

  end;  % while (restart)

  [i,j] = find(matrix==1);
  if (length(i)~=N)
    matrix
    error('  RANDBINM: Error in constructing random matrix');
  end;

  return;
