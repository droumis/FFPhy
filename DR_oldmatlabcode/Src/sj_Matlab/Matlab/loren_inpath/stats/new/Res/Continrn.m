% CONTINRN: Generates a pseudorandom [r x c] contingency table
%           with either fixed or floating marginal (row and column) totals
%
%     Syntax: matrix = continrn(rowtot,coltot,{fixed})
%
%            rowtot = vector of row totals.
%            coltot = vector of column totals.
%            fixed =  optional flag indicating whether marginal totals are to
%                       be treated as fixed (=TRUE) or probabilistic (=FALSE)
%                       [default = TRUE].
%            -----------------------------------------------------------------
%            matrix = uniformly random contingency table.
%

% Romesburg, HC & K Marshall. 1985. CHITEST: A Monte-Carlo computer program
%   for contingency table tests.  Computers & Geosciences 11:69-78.

% For fixed marginal totals, uses method of:
%   Boyett,JM. 1979.  Algorithm AS 144: Random RxC tables with given row and
%     column totals.  Appl. Statist. 28:329-332.
%   Agresti,A, D Wackerly, JM Boyett.  1979.  Exact conditional test for cross-
%     classifications: approximation of attained significance levels.
%     Psychometrika 44:75-83.

% See also:
%   Patefield,WM. 1981. Algorithm AS 159: An efficient method of generating
%     random r x c tables with given row and column totals.
%     Appl. Statist. 30:91-97.
%   Patefield,WM. 1982. Exact tests for trends in ordered contingency tables.
%     Appl Statist. 31:32-43.

% RE Strauss, 4/22/98

function matrix = continrn(rowtot,coltot,fixed)
  TRUE = 1; FALSE = 0;

  if (nargin<3)
    fixed = TRUE;
  end;

  nrow = length(rowtot);
  ncol = length(coltot);
  ncell = nrow * ncol;
  total = sum(rowtot);
  if (total ~= sum(coltot))
    error('  CONTINRN error: sums of row totals and column totals not equal');
  end;

  if (fixed)                          % Fixed marginal totals
    matrix = zeros(nrow,ncol);          % Allocate random-table matrix
    obs = ones(total,1);                % Col vector of ones
    cumcoltot = cumsum(coltot);         % Cumulative sums of col totals
    cumrowtot = cumsum(rowtot);

    b = 1;
    for j = 1:ncol                      % Label obs with col ids
      e = cumcoltot(j);
      obs([b:e]) = obs([b:e])*j;
      b = e+1;
    end;
    p = randperm(total);                % Random permutation of col-id vector
    obs = obs(p);

    b = 1;
    for i = 1:nrow
      e = cumrowtot(i);
      o = obs([b:e]);
      b = e+1;

      for j = 1:ncol
        matrix(i,j) = sum(o==j);
      end;
    end;
  else                                % Floating marginal totals
    matrix = zeros(ncell,1);            % Allocate random-table matrix (as vector)
    prob =   zeros(nrow,ncol);          % Allocate probability matrix

    for i = 1:nrow
      for j = 1:ncol
        prob(i,j) = rowtot(i)*coltot(j)/total^2;  % Expected proportions
      end;
    end;
    prob = cumsum(prob(:));             % Reshape as vector and accumulate

    for n = 1:total                     % Randomly allocate obs to cells
      r = rand;                         %   in proportion to expected values
      i = (sum(prob<=r))+1;
      matrix(i) = matrix(i)+1;
    end;

    matrix = reshape(matrix,nrow,ncol); % Reshape from vector to matrix
  end;

  return;

