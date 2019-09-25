% TRISQMAT: Given a column vector representing the lower-triangular (by
%           column) or upper-triangular (by row) portion of a square matrix,
%           or a 2-column matrix representing both, creates the square matrix.
%           Optionally accepts a scalar or vector for the diagonal.
%
%       Usage: sqmat = trisqmat(tri,{dg})
%
%           tri =   column vector(s) for lower-triangular and (optionally)
%                     upper-triangular values.
%           dg =    optional scalar or vector for the diagonal [default = 0].
%           -----------------------------------------------------------------
%           sqmat = corresponding square matrix.
%

% RE Strauss, 10/14/97

function sqmat = trisqmat(tri,dg)
  [r,c] = size(tri);

  if (r==1 & c>2)                     % Transpose row vector
    tri = tri';
    [r,c] = size(tri);
  elseif (r==2 & c>2)
    tri = tri';                       % Transpose 2-row matrix
    [r,c] = size(tri);
  elseif (c>2)
    error('  Invalid input-matrix dimensions');
  end;

  found_p = 0;                        % Find size of square matrix
  P = 1;
  while (P<1000 & ~found_p)
    P = P+1;
    if (sum(1:(P-1)) == r)
      found_p = 1;
    end;
  end;

  if (~found_p | c>3)                 % Check input parameters
    error('  Error: Triangular-matrix columns too many or wrong length');
  end;

  if (nargin > 1)                     % Vector for diagonal
    d = max(size(dg));
    if (d==1)                           % Expand scalar
      dg = dg * ones(P,1);
    else
      if (min(size(dg)) > 1)
        error('  Error: value for diagonal must be vector or scalar');
      end;
    end;
  else
    dg = zeros(P,1);                    % Default diagonal
  end;

  sqmat = eye(P)*diag(dg);            % Initialize square-matrix diagonal

  k = 0;                              % Stash lower-triangular values
  for j = 1:(P-1)                     %   in upper triangle also
    for i = (j+1):P
      k = k+1;
      t = tri(k,1);
      sqmat(i,j) = t;
      sqmat(j,i) = t;
    end;
  end;

  if (c>1)                            % If a second tri vector given,
    k = 0;                            %   stash it in upper triangle
    for i = 1:(P-1)
      for j = (i+1):P
        k = k+1;
        sqmat(i,j) = tri(k,2);
      end;
    end;
  end;

  return;
