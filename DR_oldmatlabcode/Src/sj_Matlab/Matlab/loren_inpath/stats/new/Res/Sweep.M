% SWEEP: Sweep operator for symmetric matrices.  Sweeps the matrix on row
%         and column k.
%
%     Usage:  H = sweep(G,k)
%
%         G = [p x p] symmetric matrix.
%         k = vector of indices of row(s) & col(s) to be swept from G.  
%               The sweep operator is commutative for any permutation of k.
%         -----------------------------------------------------------------
%         H = residual matrix.
%

% Little, RJA and DB Rubin. 1987. Statistical Analysis with Missing Data.
%   Wiley.  P. 112-119.

% RE Strauss, 7/5/95

function H = sweep(G,kk)
  [n,p] = size(G);
  if (n~=p)
    error('  SWEEP: input matrix must be symmetric.');
  end;

  H = G;
  for k = kk                        % Rows/cols to be swept
    H(k,k) = -1/G(k,k);               % Diagonal element k
    for j = 1:p                       % Elements in row & col k
      if (j ~= k)
        H(j,k) = -G(j,k)*H(k,k);
        H(k,j) = H(j,k);
      end;
    end;
    for j = 1:p                       % Elements not in row nor col k
      if (j ~= k)
        for l = 1:p
          if (l ~= k)
            H(j,l) = G(j,l) - H(j,k)*G(k,l);
            H(l,j) = H(j,l);
          end;
        end;
      end;
    end;
    G = H;                            % Save current sweep
  end;

  return;
