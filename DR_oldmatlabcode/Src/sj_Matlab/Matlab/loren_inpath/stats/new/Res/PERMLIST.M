% PERMLIST: Enumerates the complete set of all N! permutations of N objects,
%           returning them in an N!xN matrix P.  Practical up to about N=7.
%
%       Usage: P = permlist(N)
%

% Based on Algorithm AS179, Berry (1982), Appl.Stat. 31:169-173.

% RE Strauss, 2/11/95

function P = permlist(N)
  TRUE = 1; FALSE = 0;
  first = TRUE;
  br = FALSE;

  K = N;
  A = zeros(K,N);
  P = zeros(prod(1:N),N);
  p = 0;

  for r=1:K-1                       % Step 1
    A(r,1) = r;
    A(r,2:N-r+1) = (r+1)*ones(1,N-r);
  end;

  while (1)                         
    i = K-1;                        % Step 2

    if (first)
      first = FALSE;                
      A(K,:) = A(1,:);              % Step 5
      for r=2:K-1
        cc = 0;
        for c=1:K
          if (A(K,c)==r)
            cc = cc+1;
            A(K,c) = A(r,cc);
          end;
        end;
      end;

    else; % not first
      while (~br)
        found = FALSE;              % Step 3
        nz = nnz(A(i,:));           % Number of nonzero elements in row
        for j=1:nz-1
          if (all([A(i,j)==i, A(i,j+1)==i+1]))
            A(i,j) = i+1;
            A(i,j+1) = i;
            limit = j-2;
            found = TRUE;
            break;
          end;
        end;

        if (all([found, limit>0]))  % Left-shift marks out of sequence
          flag = TRUE;
          while (flag)
            flag = FALSE;
            for j=1:limit
              if (all([A(i,j)==i+1, A(i,j+1)==i]))
                A(i,j) = A(i,j+1);
                A(i,j+1) = A(i,j)+1;
                flag = TRUE;
              end;
            end;
          end;
        end;

        if (found)
          A(K,:) = A(1,:);           % Step 5
          for r=2:K-1
            cc = 0;
            for c=1:K
              if (A(K,c)==r)
                cc = cc+1;
                A(K,c) = A(r,cc);
              end;
            end;
          end;
          br = TRUE;
        end;

        if (~br)
          if (all([~found, i==1]))  % Step 4
            return;
          else
            A(i,1) = i;
            A(i,2:N-i+1) = (i+1)*ones(1,N-i);
            i = i-1;
          end;
        end;
      end; %while (~br)
    end;

    br = FALSE;
    p = p+1;                        % Step 6
    P(p,:) = A(K,:);
  end; %while (1)

