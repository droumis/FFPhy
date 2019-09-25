% BASEDIFFP: Converts [n x m] matrix of base sequences (length m) for n taxa 
%             into an [n(n-1)/2 x 10] matrix of pairwise contrasts, for basediff().

% RE Strauss, 1/23/01

function t = basediffp(T)
  [n,m] = size(T);

  u = uniquef(T(:),1);
  lenu = length(u);
  if (lenu>4)
    error('  BASEDIFF: invalid input matrix.');
  end;

  TT = T;
  t = zeros(n*(n-1)/2,10);
  r = 0;

  for in = 1:(n-1)                    % Cycle thru pairs of taxa (rows)
    for jn = (in+1):n
      T = TT([in,jn],:);

      i = find(~isfinite(sum(T)));
      if (~isempty(i))
        T(:,i) = [];
      end;

      r = r+1;
      for i = 1:lenu
        ui = u(i);
        for j = i:lenu
          uj = u(j);
          k = find((T(1,:)==ui & T(2,:)==uj) | (T(2,:)==ui & T(1,:)==uj));
          if (~isempty(k))
            lenk = length(k);
            if (i==j)
              t(r,i+6) = lenk;
            else
              k = i+j-1;
              if (i==1)
                k = k-1;
              end;
              t(r,k) = lenk;
            end;
          end;
        end;
      end;
    end;
  end;