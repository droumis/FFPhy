% Allgrpf:  Recursive function for listing all possible groupings of N objects 
%           into k groups.  See allgrps().
%
%     Usage: grps = allgrpf(N,k,unlabeled)
%
%             N =         vector of remaining object identifiers.
%             k =         vector of sample sizes per group.
%             unlabeled = boolean flag indicating that groups are unlabeled 
%                           (ie, groups of equal sample size are equivalent and 
%                           indistinguishable).
%             -----------------------------------------------------------------
%             grps =      [ngrp x N] matrix of partitions.
%

% RE Strauss, 8/29/99

function grps = allgrpf(N,k,unlabeled)
  lenk = length(k);

  if (lenk == 1)
    g = combvals(length(N),k);
    grps = zeros(size(g));
    for i = 1:size(g,1)
      grps(i,:) = N(g(i,:));
    end;
    return;
  end;

  grps = [];
  gp = combvals(length(N),k(1));
  for igp = 1:size(gp,1)
    n = N;
    Nout = n(gp(igp,:));
    n(gp(igp,:)) = [];

    g = allgrpf(n,k(2:lenk),unlabeled);
    grps = [grps; ones(size(g,1),1)*Nout g];
  end;

  if (unlabeled)
    ngrps = size(grps,1);
    if (ngrps > 1)
      k1 = k(1);
      k2 = k(2);
      if (k1 == k2)
        if (k1 == 1)
          i = find(grps(:,1) < grps(:,2));
          grps = grps(i,:);
        else
          i = max(find(grps(:,1)==grps(1,1)));
          i1 = [1:i]';
          i2 = [(i+1):ngrps]';
          gk1 = grps(i1,1:(k1+k2));
          gk2 = grps(i2,[(k1+1):(k1+k2) 1:k1]);
          keep = ~ismember(gk2,gk1,'rows');
          i = [i1; i+find(keep)];
          grps = grps(i,:);
        end;
      end;
    end;
  end;

  return;
