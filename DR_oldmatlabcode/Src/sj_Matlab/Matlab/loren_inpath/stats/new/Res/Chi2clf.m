% Chi2clf:  Function for chi2clst().  For a single permutation of rows, cycles 
%           thru all possible numbers of partitions from mingrps to maxgrps and 
%           returns best solutions.
%
%     Usage: [chi2df,grps] = chi2clf(table,perm,mingrps,maxgrps,chi2df,grps)
%
%           table =   [r x c] matrix of observed counts.
%           perm =    permutation of values 1:N.
%           mingrps = minimum number of groups examined.
%           maxgrps = maximum number of groups examined.
%           chi2df =  [nbest x 1] vector of largest statistic values.
%           grps =    [nbest x r] matrix of corresponding group labels for rows.
%

% RE Strauss, 8/27/99

function [chi2df,grps] = chi2clf(table,perm,mingrps,maxgrps,chi2df,grps)
  [N,c] = size(table);
  table = table(perm,:);                  % Permute table

  for k = mingrps:maxgrps                 % Cycle thru numbers of groups
    obs = zeros(k,c);                       % Allocate reduced table
    partit = partion(N,k);                  % Get partitions
    for p = 1:size(partit,1)                  % Save best solutions
      part = partit(p,:);                     % Current partition
      for g = 1:k                             % Reduced table
        i = find(part==g);
        if (length(i)==1)
          obs(g,:) = table(i,:);
        else
          obs(g,:) = sum(table(i,:));
        end;
      end;

      coltot = sum(obs);                      % Marginal totals
      rowtot = sum(obs')';
      totval = sum(sum(obs));                 % Grand total
      exp = rowtot*coltot/totval;             % Expected values
      cell_x2 = (obs-exp).^2./exp;            % Pearson chi-square statistic
      totx2 = sum(sum(cell_x2));              % Observed total chi-square
      df = (k-1)*(c-1);                       % Degrees of freedom
      stat = totx2 / df;                      % Chi-square criterion

      [minstat,i] = min(chi2df);              % Current min of best statistics
      if (stat > minstat)                     % Save results if better
        chi2df(i) = stat;
        grps(i,:) = part(perm);
      end;
    end;
  end;
  
  return;
