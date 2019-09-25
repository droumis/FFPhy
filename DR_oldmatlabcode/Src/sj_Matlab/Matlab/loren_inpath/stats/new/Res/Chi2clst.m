% Chi2clst: Cluster analysis of rows of a 2-way contingency table of counts.
%           Proceeds by pooling rows so as to maximize chi2/df.  Exact, brute-
%           force solutions are practical for only 9 or so rows; for >9 rows,
%           random permutations are used and the optimal solution is not 
%           guaranteed.  Exact solutions for the single unique input sequence of 
%           rows are reasonable for N<=20 or so.
%
%     Usage: [chi2df,grps] = chi2clst(table,{mingrps},{maxgrps},{keepseq},{nbest})
%
%           table =   [r x c] matrix of observed counts.
%           mingrps = minimum number of groups examined [default = 2].
%           maxgrps = maximum number of groups examined [default = r].
%           keepseq = boolean flag indicating that sequence of rows is to be 
%                       maintained, for linearly continguous clustering.
%           nbest =   number of best solutions to be reported [default = 10].
%           ---------------------------------------------------------------------
%           chi2df =  [nbest x 1] vector of largest statistic values, sorted 
%                       high to low.
%           grps =    [nbest x r] matrix of corresponding group labels for rows.
%

% RE Strauss, 8/27/99

function [chi2df,grps] = chi2clst(table,mingrps,maxgrps,keepseq,nbest)
  if (nargin < 2) mingrps = []; end;
  if (nargin < 3) maxgrps = []; end;
  if (nargin < 4) keepseq = []; end;
  if (nargin < 5) nbest = []; end;

  [N,c] = size(table);
  if (min([N,c]) < 2)
    error('CHI2CLST: 2-way contingency table required');
  end;

  if (isempty(mingrps))                   % Default arguments
    mingrps = 2;
  end;
  if (isempty(maxgrps))
    maxgrps = N;
  end;
  if (isempty(keepseq))
    keepseq = 0;
  end;
  if (isempty(nbest))
    nbest = 10;
  end;

  if (mingrps < 2)
    mingrps = 2;
  end;

  chi2df = zeros(nbest,1);                % Allocate output matrices
  grps = zeros(nbest,N);

  perm = [1:N];
  [chi2df,grps] = chi2clf(table,perm,mingrps,maxgrps,chi2df,grps);  % Initial solution

  if (~keepseq)
    if (N > 9)                            % If >9 rows in table,
      iter = prod(1:9);                     % Randomly iterate equivalent of 9 rows
      for it = 1:iter
        perm = randperm(N);                 % Randomly permuted table
        [chi2df,grps] = chi2clf(table,perm,mingrps,maxgrps,chi2df,grps);
      end;
    else
      iter = prod(1:N)-1;
      for it = 1:iter
        perm = permnext(perm);              % Systematically permuted table
        [chi2df,grps] = chi2clf(table,perm,mingrps,maxgrps,chi2df,grps);
      end;
    end;
  end;

  i = find(chi2df==0);                    % Remove excess statistics
  if (~isempty(i))
    chi2df(i) = [];
    grps(i,:) = [];
  end;

  [chi2df,i] = sort(-chi2df);             % Sort statistics high to low
  chi2df = -chi2df;
  grps = grps(i,:);

  return;
