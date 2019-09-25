% HomoPropAIC: Estimtes homogeneous subsets of proportions based on Akaike 
%              information criteria.
%
%     Usage: [hs,crit,ordgrps] = homopropaic(p,n,{grps},{kind},{nbest})
%
%         p =     vector (length k) of proportions.
%         n =     corresponding vector of sample sizes.
%         grps =  optional corresponding vector of group labels 
%                     [default = 1,2,...,k].
%         kind =  optional flag indicating the information statistic to be used:
%                     0 = AIC (Akaike information criterion) [default];
%                     1 = BIC (Schwarz criterion, with penalty term of log(N));
%                     2 = CAIC (Bozdogan criterion, with penalty term of log(N)+1).
%         nbest =    optional number of best subsets to return [default = 1].
%         -------------------------------------------------------------------------
%         hs =      [nbest x k] matrix of best homogeneous subsets.
%         crit =    [nbest x 1] vector of criterion values.
%         ordgrps = group identifiers ordered by increasing proportion.
%         

% Dayton CM. 1998. Information criteria for the paired-comparisons problem.
%   American Statistician 52:144-151.
% Dayton CM. 2001. SUBSET: best subsets using information criteria.
%   Journal of Statistical Software 6(2):1-10.
% Schwarz, G. 1978. Estimating the dimension of a model.
%   Annals of Statistics 6:461-464.
% Bozdogan, H. 1987. Model-selection and Akaike's information criterion (AIC):
%   the general theory and its analytical extensions.  Psychometrika 51:345-370.

% RE Strauss, 2/16/02

function [hs,crit,ordgrps] = homopropaic(p,n,grps,kind,nbest)
  if (nargin < 3) grps = []; end;
  if (nargin < 4) kind = []; end;
  if (nargin < 5) nset = []; end;

  p = p(:);
  n = n(:);
  k = length(p);

  if (isempty(grps))
    grps = [1:k]';
  end;
  if (isempty(kind))
    kind = 0;
  end;
  if (isempty(nbest))
    nbest = 1;
  end;

  if (length(n)~=k | length(grps)~=k)
    error('  HomoPropAIC: input vectors not compatible');
  end;

  hs = inf*ones(nbest,1);               % Initialize output matrices
  crit = zeros(nbest,k);

  [p,n,ordgrps] = sortmat(p,n,grps);    % Sort by increasing proportion

  bsets = dec2bin(0:((2.^k)-2));        % Binary representations of possible subsets
  bsets = bsets(1:2.^(k-1),:);
  nsets = size(bsets,1);
  gsets = ones(nsets,k);
bsets

  for i = 1:nsets                       % Convert to ordered list of partitions
    for j = 2:k
      if (bsets(i,j)~=bsets(i,j-1))
        gsets(i,j:end) = gsets(i,j:end)+1;
      end;
    end;
  end;
    
  for i = 1:nsets                       % Cycle thru all possible partitions
    loglik = 
  end;

  return;


