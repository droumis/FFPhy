% SIZEVECT: Estimates within-group PC1 loadings, by group, and converts them to 
%           allometric coefficients.  Assumes that the data matrix has been 
%           log-transformed.
%
%     Usage: load = sizevect(logX,{grps},{noconv})
%
%           logX =    [n x p] matrix of log-transformed morphometric data.
%           grps =    group-membership vector (length n) for k groups.  If not 
%                       passed, a single group is assumed.
%           noconv =  optional boolean flag indicating, if true, that the PC1 
%                       loadings are to be left as eigenvector coefficients 
%                       and not transformed to allometric coefficients 
%                       [default = 0].
%           -----------------------------------------------------------------
%           load =    [p x k] matrix of within-group coefficients; cols are in 
%                       the same sequence as the unique identifiers in 'grps'.
%

% RE Strauss, 1/15/00

function load = sizevect(logX,grps,noconv)
  if (nargin < 2) grps = []; end;
  if (nargin < 3) noconv = []; end;

  [N,P] = size(logX);

  if (isempty(grps))
    grps = ones(N,1);
  end;
  if (isempty(noconv))
    noconv = 0;
  end;

  grpid = uniquef(grps);
  ngrps = length(grpid);
  load = zeros(P,ngrps);

  for g = 1:ngrps
    load(:,g) = pcacov(logX(grps==grpid(g),:),1,1);
  end;

  if (~noconv)
    load = allom(load);
  end;

  return;
