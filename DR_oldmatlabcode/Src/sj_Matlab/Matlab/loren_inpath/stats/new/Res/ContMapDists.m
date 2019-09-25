% CONTMAPDISTS: Least-squares mapping of distributions of a single continuous 
%               character onto a dichotomous tree via randomized resampling of 
%               terminal distributions.
%
%     Usage:[Xnode,stats,mappedstats] = ...
%                        contmapdists(anc,X,iter,{doplots},{rooted},{brlen},{k})
%
%           anc =     ancestor function describing tree topology for n taxa; the 
%                       root must have ancestor 0.
%           X =       [n x m] data matrix for n taxa (rows), with a maximum of 
%                       m observations per taxon; use NaN values to fill out the 
%                       matrix for varying sample sizes per taxon.
%           iter =    number of resampling iterations.
%           doplots = optional boolean flag indicating that plots are to be 
%                       produced for all terminal and nodal distributions, to a 
%                       common scale [default = 0].
%           rooted =  optional boolean flag indicating whether tree root is (=1) 
%                       or is not (=0) to be optimized; if not, the tree is
%                       treated as an arbitrarily rooted network [default = 1].
%           brlen =   optional branch lengths for weighted least-squares 
%                       mapping [default = null].
%           k =       optional value of Minkowski k [default = 2 for  
%                       least-squares solution].
%           ---------------------------------------------------------------------
%           Xnode =   [n-1 x iter] data matrix containing predicted observations
%                       for n-1 nodes, with 'iter' observations per distribution.
%           stats =   [2n-1 x 7] matrix of distributional statistics (from 
%                       univar) for the terminal + nodal distributions.
%           mappedstats = [2n-1 x 9] matrix of distributional statistics (from 
%                       univar) mapped onto the tree, rather than discovered 
%                       from randomized distributions.
%

% RE Strauss, 5/22/01

function [Xnode,stats,mappedstats] = contmapdists(anc,X,iter,doplots,rooted,brlen,k)
  if (nargin < 4) doplots = []; end;
  if (nargin < 5) rooted = []; end;
  if (nargin < 6) brlen = []; end;
  if (nargin < 7) k = []; end;

  get_stats = 0;
  get_mappedstats = 0;
  nstats = 7;
  if (nargout > 1)
    get_stats = 1;
    get_mappedstats = 1;
  end;

  if (isempty(doplots))
    doplots = 0;
  end;

  [n,m] = size(X);
  nnodes = n-1;
  sampsize = sum(isfinite(X'))';

  Xn = zeros(n-1,iter);
  x = zeros(n,1);

  for it = 1:iter
    i = ceil(rand(n,1).*sampsize);
    for j = 1:n
      x(j) = X(j,i(j));
    end;
    Xnode(:,it) = contmap(anc,x,rooted,brlen,k);
  end;

  if (get_stats)
    bs = ones(1,nstats);
    s = univar(X',bs)';
    stats = [s; univar(Xnode',bs)'];
  end;
  if (get_mappedstats)
    ms = contmap(anc,s,rooted,brlen,k);
    mappedstats = [s; ms];
  end;

  if (doplots)
    x = X';
    xn = Xnode';
    histdata = [x(:); xn(:)];
    g = makegrps(1:(2*n-1),[sampsize; iter*ones(nnodes,1)]);

    [nplot,freqs] = histgram(histdata,g,1:(2*n-1),[],[1 1],[],'rel');
  end;

  return;


