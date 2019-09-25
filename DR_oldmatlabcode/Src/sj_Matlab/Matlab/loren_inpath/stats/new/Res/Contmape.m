% CONTMAPE: Least-squares mapping of continuous characters onto a completely
%           dichotomous tree, using the analytic solution of 
%           McArdle & Rodrigo (1994) for a completely bifurcated tree.
%           If branch lengths are provided for weights, the weights are 
%           specified as 1/(b^2).  The mapping thus corresponds to a Brownian-walk 
%           model if the squared length of each branch is equal to the expected 
%           increase in variance of the character along that branch.
%           Called by function contmap(), but can be used stand-alone.
%
%     Usage: Xnode = contmape(anc,X,{rooted},{brlen})
%
%           anc =     ancestor function describing tree topology; the 
%                      root must have ancestor 0.
%           X =       [n x p] data matrix for n taxa and p characters.
%           rooted =  boolean flag indicating whether tree root is (=1) or
%                       is not (=0) to be optimized; if not, the tree is
%                       treated as an arbitrarily rooted network [default=1].
%           brlen =   optional branch lengths for weighted least-squares 
%                       mapping, in same sequence as anc [default = 1 for all 
%                       branches except ancestor].
%           ------------------------------------------------------------------
%           Xnode =   [m x p] data matrix containing predicted character values
%                       for m nodes and p characters.
%

% RE Strauss, 6/23/98
%   5/22/01 - major rewrite; incorporated weighting by 1/brlen^2;
%               set negative branch lengths to zero.
%   5/23/01 - added call to ancmove().

function Xnode = contmape(anc,X,rooted,brlen)
  if (nargin < 3) rooted = []; end;
  if (nargin < 4) brlen = []; end;

  def_brlen = 0;
  if (isempty(brlen))
    brlen = ones(size(anc));
    def_brlen = 1;
  end;
  if (isempty(rooted))
    rooted = 1;
  end;

  if (length(brlen)~=length(anc))
    error('  CONTMAP: branch-length and ancestor-fn vectors incompatible.');
  end;
  if (any(brlen < 0))
    disp('  CONTMAPE warning: negative branch lengths set to zero');
    i = find(brlen < 0);
    brlen(i) = zeros(length(i),1);
  end;

  if (isvector(X))
    X = X(:);
  end;

  [T,p] = size(X);                    % Numbers of taxa and chars
  N = T-1;                            % Number of nodes
  len_anc = length(anc);
  V = T+N;                            % Number of vertices (nodes + taxa)
  if (len_anc ~= V)
    disp( '  CONTMAPE: Data matrix and ancestor function not compatible.');
    error('            Root must have ancestor 0.');
  end;

  if (~rooted)                        % If tree not rooted, remove root
    root = find(anc == 0);            % Find arbitrary root node
    desc = find(anc == root);         % Find its descendants
    d1 = desc(1);
    d2 = desc(2);
    anc(d1) = d2;                     % Form trichotomy
    if (~def_brlen)
      brlen(d1) = brlen(d1)+brlen(d2);  % Sum the branch lengths
    end;
    indx = find(anc);
    anc = anc(indx);                  % Remove root node
    len_anc = length(anc);
    brlen = brlen(indx);
    N = N-1;
    V = V-1;
  end;

  [anc,brlen] = ancmove(anc,brlen);

  M = zeros(V,V);                     % Allocate connectance matrix
  C = zeros(N,N);                     % Allocate coefficient matrix
  Y = zeros(N,p);                     % Allocate solution matrix

  w = 1./(brlen.^2);                  % Weights
  for i = 1:len_anc                   % Generate connectance matrix with recip branch lengths
    j = anc(i);
    if (j>0)
      M(i,j) = w(i);
      M(j,i) = w(i);
    end;
  end;
  Ma = M(T+1:V,1:V);                  % Connectance submatrix
  Mai = Ma(:,T+1:V);                  % Internal connectance submatrix
  Mt = Ma(:,1:T);

  [i,j] = find(Mai>0);
  for k = 1:length(i)
    C(i(k),j(k)) = -Mai(i(k),j(k));
  end;
  C = C + diag(rowsum(Ma));

  for i = 1:N
    j = find(Mt(i,:)>0);
    if (~isempty(j))
      for k = 1:length(j)
        Y(i,:) = Y(i,:) + X(j(k),:);
      end;
    end;
  end;

  Xnode = C\Y;                        % Estimate HTU character values

  return;

