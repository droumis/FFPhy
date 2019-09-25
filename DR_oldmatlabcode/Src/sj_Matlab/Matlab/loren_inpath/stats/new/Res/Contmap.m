% CONTMAP:  Least-squares mapping of continuous characters onto a tree.
%           Uses analytic solution of McArdle & Rodrigo (1994) for a least-squares 
%           fitting of a completely dichotomous tree, or iterative optimization 
%           for trees containing zero branch lengths or for Minkowski coefficient 
%           other than 2.
%
%     Usage: Xnode = contmap(anc,X,{rooted},{brlen},{k},{useit})
%
%           anc =     ancestor function describing tree topology; the 
%                       root must have ancestor 0.
%           X =       [n x p] data matrix for n taxa and p characters.
%           rooted =  optional boolean flag indicating whether tree root is (=1) 
%                       or is not (=0) to be optimized; if not, the tree is
%                       treated as an arbitrarily rooted network [default=1].
%           brlen =   optional branch lengths for weighted least-squares 
%                       mapping, in same sequence as anc [default = 1 for all 
%                       branches except ancestor].
%           k =       optional value for Minkowski k [default=2 for least-squares 
%                       solution].
%           useit =   optional boolean flag indicating that the iterative 
%                       algorithm is to be used in all cases [default=0].
%           ---------------------------------------------------------------------
%           Xnode =   [m x p] data matrix containing predicted character values
%                       for m nodes and p characters.
%

% RE Strauss - 7/22/95
%   6/23/98 -  added McArdle & Rodrigo solution; rewrite into separate functions.
%   12/17/99 - added Minkowski k.
%   5/22/01 -  set negative branch lengths to zero.
%   5/23/01 -  add option of using iterative algorithm in all cases;
%                added call to ancmove().

function Xnode = contmap(anc,X,rooted,brlen,k,useit)
  if (nargin < 3) rooted = []; end;
  if (nargin < 4) brlen = []; end;
  if (nargin < 5) k = []; end;
  if (nargin < 6) useit = []; end;
  
  if (isempty(rooted)) rooted = 1; end;
  if (isempty(k))      k = 2; end;

  if (isempty(brlen))
    brlen = ones(size(anc));
  end;
  if (length(brlen)~=length(anc))
    error('  CONTMAP: branch-length and ancestor-fn vectors incompatible.');
  end;
  if (any(brlen < 0))
    disp('  CONTMAP warning: negative branch lengths set to zero');
    i = find(brlen < 0);
    brlen(i) = zeros(length(i),1);
  end;

  exact_soln = 0;
  if (all(brlen>0) & k==2 & ~useit)
    exact_soln = 1;
  end;

  [T,n_chars] = size(X);              % Numbers of taxa and chars
  if (T==1)                           % If input is row vector, transpose
    X = X';
  end;

  if (exact_soln)
    Xnode = contmape(anc,X,rooted,brlen);     % Exact least-squares solution
  else
    Xnode = contmapi(anc,X,rooted,brlen,k);   % Iterative solution
  end;

  return;

