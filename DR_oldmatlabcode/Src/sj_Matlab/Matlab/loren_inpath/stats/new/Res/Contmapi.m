% CONTMAPI: Iterative mapping of continuous characters onto a tree, which 
%           may be rooted or unrooted and may contain polytomies 
%           (zero branch lengths).
%           Called by function contmap(), but can be used stand-alone.
%
%     Usage: Xnode = contmapi(anc,X,rooted,brlen,k)
%
%           anc =     ancestor function describing tree topology; the 
%                       root must have ancestor 0.
%           X =       [n x p] data matrix for n taxa and p characters.
%           rooted =  boolean flag indicating whether tree root is (=1) or
%                       is not (=0) to be optimized; if not, the tree is
%                       treated as an arbitrarily rooted network [default = 1].
%           brlen =   optional branch lengths for weighted least-squares 
%                       mapping, in same sequence as anc [default = 1 for all 
%                       branches except ancestor].
%           k =       Minkowski k [default=2 for least-squares solution].
%           ---------------------------------------------------------------------
%           Xnode =   [m x p] data matrix containing predicted character values
%                       for m nodes and p characters.
%           treelen = final tree length (in units of weighted squared character 
%                       values).
%

% RE Strauss, 6/23/98
%   12/17/99 -  added Minkowski k
%   1/4/00 -    changed fminu() to fmins().
%   5/22/01 -   changed fmins() to fminsearch();
%               set negative branch lengths to zero.
%   5/23/01 - added call to ancmove().

function Xnode = contmapi(anc,X,rooted,brlen,k)
  if (nargin < 3) rooted = []; end;
  if (nargin < 4) brlen = []; end;
  if (nargin < 5) k = []; end;

  def_brlen = 0;
  if (isempty(brlen))
    brlen = ones(size(anc));
    def_brlen = 1;
  end;
  if (isempty(rooted))
    rooted = 1;
  end;
  if (isempty(k))
    k = 2;
  end;

  if (length(brlen)~=length(anc))
    error('  CONTMAPI: branch-length and ancestor-fn vectors incompatible.');
  end;
  if (any(brlen < 0))
    disp('  CONTMAPI warning: negative branch lengths set to zero');
    i = find(brlen < 0);
    brlen(i) = zeros(length(i),1);
  end;

  if (isvector(X))
    X = X(:);
  end;

  [T,p] = size(X);                    % Numbers of taxa and chars
  m = T-1;                            % Number of nodes
  len_anc = length(anc);
  V = T+m;                            % Number of vertices
  if (len_anc ~= V)
    disp( '  CONTMAPI: Data matrix and ancestor function not compatible');
    error('            Root must have ancestor 0');
  end;

  Xnode = zeros(m,p);                 % Allocate results matrix
  for i = 1:m                         % Initialize node values
    tips = treetips(anc,T+i);           % Find taxa within node
    Xnode(i,:) = mean(X(tips,:));       % Find mean of taxon values
  end;

  rootnode = find(anc==0);            % Find root
  [links,anc,curr_node,level] = treedivd([],anc,rootnode,0);  % Get list of branches
  links = links(:,1:2);               % Save ancestors & descendants
  brlen = brlen(links(:,2))';         % Put branch lengths in same sequence
 
  if (~rooted)                        % If tree only rooted for convenience,
    a = links(:,1);                     % Isolate ancestors & descendants of branches
    d = links(:,2);
    nlinks = size(links(1));
    i = find(a==rootnode);              % Find the two branches from root
    d1 = d(i(1));                       % Save identities of descendants
    d2 = d(i(2));
    b1 = brlen(i(1));                   % Save branch lengths
    b2 = brlen(i(2));
    links = [links; d1 d2];             % Add new link to list connecting descendants
    brlen = [brlen; b1+b2];             % Sum the two branches from root
    links(i,:) = [];                    % Remove two branches from root
    brlen(i) = [];
  end;

%  [anc,brlen] = ancmove(anc,brlen);

  for c = 1:p                         % Fit chars one at a time
    x = X(:,c);                         % Isolate taxon & node values
    xnode = Xnode(:,c);                 % Initial guesses at nodes
    options = [];
    xnode = fminsearch('contmapf',xnode,optimset('Display','off'),x,brlen,k,links);  % Optimize nodes

    if (rooted)                         % If tree rooted,
      Xnode(:,c) = xnode;                 % Stash final node estimates
    else                                % else if unrooted,
      rootnode = rootnode - T;
      xnode(rootnode) = [];               % Remove root node
      Xnode(1:(m-1),c) = xnode;           % Stash final node estimates
    end;
  end;

  if (~rooted)                        % If tree not rooted,
    Xnode(m,:) = [];                  %   remove last line of node values
  end;

  return;

