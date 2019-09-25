% CONTMAPD:  Mapping of distributions by randomization, given observed 
%            terminal data distributions.  Only a single character can be 
%            mapped.
%
%     Usage: Xnode = contmap(anc,X,{rooted},{brlen})
%
%           anc =     ancestor function describing tree topology; the 
%                       root must have ancestor 0.
%           X =       [n x k] data matrix for n taxa and k terminal taxa.  
%                       If n differs for different taxa, NaN's can be used 
%                       to complete the rectangular matrix.
%           rooted =  boolean flag indicating whether tree root is (=TRUE) or
%                       is not (=FALSE) to be optimized; if not, the tree is
%                       treated as an arbitrarily rooted network 
%                       [default = TRUE].
%           brlen =   optional branch lengths for weighted least-squares mapping.
%           ---------------------------------------------------------------------
%           Xnode =   [m x p] data matrix containing predicted distributions
%                       for m nodes and p characters.
%

% RE Strauss, 4/26/97

<<< Incomplete >>>

function Xnode = contmapd(anc,X,rooted,brlen)
  if (nargin < 3) rooted = []; end;
  if (nargin < 4) brlen = []; end;

  if (isempty(rooted))
    rooted = 1;
  end;
  if (isempty(brlen))
    weighted = 0;
    brlen = ones(length(anc));
  else
    weighted = 1;
    if (~isempty(find(brlen <= 0)))
      error('  Cannot handle zero or negative branch lengths');
    end;
  end;

  [T,n_chars] = size(X);              % Numbers of taxa and chars
  if (T==1)                           % If input is row vector, transpose
    X = X';
    t = T;
    T = n_chars;
    n_chars = t;
  end;

  N = T-1;                            % Number of nodes
  len_anc = length(anc);
  V = T+N;                            % Number of vertices
  if (len_anc ~= V)
    disp('  Data matrix and ancestor function not compatible');
    error('  Root must have ancestor 0');
  end;

  if (~rooted)                        % If tree not rooted, remove root
    root = find(anc == 0);            % Find arbitrary root node
    desc = find(anc == root);         % Find its descendants
    d1 = desc(1);
    d2 = desc(2);
    anc(d1) = d2;                     % Form trichotomy
    if (weighted)
      brlen(d1) = brlen(d1)+brlen(d2);  % Sum the branch lengths
    end;
    indx = find(anc);
    anc = anc(indx);                  % Remove root node
    brlen = brlen(indx);
    N = N-1;
    V = V-1;
  end;

  Y = zeros(N,n_chars);               % Allocate solution vector
  c = zeros(N);                       % Allocate coefficient matrix
  d = zeros(N,1);                     % Allocate diagonal

  for n = 1:N                         % Cycle thru internal nodes
    for v = 1:V                       % Cycle thru all vertices
      if (n+T == anc(v) | anc(n+T) == v);
        c(n,v) = 1/min([brlen(n),brlen(v)]);
        d(n) = sum(c(n,:));
        if (v <= T)
          Y(n,:) = Y(n,:) + X(v,:);
        end;
      end;
    end;
  end;

  C = -c((1:N),(T+1:V)) + diag(d);    % Create coefficient matrix
  Xnode = C\Y;                        % Estimate HTU character values

  return;

