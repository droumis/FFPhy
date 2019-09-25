% TREELEN:  Calculates the Minkowski length of a phylogenetic tree, separately for
%           each character, given the ancestor function and data matrix for all
%           taxa and nodes.
%
%       Usage: treelength = treelen(anc,X,{mink},{rooted},{brlen});
%
%           anc =        ancestor function
%           X =          [(n+m) x p] matrix of character states for n taxa
%                          and m nodes
%           mink =       Minkowski exponent (eg, 1=Manhattan, 2=Euclidean)
%                          [default = Euclidean].  If mink is a vector, then
%                          a corresponding vector of tree lengths is returned.
%           rooted =     boolean flag; if 1 (=TRUE), the root is treated as a node;
%                          if 0 (=FALSE), the tree is treated as an arbitrarily
%                          rooted network [default = 1]
%           brlen =      optional vector of branch lengths, corresponding to
%                          edges specified by the ancestor function, for
%                          calculation of tree lengths weighted by the
%                          reciprocals of branch lengths; the branch length
%                          for the root node is assumed to be zero.
%           treelength = [k x p] vector of tree lengths for each of p
%                          characters, where k is the number of Minkowski exponents
%                          passed
%

% RE Strauss, 8/11/95
%   9/20/99 - update handling of null input arguments.
%   5/26/01 - handle input data as vector.

function treelength = treelen(anc,X,mink,rooted,brlen)
  if (nargin < 3) mink = []; end;
  if (nargin < 4) rooted = []; end;
  if (nargin < 5) brlen = []; end;

  if (isvector(X))                      % Data vector should be column
    X = X(:);
  end;

  lenanc = length(anc);
  lenmink = length(mink);
  [N,P] = size(X);

  if (lenanc ~= N)
    disp('  TREELEN: Number of rows of data matrix must'); 
    error('           equal length of ancestor function');
  end;

  root = find(anc == 0);              % Find root
  if (length(root) ~= 1)
    error('  TREELEN: Ancestor function must have one root node');
  end;

  base = find(anc == root);           % Find the two basal branches
  if (length(base) ~= 2)
    error('TREELEN: Can''t find basal branches in ancestor function');
  end;

  if (isempty(mink))                 % Optional arguments
    mink = 2;                           % Default to Euclidean distance
  end;
  if (isempty(rooted))
    rooted = 1;                         % Default to rooted tree
  end;
  if (isempty(brlen))
    wt = ones(1,lenanc);                % Default to equal branch weights
  else
    if (length(brlen) ~= lenanc)
      error('  TREELEN: Ancestor and branch-length vectors must be same length');
    end;
    brlen(root) = 1;
    wt = 1./brlen;                    % Wts are reciprocals of branch lengths
    wt = (lenanc-1)*wt/sum(wt);       % Standardize weights: sum to # branches
  end;

  treelength = zeros(lenmink,P);
  for c = 1:P                         % Cycle thru characters
    basesum = 0;
    basewt = 0;
    for a=1:lenanc                      % Cycle thru ancestor function
      if (a ~= root)                      % Skip root node
        top = a;
        bot = anc(a);
        d = abs(X(top,c)-X(bot,c));

        if (rooted)
          treelength(:,c) = treelength(:,c) + (wt(a).*d.^mink(:));
        else % (not rooted)
          if (any(find(bot == base)))     % Regular branch
            treelength(:,c) = treelength(:,c) + (wt(a).*d.^mink(:));
          else                            % Examine basal branches separately
            basesum = basesum + d;
            basewt = basewt + wt(a)/2;
          end;
        end;  % if (rooted)
      end;  % if (a~=root)
    end;  % for a=1:lenanc

    if (~rooted)
      treelength(:,c) = treelength(:,c) + (basewt.*basesum.^mink(:));
    end;
  end;  % for c=1:P

  for r = 1:lenmink                     % Take roots to convert to original units
    treelength(r,:) = treelength(r,:).^(1/mink(r));
  end;

  return;

