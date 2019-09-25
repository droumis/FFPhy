<<< UNFINISHED >>>

% TREESTAT:  Calculates tree length and pairwise patristic distances and 
%            consistency indices, given an ancestor function and character 
%            values for terminal taxa and internal nodes.
%
%     Usage: [treelength,patrist,consist] = treestat(anc,X,{rooted})
%
%           anc =        ancestor function describing tree topology
%           X =          [n+m x p] data matrix for n terminal taxa, m 
%                          internal nodes, and p characters
%           rooted =     boolean flag indicating whether tree root is (=TRUE) or
%                          is not (=FALSE) rooted; if not, the tree is treated 
%                          as an arbitrarily rooted network [default = TRUE]
%           treelength = [p x 2] vector containing the root-square (col1) and
%                          linear (absolute) (col2) tree lengths for each character
%           patrist =    [n x n] symmetric matrix of pairwise patristic
%                          distances among taxa
%           consist =    [n x n] symmetric matrix of pairwise consistency-index
%                          ratios among taxa
%

% RE Strauss, 4/6/95

function [treelength,patrist,consist] = treestat(anc,X,rooted)
  TRUE = 1; FALSE = 0;

  if (nargin < 3)
    rooted = TRUE;
  end;
  if (nargout < 2)
    calc_treelength = FALSE;
  else
    calc_treelength = TRUE;
  end;
  if (nargout < 3)
    calc_patrist = FALSE;
  else
    calc_patrist = TRUE;
  end;
  if (nargout < 4)
    calc_consist = FALSE;
  else
    calc_consist = TRUE;
  end;

  [V,n_chars] = size(X);              % Numbers of vertices and chars
  T = (V+1)/2;                        % Number of terminal taxa
  N = V-T;                            % Number of nodes
  len_anc = length(anc);
  if (len_anc ~= V)
    error('  Data matrix and ancestor function not compatible');
  end;

  return;

