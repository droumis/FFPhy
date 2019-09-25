% TREEASYM: For a given number of taxa T, calculates the Colless I asymmetry
%           values corresponding to a set of trees specified by Rohlf's M.
%
%     Usage: [I,In] = treeasym(T,M)
%
%           T =  number of terminal taxa.
%           M =  vector (length t) of Rohlf's topology numbers, varying in value 
%                  from 0 to treenum(T)-1.
%           -------------------------------------------------------------------
%           I =  Colless' I (raw).
%           In = Colless' I (normalized). 
%

% RE Strauss, 8/2/98

function [I,In] = treeasym(T,M)
  n_tree = treenum(T);
  if (max(M) > n_tree)
    error('TREEASYM: value of Rohlfs M out of bounds.');
  end;

  I =  zeros(size(M));
  In = zeros(size(M));

  for i=1:length(M)
    anc = treemanc(T,M(i));
    [I(i),In(i)] = treecoli(anc);
  end;

  return;
