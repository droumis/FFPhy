% COLLAPSE: Given a tree containing zero branch lengths, collapses the ancestor function 
%           into a reduced number of nodes.
%
%     Usage: [new_anc,new_brlen] = collapse(anc,brlen)
%

function [anc,brlen] = collapse(anc,brlen)


  if (all(brlen>0))
    return;
  end;

  while (any(brlen==0))
    i = find(brlen==0);
    i = i(1);
    j = anc(i);
i_j = [i j]
    anc(i) = anc(j);
    brlen(i) = brlen(j);
    anc(j) = [];
    brlen(j) = [];
  end;
  
  return;
