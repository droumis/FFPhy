% TREETIPS: Given an ancestor function, returns a vector of the terminal taxa 
%           (tips) contained within a given clade, as specified by its 
%           ancestral node.
%
%     Usage: tips = treetips(anc,node)
%
%           anc =  ancestor function, with ancestor of root node specified as 0.
%           node = the node identifying the clade.
%           -----------------------------------------------------------------
%           tips = column vector containing a list of the terminal taxa 
%                    contained within the specified clade; set to null if
%                    'node' is a terminal taxon.
%

function tips = treetips(anc,node)
  ntaxa = (length(anc)+1)/2;

  tips = [];
  if (node > ntaxa)
    for t = 1:ntaxa
      c = t;
      while (anc(c)>0)
        if (anc(c)==node)
          tips = [tips; t];
        end;
        c = anc(c);
      end;
    end;
  end;

  return;

