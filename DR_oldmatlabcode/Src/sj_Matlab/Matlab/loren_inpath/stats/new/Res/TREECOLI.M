% TREECOLI: Calculates the raw and normalized forms of Colless' I, a measure
%           of tree asymmetry (imbalance), from the ancestor function.
%
%     Usage:  [I,In] = treecoli(anc)
%
%           anc = vector specifying ancestor function, with ancestor of 
%                   root specified as 0.
%           -----------------------------------------------------------
%           I =   Colless' I (raw).
%           In =  colless' I (normalized).
%

% Rogers, JS. 1993. Response of Colless's tree imbalance to number of terminal taxa.  
%   Syst. Biol. 42:102-105.

% RE Strauss, 8/2/98

function [I,In] = treecoli(anc)
  lenanc = length(anc);
  n_taxa = ceil(lenanc/2);            % Number of taxa
  n_anc = n_taxa-1;                   % Number of nodes
  tally = zeros(n_anc,2);

  while (any(anc))
    for i=1:lenanc
      a = anc(i)-n_taxa;
      if (a>0)
        if (tally(a,1)==0)
          j = 1;                        % 1st set of descendants for node
        else
          j = 2;                        % 2nd set of descendants for node
        end;

        if (i<=n_taxa)
          tally(a,j) = 1;                  % Add 1 for terminal taxa
          anc(i) = 0;
        elseif (all(tally(i-n_taxa,:)))
          tally(a,j) = sum(tally(i-n_taxa,:));  % Add #descendants for nodes
          anc(i) = 0;
        end;
      end;
    end;
  end;

  I = sum(abs(tally(:,1)-tally(:,2)));    % Raw

  Imax = (n_taxa-1)*(n_taxa-2)/2;
  In = I / Imax;                          % Normalized 0-1

  return;
