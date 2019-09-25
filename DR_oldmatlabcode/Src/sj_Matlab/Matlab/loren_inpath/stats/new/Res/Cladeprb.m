% CLADEPRB: Calculates the proportion of total rooted trees for T labeled 
%           terminal taxa (excluding an outgroup, if used to root the tree) 
%           that contain a given labeled clade of t taxa.  
%           This proportion can be used as the probability of drawing a 
%           random tree containing the given clade.
%
%     Usage: [p,t,nw,no,nc,N] = cladeprb(t,T)
%
%             t =  number of labeled taxa in given clade [scalar or vector].
%             T =  number of labeled terminal taxa [scalar or vector], 
%                    range 2 to T-1.
%             -------------------------------------------------------------------
%             p =  proportion of possible trees for T taxa containing the 
%                    given clade.
%             t =  number of labeled taxa in given clade.
%             nw = number of possible rearrangements of t taxa within the clade.
%             no = number of possible trees containing the given clade;
%                    =number of possible rearrangements of the T-t taxa 
%                    outside of the given clade.
%             nc = total number of possible trees containing the given clade.
%             N =  total number of possible rooted trees for T-1 taxa.
%

% Rohlf 1983; Edwards 1964; Cavalli-Sforza & Edwards 1967

% RE Strauss, 8/17/98

function [p,t,nw,no,nc,N] = cladeprb(t,T)
  t = t(:);
  T = T(:);

  nt = length(t);
  nT = length(T);

  if (nt>1 & nT==1)
    nT = nT * ones(size(nt));
  elseif (nt==1 & nT>1)
    nt = nt * ones(size(nT));
  end;

  if (any(t<2) | any(t>(T-1)))
    disp('  CLADEPRB warning: clade size (t) out of range');
    i = find(t < 2);
    if (~isempty(i))
      t(i) = [];
    end;
    i = find(t > T-1);
    if (~isempty(i))
      t(i) = [];
    end;
  end;

  T = T+1;                % Include root as terminal taxon
  nw = treenum(t);
  no = treenum(T-t);
  nc = nw.*no;
  N = treenum(T-1);
  p = nc./N;

  return;
