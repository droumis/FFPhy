% TREENUM:  Calculates the numbers of possible rooted and unrooted binary
%           trees for T labelled terminal (pendant) taxa.
%
%     Usage: [r,u] = treenum(T)
%
%           T = number of terminal taxa [scalar or vector].
%           -----------------------------------------------
%           r = number of rooted-tree topologies.
%           u = number of unrooted-tree topologies.
%

% Rohlf, 1983 (rooted)
% Cavalli-Sforza and Edwards, 1967 (unrooted)

function [r,u] = treenum(T)
  T = T(:);
  nT = length(T);

  if (nT == 1)
    r = prod(1:2:(2*T-3));
    u = prod(1:2:(2*T-5));
  else
    r = zeros(size(T));
    u = zeros(size(T));
    for i = 1:nT
      r(i) = prod(1:2:(2*T(i)-3));
      u(i) = prod(1:2:(2*T(i)-5));
    end;
  end;

  return;
