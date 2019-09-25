% STEPMOVA: Stepwise MANOVA to determine the best subset of variables
%
%     Usage: [incl,F] = stepmova(X,grps)
%
%        X =    [n x p] data matrix (obs x vars).
%        grps = row or column vector of group identifiers.
%        ---------------------------------------------------------------------
%        incl = list of indices of variables in order of inclusion.
%        F =    corresponding F statistics.
%

% RE Strauss, 4/28/98
%   11/29/99 - changed calling sequence.

function [incl,F] = stepdiscX,grps)
  [n,p] = size(X);

  vars_incl = zeros(1,p);           % Positions in which vars introduced
  Fmax = zeros(1,p);                % F statistic for that series of vars

  for step = 1:p
    vi = find(vars_incl>0);             % Indices of vars included
    vni = find(vars_incl==0);           % Indices of vars not included
    nvni = length(vni);                 % Number of vars not included

    F = zeros(nvni,1);
    for v = 1:nvni
      [lambda,F(v),pr,df] = manova(X(:,[vi vni(v)]),grps);
    end;

    [Fm,i] = max(F);
    Fmax(vni(i)) = Fm;
    vars_incl(vni(i)) = step;
  end;
vF = [vars_incl' Fmax']

  [y,i] = sort(vars_incl);
  incl = 1:p;
  incl = incl(i);
  F =    Fmax(i);

  return;
