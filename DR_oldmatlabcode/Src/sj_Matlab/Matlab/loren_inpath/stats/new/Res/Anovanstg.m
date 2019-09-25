% ANOVANSTG: Design matrices and auxiliary information for anovanst().
%
%     Usage: [subgrps,G,SG,n4] = anovanstg(grps,subgrps)
%
%         grps =    group-identification vector.
%         subgrps = subgroup-identification vector.
%         -----------------------------------------------------------------
%         subgrps = revised subgroup vector with unique id's across groups.
%         G =       design matrix for groups.
%         SG =      design matrix for subgroups.
%         ngrps =   number of groups.
%         n4 =      required parameter for nested anova.
%

% RE Strauss, 7/14/00 - separated from anovanst.

function [subgrps,G,SG,ngrps,n4] = anovanstg(grps,subgrps)
  G = design(grps);
  ngrps = size(G,2);

  id = 0;
  newsubgrps = zeros(size(subgrps));
  n4 = 0;

  for g = 1:ngrps                     % Ensure unique subgrp ids across groups
    gobs = find(G(:,g));
    sg = uniquef(subgrps(gobs),1);
    n4p = 0;

    for i = 1:length(sg)
      id = id+1;
      cursg = gobs(find(subgrps(gobs)==sg(i)));
      lencur = length(cursg);
      newsubgrps(cursg) = id*ones(lencur,1);
      n4p = n4p + lencur*lencur;
    end;

    n4 = n4 + n4p/length(gobs);
  end;
  subgrps = newsubgrps;

  SG = design(subgrps);

  return;


