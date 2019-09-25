% MSTINFL: Finds the influence of individual data points on an MST 
%          k-groups cluster analysis, based on the Rand index (using the 
%          Hubert & Arabie correction), HA.
%             Because the Hubert & Arabie statistic is used here as a measure 
%          of internal influence, it is expressed here as ha=1-HA, which 
%          varies from 0 for a point having no influence on the partitioning, 
%          to 1+ for a point having an extreme impact on the partitioning.
%
%          Cheng, R. & G.W. Milligan. 1996. Measuring the influence of 
%             individual data points in a cluster analysis.  
%             J.Classif. 13:315-335.
%
%     Syntax: ha = mstinfl(crds,ngrps)
%
%         crds -   [n x p] matrix of point coordinates.
%         ngrps -  vector (of length m) containing number of groups 
%                     for which influences are desired 
%                     [default = 2:floor(n/2)].
%         ------------------------------------------------------------
%         ha -     [n x m] matrix of Hubert & Arabie corrected Rand 
%                     indices.
%

% RE Strauss, 12/31/96
%   1/25/00 - misc but minor improvements.

function ha = mstinfl(crds,ngrps)
  if (nargin < 2) ngrps = []; end;

  [n,p] = size(crds);

  if (isempty(ngrps))
    ngrps = [2:floor(n/2)]';
  end;

  m = length(ngrps);                    % Allocate output matrices
  ha = zeros(n,m);

  for mm = 1:m                          % Cycle thru number of groups
    ng = ngrps(mm);                       % Current number of groups
    rgrps = mstgrp(crds,ng);              % Reference clustering

    for obs = 1:n                           % Remove each obs in turn
      xcrds = crds;                         % from coord matrix
      xcrds(obs,:) = [];
      refgrps = rgrps;                      % and from reference clustering
      refgrps(obs) = [];
      iref = uniquef(refgrps);              % Reference groups
      nref = length(iref);

      redgrps = mstgrp(xcrds,ng);           % Clustering of reduced obs set
      ired = uniquef(redgrps);              % Reduced groups
      nred = length(iref);

      clst = zeros(nref,nred);              % Construct partitions matrix
      for i = 1:nref                          % Reference clusters
        for j = 1:nred                        % Reduced clusters
          clst(i,j) = sum(refgrps==i & redgrps==j);
        end;
      end;
      colsum = sum(clst);
      rowsum = sum(clst')';
      totsum = sum(colsum);

      term1 = 0; term2 = 0; term3 = 0;      % Calc H&A Rand index for the 
      for i = 1:nref                        %   two partitions
        term2 = term2 + comb(rowsum(i),2);
        for j = 1:nred
          term1 = term1 + comb(clst(i,j),2);
        end;
      end;
      for j = 1:nred
        term3 = term3 + comb(colsum(j),2);
      end;
      term4 = comb(totsum,2);

      c = term2*term3/term4;
      ha(obs,mm) = 1 - ((term1-c)/(0.5*(term2+term3)-c));
    end;
  end;
  
  return;

