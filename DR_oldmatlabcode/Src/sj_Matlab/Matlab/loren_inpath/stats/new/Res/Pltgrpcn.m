% PLTGRPCN: Plots group centroids for function plotgrps()
%

% RE Strauss, 5/31/99
%   9/3/99 - changed plot color for Matlab v5.

function pltgrpcn(x,y,g,centr_labels,grpid,rng)
  if (length(x)>1)
    m = mean([x y]);
  else
    m = [x y];
  end;
  plot(m(1),m(2),'k+');

  if (centr_labels)                     % Print centroid labels
    xpos = m(1) + 0.015 * rng(1);
    ypos = m(2) + 0.040 * rng(2);
    if (isempty(grpid))
      text(xpos,ypos,num2str(g));
    else
      text(xpos,ypos,grpid(g,:));
    end;
  end;

  return;


