% PLTGRPMD: Plots group medians (via hull peeling) for function plotgrps()
%

% RE Strauss, 5/31/99
%   9/2/99 -  changed plot colors for Matlab v5.
%   5/20/00 - changed output arguments of hullpeel.

function pltgrpcn(x,y,g,centr_labels,grpid,rng)
  if (length(x)>1)
    [pp,hp,ha,m] = hullpeel([x y],0);
  else
    m = [x y];
  end;
  plot(m(1),m(2),'kx');

  if (centr_labels)                     % Print mediod labels
    xpos = m(1) + 0.015 * rng(1);
    ypos = m(2) + 0.040 * rng(2);
    if (isempty(grpid))
      text(xpos,ypos,num2str(g));
    else
      text(xpos,ypos,grpid(g,:));
    end;
  end;

  return;
