% PLTGRPCH: Plots convex hulls for function plotgrps()
%

% RE Strauss, 5/31/99
%   9/2/99 - changed plot colors for Matlab v5.

function pltgrpch(x,y,g,grp_labels,grpid,plot_points,symbol,delta)
  hull_pts = hull([x y]);
  plot(hull_pts(:,1),hull_pts(:,2),'k');

  if (grp_labels)                     % Print group labels
    [h,i] = max(hull_pts(:,1));         % Find right-most hull pt
    xpos = hull_pts(i,1) + delta;
    ypos = hull_pts(i,2);
    if (isempty(grpid))
      text(xpos,ypos,num2str(g));
    else
      text(xpos,ypos,grpid(g,:));
    end;
  end;

  return;



