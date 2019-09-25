% PLTGRPMA: Plot major-axes or ellipsebound(es for function plotgrps()
%

% RE Strauss, 5/31/99
%   9/2/99 - changed plot colors for Matlab v5.
%   4/6/02 -  change ellips() to ellipsebound().

function XY = pltgrpma(x,y,XY,plot_majaxes,plot_cis,plot_ellipses, ...
                        plot_regr,g,plot_hulls,grp_labels,grpid)
  b = majaxis(x,y);
  slope = b(1,1);
  intcpt = b(1,2);

  if (plot_cis | plot_ellipses)
    scores1 = project(x,y,slope,0);     % Project onto major axis
    scores2 = project(x,y,-1/slope,0);  % Project onto minor axis
    h = mean(x);
    k = mean(y);
    theta = atan(slope);

    n = length(x);
    critval = abs(tinv(0.025,n-1));
  end;

  if (plot_majaxes)                     % Plot major axes
    xmargin = 0.05*range(x);
    ymargin = 0.05*range(y);

    xmin = min(x);
    xmax = max(x);
    ymin = min(y);
    ymax = max(y);

    xlow =  xmin - xmargin;
    xhigh = xmax + xmargin;

    ylow = xlow * slope + intcpt;
    yhigh = xhigh * slope + intcpt;
        
    if (yhigh > ymax+ymargin)
      yhigh = ymax + ymargin;
      xhigh = (yhigh - intcpt) / slope;
    end;

    if (ylow < ymin-ymargin)
      ylow = ymin - ymargin;
      xlow = (ylow - intcpt) / slope;
    end;

    if (plot_regr)
      plot([xlow,xhigh],[ylow,yhigh],':k');
    else
      plot([xlow,xhigh],[ylow,yhigh],'k');
    end;

    XY = [XY; xlow ylow; xhigh yhigh];  % Extend axis ranges    
  end;

  if (plot_cis)                         % Plot centroid confidence ellipses
    a = critval*std(scores1)/sqrt(n);
    b = critval*std(scores2)/sqrt(n);
    [xx,yy] = ellipsebound((a,b,h,k,theta);
    plot(xx,yy,'k');
  end;

  if (plot_ellipses)                    % Plot data ellipses
    a = critval*std(scores1);
    b = critval*std(scores2);
    [xx,yy] = ellipsebound((a,b,h,k,theta);
    plot(xx,yy,'k');

    if (grp_labels & ~plot_hulls)
      delta = 0.02*range(x);

      [h,i] = max(xx);                  % Find right-most ellipse pt
      xpos = xx(i) + delta;
      ypos = yy(i);
      if (isempty(grpid))
        text(xpos,ypos,num2str(g));
      else
        text(xpos,ypos,grpid(g,:));
     end;
    end;

    XY = [XY; min(xx) min(yy); max(xx) max(yy)]; % Extend axis ranges
  end;

  return;
