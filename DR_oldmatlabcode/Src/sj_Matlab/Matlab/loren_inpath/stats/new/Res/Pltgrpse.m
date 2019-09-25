% PLTGRPSE: Plot standard-error bars parallel to axes for function plotgrps()
%

% RE Strauss, 5/31/99
%   9/2/99 - changed plot colors for Matlab v5.

function pltgrpse(x,y,g,plot_centroids,centr_labels,grpid,rng)
  n = length(x);
  xm = mean(x);
  ym = mean(y);

  xse = std(x)/sqrt(n);
  yse = std(y)/sqrt(n);

  critval = abs(tinv(0.025,n-1));
  xl = xm - critval*xse;
  xu = xm + critval*xse;
  yl = ym - critval*yse;
  yu = ym + critval*yse;

  d = 0.02;
  rx = range(x);
  ry = range(y);
  xdl = xm - d*rx;
  xdu = xm + d*rx;
  ydl = ym - d*ry;
  ydu = ym + d*ry;

  plot([xm xm],[yl yu],'k');
  plot([xl xu],[ym ym],'k');
  plot([xdl xdu],[yl yl],'k');
  plot([xdl xdu],[yu yu],'k');
  plot([xl xl],[ydl ydu],'k');
  plot([xu xu],[ydl ydu],'k');

  if (centr_labels & ~plot_centroids)
    xpos = xm + 0.015 * rng(1);
    ypos = ym + 0.040 * rng(2);
    if (isempty(grpid))
      text(xpos,ypos,num2str(g));
    else
      text(xpos,ypos,grpid(g,:));
    end;
  end;

  return;