% CORRPLOT: Produces a matrix of plots reflecting the structure of a
%           correlation matrix, with bivariate plots on the off-diagonal
%           and univariate histograms on the diagonal.  Bivariate plots
%           can optionally portray centroids, hulls, major axes, or
%           confidence ellipses of data or centroids.
%
%     Syntax:  corrplot(x,{option})
%
%           x -      [n x p] data matrix for n observations and p variables
%           option - optional row vector of boolean flags specifying options,
%                    in any combination:
%                       1) preserve scaling units for all variables
%                            (default: axes are individually scaled)
%                       2) plot centroids
%                       3) plot convex hulls
%                       4) plot major axes
%                       5) plot 95% confidence ellipses of centroids
%                       6) plot 95% confidence ellipses of data
%                       7) square subplots (default: rectangular)
%

% RE Strauss, 10/23/97
%   8/20/99 - changed plot colors and other characteristics for Matlab v5.
%   2/29/00 - use histgram() function for histograms.
%   4/6/02 -  change ellips() to ellipsebound().
%   9/9/02 -  corrected syntax errors in plot_cis and plot_ellipses blocks.

function corrplot(x,option)
  n_option = 7;

  if (nargin < 2)
    option = zeros(1,n_option);
  elseif (length(option) < n_option)
    option = [option, zeros(1,n_option-length(option))];
  end;

  if (option(1))
    preserve_scaling = 1;
  else
    preserve_scaling = 0;
  end;
  if (option(2))
    plot_centroids = 1;
  else
    plot_centroids = 0;
  end;
  if (option(3))
    plot_hulls = 1;
  else
    plot_hulls = 0;
  end;
  if (option(4))
    plot_majaxes = 1;
  else
    plot_majaxes = 0;
  end;
  if (option(5))
    plot_cis = 1;
  else
    plot_cis = 0;
  end;
  if (option(6))
    plot_ellipses = 1;
  else
    plot_ellipses = 0;
  end;
  if (option(7))
    square_subplots = 1;
  else
    square_subplots = 0;
  end;

  [n,p]=size(x);

  meanx = mean(x);                        % Basic statistics
  stdx = std(x);
  minx = min(x);
  maxx = max(x);

  minminx = min(minx);
  maxmaxx = max(maxx);

  lowpt = zeros(1,2);                     % Pts to plot regression line
  highpt = zeros(1,2);
  l = zeros(1,2);
  u = zeros(1,2);
  lbound = zeros(1,2);
  ubound = zeros(1,2);

  figure;
  for i=1:p
    for j=1:p
      subplot (p,p,((i-1)*p+j));

      % Histogram
      if (i==j)
        [n,f] = histgram(x(:,i),[],[],[],[],1);
%        axis('square');
        putylab(' ');
%        if (preserve_scaling)             % Axis ranges
%          gap = 0.10*(maxmaxx-minminx);
%          maxf = max(f(:,2));
%          axrng = [minminx-gap maxmaxx+gap 0 maxf+(0.1*maxf)];
%          axis(axrng);                      % Set axes
%        end;

      % Bivariate plot
      else
        hold on;
        plot(x(:,i),x(:,j),'.k');         % Plot points

        if (plot_centroids)
          plot(meanx(i),meanx(j),'+k');   % Plot centroid
        end;

        if (plot_hulls)                   % Plot convex hull
          hull_pts = hull([x(:,i),x(:,j)]);
          plot(hull_pts(:,1),hull_pts(:,2),'k');
        end;

        if (plot_majaxes | plot_cis | plot_ellipses)
          b = majaxis(x(:,i),x(:,j));     % Major axis
          slope = b(1,1);
          intcpt = b(1,2);
          if (plot_cis | plot_ellipses)
            scores1 = project(x(:,i),x(:,j),slope,0);    % Project onto major axis
            scores2 = project(x(:,i),x(:,j),-1/slope,0); % Project onto minor axis
            h = meanx(i);
            k = meanx(j);
            theta = atan(slope);
          end;
        end;

        if (plot_cis)
          a = 1.96*std(scores1)/sqrt(n);
          b = 1.96*std(scores2)/sqrt(n);
          [xx,yy] = ellipsebound(a,b,h,k,theta);  % 95% confidence ellipse
          plot(xx,yy,'k');
        end;

        if (plot_ellipses)
          a = 1.96*std(scores1);
          b = 1.96*std(scores2);
          [xx,yy] = ellipsebound(a,b,h,k,theta);  % 95% confidence ellipse
          plot(xx,yy,'k');

          l(1) = min(minx(i),min(xx));
          u(1) = max(maxx(i),max(xx));
          l(2) = min(minx(j),min(yy));
          u(2) = max(maxx(j),max(yy));
        else
          l(1) = minx(i);
          u(1) = maxx(i);
          l(2) = minx(j);
          u(2) = maxx(j);
        end;

        gap = 0.10*(u-l);                 % Leave gap between data and axes
        lbound = l - gap;                 %   on bivariate plots
        ubound = u + gap;

        if (plot_majaxes)
          lowpt(1) = lbound(1);           % Draw major axis
          lowpt(2) = lbound(1)*slope + intcpt;
          highpt(1) = ubound(1);
          highpt(2) = ubound(1)*slope + intcpt;
          plot([lowpt(1) highpt(1)],[lowpt(2) highpt(2)],'k');
        end;

        axis([lbound(1),ubound(1),lbound(2),ubound(2)]);
        if (preserve_scaling)
          axis('equal');
        end;
        hold off;
      end;

    if (square_subplots)
      axis('square');
    end;
    set(gca,'box','on');
    set(gca,'XTickLabel', []);
    set(gca,'YTickLabel', []);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    end;
  end;
  return;


