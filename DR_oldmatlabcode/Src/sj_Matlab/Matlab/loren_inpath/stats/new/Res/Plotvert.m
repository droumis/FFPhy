% PLOTVERT: For bivariate data with multiple values of y for each x, plots the 
%           data as a series of vertical lines indicating the ranges for each x, 
%           with consecutive means or medians connected.  Does not allow for missing data.
%
%     Usage: plotvert(x,y,{center},{yextent},{xtransform},{ytransform},{'color'},{lineonly})
%
%           x,y = matching vectors of x,y coordinates.
%           center =  optional flag indicating the measure of center:
%                         0 = mean [default]
%                         1 = median
%           yextent = optional flag indicating the vertical extent of line  
%                       at each value of x:
%                         0 = center only [default]
%                         1 = one standard deviation in each direction
%                         2 = two standard deviations in each direction
%                         3 = central 50% of distribution
%                         4 = central 95% of distribution
%                         5 = range
%                         6 = standard error of mean or median
%                         7 = 95% CI of mean or median
%           xtransform,ytransform = optional flags indicating that x or y data 
%                       are to be transformed after calculating the y-statistics 
%                       but before plotting:
%                         0 = no transformation [default]
%                         1 = log(x)
%                         2 = log10(x)
%                         3 = exp(x)
%                         4 = 10.^x
%           color =     optional line color [default = 'k'].
%           lineonly =  optional boolean variable indicating, if true, that the 
%                         color supplied is to be used only for the connecting 
%                         lines, not the vertical bars [default = 0].
%

% RE Strauss, 4/16/98
%   7/1/00 -  substituted switch for elseif statements; 
%             added tips to vertical lines.
%   6/3/01 -  added 'lineonly' option.
%   5/30/02 - added 'center' option; 
%             added options for SE's and CI's of mean or median.

function plotvert(x,y,center,yextent,xtf,ytf,color,lineonly)
  if (nargin < 3) center = []; end;
  if (nargin < 4) yextent = []; end;
  if (nargin < 5) xtf = []; end;
  if (nargin < 6) ytf = []; end;
  if (nargin < 7) color = []; end;
  if (nargin < 8) lineonly = []; end;

  if (isempty(center))
    center = 0;
  end;
  if (isempty(lineonly))
    lineonly = 0;
  end;
  if (isempty(yextent))
    yextent = 0;
  end;
  if (isempty(xtf))
    xtf = 0;
  end;
  if (isempty(ytf))
    ytf = 0;
  end;

  if (isempty(color))                 
    colorv = 'k';
    colorh = 'k';
  else
    colorh = color;
    colorv = color;
    if (lineonly)
      colorv = 'k';
    end;
  end;

  if (~isvector(x) | ~isvector(y))
    error('  PLOTVERT: input vectors must be vectors.');
  end;
  x = x(:);
  y = y(:);

  if (size(x) ~= size(y))
    error('  PLOTVERT: input vectors must match in size');
  end;

  if (xtf==1 | xtf==2)
    if any(x<0)
      disp('  PLOTVERT warning: some x are negative');
    end;
  end;
  if (ytf==1 | ytf==2)
    if any(y<0)
      disp('  PLOTVERT warning: some y are negative');
    end;
  end;

  xunique = uniquef(x,1);
  nx = length(xunique);

  ystats = zeros(nx,3);

  for i = 1:nx
    yy = y(x==xunique(i));
    
    switch(center)
      case 0,
        m = mean(yy);
        
      case 1,
        m = median(yy);
        
      otherwise,
        error('  PLOTVERT: invalid value for center');
    end;

    switch(yextent)
      case 0,
        ystats(i,:) = [m m m];

      case 1,
        s = std(yy);
        ystats(i,:) = [m-s m m+s];

      case 2,
        s = std(yy);
        ystats(i,:) = [m-2*s m m+2*s];

      case 3,
        lims = centdist(yy,50);
        ystats(i,:) = [lims(1) m lims(2)];

      case 4,
        lims = centdist(yy,95);
        ystats(i,:) = [lims(1) m lims(2)];

      case 5,
        ystats(i,:) = [min(yy) m max(yy)];
        
      case 6,
        se = std(yy)/sqrt(length(yy));
        if (center==1)
          se = 1.2533*se;
        end;
        ystats(i,:) = [m-se m m+se];
        
      case 7,
        
        if (center==0)
          se = std(yy)/sqrt(length(yy));
          ci = 1.96*se;
          ystats(i,:) = [m-ci m m+ci];
        else
          [m,ci,rcl] = medianci(yy,0.95);
          ystats(i,:) = [ci(1) m ci(2)];
        end;

      otherwise,
        error('  PLOTVERT: invalid value for yextent');
    end;
  end;

  switch(xtf)                           % Transformations of x
    case 0,
      xvals = xunique;

    case 1,
      xvals = NaN*ones(size(xunique));
      pos = (xunique>0);
      xvals(pos) = log(xunique(pos));

    case 2,
      xvals = NaN*ones(size(xunique));
      pos = (xunique>0);
      xvals(pos) = log10(xunique(pos));

    case 3,
      xvals = exp(xunique);

    case 4,
      xvals = 10.^(xunique);

    otherwise,
      error('  PLOTVERT: invalid valaue for xtransform');
  end;

  switch(ytf)                           % Transformations of y
    case 0,
      ys = ystats;

    case 1,
      neg = any(ystats'<=0);                % For log transformations, if any ystats are
      if any(neg)                           % negative, replace negative values by the 
        ys = ystats(:);                     % smallest positive statistic values
        ysmin = min(ys(ys>0));
        for i = 1:3
          neg = (ystats(:,i)<=0);
          if (length(neg)>0)
            ystats(neg,i) = ysmin * ones(sum(neg),1);
          end;
        end;
      end;
      ys = log(ystats);

    case 2,
      neg = any(ystats'<=0);                % For log transformations, if any ystats are
      if any(neg)                           % negative, replace negative values by the 
        ys = ystats(:);                     % smallest positive statistic values
        ysmin = min(ys(ys>0));
        for i = 1:3
          neg = (ystats(:,i)<=0);
          if (length(neg)>0)
            ystats(neg,i) = ysmin * ones(sum(neg),1);
          end;
        end;
      end;
      ys = log10(ystats);

    case 3,
      ys = exp(ystats);

    case 4,
      ys = 10.^(ystats);
 
    otherwise,
      error('  PLOTVERT: invalid valaue for xtransform');
  end;
  
  ux = uniquef(xvals);
  lenux = length(ux);
  xint = min(ux(2:lenux)-ux(1:(lenux-1)));
  xrng = max(ux)-min(ux);
%   tiplen = min([0.4*xint, 0.04*xrng])/2;
  tiplen = 0.02*xrng/2;

  plot(xvals,ys(:,2),colorh);
  hold on;
  for i = 1:nx
    xv = xvals(i);
    ys1 = ys(i,1);
    ys3 = ys(i,3);
    plot([xv,xv],[ys1,ys3],colorv);
    if (yextent)
      plot([xv-tiplen,xv+tiplen],[ys1 ys1],colorv);
      plot([xv-tiplen,xv+tiplen],[ys3 ys3],colorv);
    end;
  end;
  hold off;
  putbnd([xvals;xvals;xvals],ys(:));

  return;
