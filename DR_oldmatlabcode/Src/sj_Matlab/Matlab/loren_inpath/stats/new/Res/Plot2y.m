% PLOT2Y: Produces a 2D plot for a single abscissa (X) but two ordinates (Y), 
%         scaled independently using the left and right vertical axes.
%
%     Usage: plot2y(x,y1,y1specs,y2,y2specs,{makesquare})
%
%         x =               vector (length n) of abscissa values.
%         y1,y2 =           vectors (length n) of values for ordinate variables.
%         y1specs,y2specs = string variables with plot characteristics.
%         makesquare =      optional boolean flag indicating, if true, that plot 
%                             should be square rather than rectangular 
%                             [default = 0].
%

% RE Strauss, 6,26,00

function plot2y(x,y1,y1specs,y2,y2specs,makesquare)
  if (nargin < 6) makesquare = []; end;

  if (isempty(makesquare))
    makesquare = 0;
  end;

  if (~isvector(x) | ~isvector(y1) | ~isvector(y2))
    error('  PLOT2Y: input variables must be vectors.');
  end;

  n = length(x);
  if (length(y1)~=n | length(y2)~=n)
    error('  PLOT2Y: input vectors not of identical length.');
  end;

  x = x(:);                             % Should be column vectors
  y1 = y1(:);
  y2 = y2(:);

  y2plot = setrange(y2,min(y1),max(y1));  % Rescaled y2 values for plotting

  plot(x,y1,y1specs);                   % Plot y1
  if (makesquare)                       % Optional square plot
    axis('square');
  end;
  putbnd(x,y1);                         % Set plot bounds
  hold on;
  plot(x,y2plot,y2specs);               % Plot y2c

  ylim = get(gca,'YLim')'              % Values and positions of y1 ticks
  yrng = ylim(2)-ylim(1)
  ytick = get(gca,'YTick')'
  yprop = (ytick-ylim(1))./yrng

  y2tick = setrange(ytick,min(y2),max(y2))  % Values of y2 ticks

%  r = range(y2tick);                      % Number of decimal positions
%r
%  if (r>10)
%    dp = 0;
%  else
%    dp = 1;
%    while (r<1)
%      dp = dp+1;
%      r = r*10;
%    end;
%dp
%  end;

  y2tick = tostr(y2tick,3)
  for i = 1:length(y2tick)
    puttext(1.01,yprop(i),y2tick(i,:),[],12);
  end;

  hold off;

  return;
