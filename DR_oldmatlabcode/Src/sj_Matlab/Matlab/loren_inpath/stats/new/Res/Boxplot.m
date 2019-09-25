% BOXPLOT: Boxplot function
%
%     Usage: boxplot(X,Y,{nparm},{connect},{width})
%
%           X =      abscissa values at which individual boxplots will be plotted.
%           Y =      vector of values for which boxplot is to be plotted.
%           nparm =  optional boolean flag indicating summary statistics to be plotted:
%                      0 = mean, stdev, range [default]
%                      1 = quartiles
%           connect = optional boolean flag indicating that centers of continguous
%                       boxplots are to be connected by lines [default = 0].
%           width =  optional boxplot width, in units of X [default = 0.5].
%

% RE Strauss, 5/6/98
%   9/7/99 - changed plot colors for Matlab v5.

function boxplot(X,Y,nparm,connect,width)
  if (nargin < 3) nparm = []; end;
  if (nargin < 4) connect = []; end;
  if (nargin < 5) width = []; end;

  if (isempty(nparm))
    nparm = 0;
  end;
  if (isempty(connect))
    connect = 0;
  end;
  if (isempty(width))
    width = 0.5;
  end;

  [r,c] = size(Y);
  if (min([r,c]) > 1)
    error('BOXPLOT: input data must be vector rather than matrix');
  end;
  if (c > 1)                     % Transpose input vector if row matrix
    Y = Y';
  end;

  xval = uniquef(X,1);           % Find unique X values, in sequence
  nx = length(xval);
  concrds = zeros(nx-1,4);

  hold on;
  for i = 1:nx
    [lc,rc] = boxplotb(Y(X==xval(i)),xval(i),width,nparm);
    concrds(i,:) = [lc rc];
  end;

  if (connect)
    for i = 1:nx-1
      plot([concrds(i,3);concrds(i+1,1)],[concrds(i,4);concrds(i+1,2)],'k');
    end;
  end;
  hold off;

  [xm,i] = min(X);
  X(i) = X(i) - width;
  [xm,i] = max(X);
  X(i) = X(i) + width;

  putbnd(X,Y);
  puttick(xval);
  box on;

  return;

