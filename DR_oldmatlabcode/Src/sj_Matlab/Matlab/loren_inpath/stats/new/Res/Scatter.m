% SCATTER: Simple unlabeled scatter plot of the first two columns of a matrix.
%
%     Usage: scatter(x,{y},{regr},{symbol})
%
%           x =         data vector (if y is provided) or matrix containing at 
%                         least two columns.
%           y =         optional corresponding vector.
%           regr =      optional flag indicating that a 
%                         regression is to be performed and plotted:
%                           0 = no regression [default];
%                           1 = predictive regression;
%                           2 = major-axis regression.
%           symbol =    optional symbol & color string for plot statement 
%                         [default = 'ok'].
%

% RE Strauss, 1/2/99
%   5/23/99 -  allow input matrix or two input vectors.
%   8/19/99 -  change plot colors for Matlab v5.
%   2/24/00 -  added optional regression.
%   3/19/00 -  fix plotting of a single point.
%   10/11/00 - convert input row vectors to col vectors.

function scatter(x,y,regr,symbol)
  if (nargin < 2) y = []; end;
  if (nargin < 3) regr = []; end;
  if (nargin < 4) symbol = []; end;

  if (isscalar(y))                      % If y not passed, shift arguments right
    symbol = regr;
    regr = y;
    y = [];
  end;

  [isvect,ncells,iscol] = isvector(x);  % Convert row vectors to col vectors
  if (isvect & ~iscol)
    x = x';
  end;
  [isvect,ncells,iscol] = isvector(y);
  if (isvect & ~iscol)
    y = y';
  end;

  if (isempty(y))
    if (isvector(x) & length(x)>2)
      y = x(:);
      x = [1:length(x)]';
    else  
      y = x(:,2);
      x = x(:,1);
    end;
  else
    if (any(size(x)~=size(y)))
      error('  SCATTER: input vectors not compatible')
    end;
    x = x(:,1);
    y = y(:,1);
  end;

  if (isempty(symbol))
    symbol = 'ok';
  end;
  if (isempty(regr))
    regr = 0;
  end;

  if (regr & size(x,1)<2)
    disp('  SCATTER warning: too few points for regression');
    regr = 0;
  end;

  figure;
  plot(x,y,symbol);
  putbnd(x,y);

  if (regr)
    xmin = min(x);
    xmax = max(x);

    switch (regr)
      case 1,
        [b,stats,pred] = linregr(x,y,[xmin xmax]);

      case 2,
        b = majaxis(x,y);
        pred = [xmin 1; xmax 1;]*b(1,:)';

      otherwise
        error('  SCATTER: invalid regression flag');
    end;

    hold on;
    plot([xmin xmax],pred,'k');
    hold off;
  end;

  return;
