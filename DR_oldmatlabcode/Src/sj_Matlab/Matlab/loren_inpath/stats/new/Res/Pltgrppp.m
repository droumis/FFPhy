% PLTGRPPP: Plot points or symbols for function plotgrps()
%

% RE Strauss, 5/31/99
%   9/2/99 - changed plot colors for Matlab v5.

function pltgrppp(x,y,g,symbol)
  plot_symbols = 0;
  if (~isempty(symbol))
    plot_symbols = 1;
    if (~ischar(symbol(1)))           % Convert numeric symbols to text
      symbol = tostr(symbol);
    end;
    [rs,cs] = size(symbol);           % Ensure that symbols are in column form
    if (rs==1 & cs>1)
      symbol = symbol';
    end;
  end;
  if (plot_symbols)                   % Find displacements for symbols
    dx = 0.01  * range(x);
    dy = 0;
  end;

  if (~plot_symbols)
    plot(x,y,'ko');
  else
    for i = 1:length(x)
      text(x-dx,y+dy,symbol(g,:));
    end;
  end;

  return;


