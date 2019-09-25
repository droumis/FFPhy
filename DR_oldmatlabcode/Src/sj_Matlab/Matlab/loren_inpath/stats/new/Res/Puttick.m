% PUTTICK:  Modifies the tick labels of an existing, current plot.
%
%     Usage: puttick(xtick,{ytick},{fontsize})  OR
%            puttick(xdp,{ydp},{fontsize})
%
%           xtick =     Vector of numeric values, or matrix of string values,
%                         for which tick marks and labels are desired on the 
%                         x-axis.  Ignore if null.
%                       If xtick = 'off', tick marks and labels are suppressed.
%                       If xtick is a scalar (xdp), the current tick labels are 
%                         adjusted to that number of decimal positions.
%           ytick =     Same, for the y-axis.
%           fontsize =  optional font size for tick labels.
%

% RE Strauss, 9/18/98
%   9/21/98 - allow suppression of tick marks and labels.
%   11/7/98 - allow change of font size.
%   9/2/99 -  change 'set' parameter values for Matlab v5.
%   1/10/00 - allow for string labels and setting of number of decimal positions.
%   5/5/01 -  fix bug with xtick & ytick string comparisons.
%   6/5/01 -  compare xtick & ytick == 'off' only if string is of length 3.
%   7/20/01 - fix bug with printing all string labels.

function puttick(xtick,ytick,fontsize)
  if (nargin < 2) ytick = []; end;
  if (nargin < 3) fontsize = []; end;

  xsupr = 0;
  ysupr = 0;

  if (ischar(xtick) & isvector(xtick) & length(xtick)==3)
      if (xtick(:)' == 'off')
        xtick = [];
        xsupr = 1;
      end;
  end;
  if (ischar(ytick) & isvector(ytick) & length(ytick)==3)
    if (ytick(:)' == 'off')
      ytick = [];
      ysupr = 1;
    end;
  end;

  if (isempty(xtick))
    if (xsupr)
      set(gca,'XTick',[]);
      set(gca,'XTickLabel',[]);
    end;
  elseif (isscalar(xtick))
    xl = str2num(get(gca,'XTickLabel'));
    xl = tostr(xl,xtick);
    set(gca,'XTickLabel',xl);
  else
    if (ischar(xtick))
      set(gca,'XTick',1:size(xtick,1));
    else
      set(gca,'XTick',xtick);
    end;
    set(gca,'XTickLabel',xtick);
  end;

  if (isempty(ytick))
    if (ysupr)
      set(gca,'YTick',[]);
      set(gca,'YTickLabel',[]);
    end;
  elseif (isscalar(ytick))
    yl = str2num(get(gca,'YTickLabel'));
    yl = tostr(yl,ytick);
    set(gca,'YTickLabel',yl);
  else
    if (ischar(ytick))
      set(gca,'Ytick',1:size(ytick,1));
    else
      set(gca,'YTick',ytick);
    end;
    set(gca,'YTickLabel',ytick);
  end;

  if (~isempty(fontsize))
    set(gca,'FontSize',fontsize);
  end;

  return;



