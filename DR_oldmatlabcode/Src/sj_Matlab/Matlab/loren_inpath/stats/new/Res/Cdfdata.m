% CDFDATA: Produces a cumulative distribution function corresponding to a vector
%          of data.
%
%     Usage: [xs,fxs] = cdfdata(x,{doplot})
%
%         x =      data vector.
%         doplot = optional boolean flag indicating, if true, that a plot of the 
%                    step-function is to be produced [default = 0].
%         ----------------------------------------------------------------------
%         xs =     sorted data vector.
%         fxs =    corresponding vector of cumulative frequencies.
%

% RE Strauss, 11/16/01

function [xs,fxs] = cdfdata(x,doplot)
  if (nargin < 2) doplot = []; end;

  if (isempty(doplot))
    doplot = 0;
  end;

  xs = sort(x);
  fxs = linspace(0,1,length(x));

  if (doplot)
    hold on;
    plot([xs(1),xs(1)],[0 fxs(1)],'k');
    for i = 2:length(xs)
      plot([xs(i-1),xs(i)],[fxs(i-1),fxs(i-1)],'k');
      plot([xs(i),xs(i)],[fxs(i-1),fxs(i)],'k');
    end;
    hold off;
    putbnds(xs,fxs);
  end;

  return;
