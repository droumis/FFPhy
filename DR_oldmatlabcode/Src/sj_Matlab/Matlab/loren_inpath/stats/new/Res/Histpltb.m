% HISTPLTB: Plots histogram bars
%
%     Usage: histpltb(n,centers,C)
%
%           n =       vector of counts per bin.
%           centers = corresponding centers of bins.
%           C =       bar color.
%           

% RE Strauss, 11/25/99

function histpltb(n,centers,C)
  half_intrvl = (centers(2)-centers(1))/2;

  hold on;
  for ib = 1:length(n)                    % Plot bars
    xlow = centers(ib)-half_intrvl;
    xhigh = centers(ib)+half_intrvl;
    ylow = 0;
    yhigh = n(ib);

    xb = [xlow xlow  xhigh xhigh xlow];
    yb = [ylow yhigh yhigh ylow  ylow];
    fill(xb,yb,C);

    if (C~='w')
      plot(xb(1:4),yb(1:4),'w');
    else
      plot(xb(1:4),yb(1:4),'k');
    end;
  end;

  return;

