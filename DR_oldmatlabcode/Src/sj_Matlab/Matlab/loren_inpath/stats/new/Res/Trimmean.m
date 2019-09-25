% TRIMMEAN: Estimate the trimmed mean as a function of t, the proportion of 
%           the sorted data trimmed at each tail.  Deletes points in pairs, 
%           one at each end.  Plots the function in refence to the median.
%
%     Syntax: [m,t] = trimmean(X,{noplot},{maxt},{nval})
%
%           X -      data vector.
%           noplot - boolean flag to suppress plot [default = 0];
%           maxt -   maximum proportion to be trimmed at each end
%                      [default = 0.25].
%           nval -   number of values sampled for function [default = 51].
%           ----------------------------------------------------------------
%           m -      column vector (of length 'nval') containing trimmed 
%                      means as a function of t.
%           t -      column vector (of length 'nval') containing proportions 
%                      of data trimmed from each tail.
%

% RE Strauss, 2/10/97
%   9/7/99 - miscellaneous changes for Matlab v5.

function [m,t] = trimmean(X,noplot,maxt,nval)
  if (nargin < 2) noplot = []; end;
  if (nargin < 3) maxt = []; end;
  if (nargin < 4) nval = []; end;

  default_noplot = 0;
  default_maxt = 0.25;
  default_nval = 51;

  if (min(size(X))>1)
    error('  TRIMMEAN: Input data must be vector');
  end;
  N = length(X);

  if (isempty(noplot))
    noplot = default_noplot;
  end;
  if (isempty(maxt))
    maxt = default_maxt;
  end;
  if (isempty(nval))
    nval = default_nval;
  end;

  if (maxt > 1)
    maxt = maxt/100;
  end;
  maxt = min(maxt,0.5);

  X = sort(X);
  t = [0:(maxt/(nval-1)):maxt]';
  d = floor(N * t);             
  m = zeros(nval,1);

  m(1) = mean(X);
  for i = 2:nval
    if (d(i) == d(i-1))
      m(i) = m(i-1);
    else
      dd = d(i)+1;
      x = X(dd:(N-dd+1));
      if (~isempty(x))
        m(i) = mean(x);
      else
        m(i) = m(i-1);
      end;
    end;
  end;

  if (~noplot)
    med = median(X);
    clf;
    plot(t,m,'k');
    hold on;
    plot([0 maxt]',[med med]','k--');
    putbnd([0; t],[med; m]);
    hold off;
    putxlab('Trimming proportion (at each end)');
    putylab('Trimmed mean (solid) & median (dashed)');
  end;

  return;
