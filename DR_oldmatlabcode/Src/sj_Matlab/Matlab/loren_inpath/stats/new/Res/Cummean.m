% CUMMEAN: Finds cumulative means and variances for a vector.
%
%     Usage: [cum_mean,cum_var] = cummean(x,plotflag)
%
%           x =        input vector.
%           plotflag = optional flag indicating, if TRUE (=1), that the 
%                       cumulative means and variances are to be plotted 
%                       [default = 0 = FALSE].
%           ------------------------------------------------------------
%           cum_mean = corresponding vector of cumulative means.
%           cum_var =  corresponding vector of cumulative variances.
%

% RE Strauss, 11/4/97
%   9/19/99 - updated handling of default input arguments.

function [cum_mean,cum_var] = cummean(x,plotflag)
  if (nargin < 2) plotflag = []; end;

  if (isempty(plotflag))
    plotflag = 0;
  end;

  lenx = length(x);
  i = [1:lenx]';
  cum_mean = cumsum(x)./i;

  cum_var = cum_mean;
  for j = 1:lenx
    cum_var(j) = var(x(1:j));
  end;

  if (plotflag)
    subplot(1,2,1);
    plot(i,cum_mean,'-w');
    putbnd(i,cum_mean);
    axis('square');
    putxlab('Cum vector length');
    putylab('Mean');
    subplot(1,2,2);
    plot(i,cum_var,'-w');
    putbnd(i,cum_var);
    axis('square');
    putxlab('Cum vector length');
    putylab('Variance');
  end;

  return;
