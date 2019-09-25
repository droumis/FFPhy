% CUMSTEP: Converts a discrete pdf into a cumulative relative cdf suitable 
%          for Kolmogorov-Smirnov analysis.
%
%     Usage:  cdf = cumstep(x,{f},{range})
%
%           x =     vector (length n) of integer X values.
%           f =     optional corresponding vector of frequencies.
%           range = 2-element vector containing optional range for abscissa 
%                   [default = observed range].
%           -----------------------------------------------------------------
%           cdf =   [(2n+1) x 2] matrix of point coordinates defining the 
%                     cumulative distribution.
%

% RE Strauss, 4/30/97
%   9/19/99 - updated handling of default input arguments.

function cdf = cumstep(x,f,range)
  if (nargin < 2) f = []; end;
  if (nargin < 3) range = []; end;

  if (isempty(f))
    f = ones(size(x));
  end;

  if (isempty(range))
    xmin = min(x);
    xmax = max(x);
  else
    if (length(range)==2)  
      xmin = range(1);
      xmax = range(2);
    else
      error('  Range vector is wrong size');
    end;
  end;

  if (min(size(x))>1 | min(size(f))>1)
    error('  Input must be vectors');
  end;

  n = length(x);
  if (length(f)~=n)
    error('  Input vectors must be same length');
  end;

  [x,i] = sort(x);            % Sort input vector, just in case
  f = f(i);

  f = f/sum(f);               % Convert freqs to cum relative freqs
  f = cumsum(f);

  lencdf = 2*n + 2;
  cdf = zeros(lencdf,2);
  cdf(1,:) = [xmin, 0];
  cdf(lencdf,:) = [xmax, 1];
  j = 1;
  
  for i = 1:n
    cdf(j+1,:) = [x(i),cdf(j,2)];
    cdf(j+2,:) = [x(i),f(i)];
    j = j+2;
  end;

  
  return;
