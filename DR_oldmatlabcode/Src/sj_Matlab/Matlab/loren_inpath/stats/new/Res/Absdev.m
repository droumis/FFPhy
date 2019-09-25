% ABSDEV: Converts a data matrix to absolute deviations from the center, by column.
%
%     Usage: dev = absdev(X,{usemedian})
%
%         X =         [n x p] data matrix.
%         usemedian = optional boolean flag: 
%                       0 = abs deviations from mean [default].
%                       1 = abs deviations from median.
%         ------------------------------------------------------
%         dev =       corresponding [n x p] matrix of deviations.
%

function dev = absdev(X,usemedian)
  if (nargin < 2)
    usemedian = [];
  end;
  if (isempty(usemedian))
    usemedian = 0;
  end;

  [n,p] = size(X);

  if (usemedian)
    dev = abs(X - ones(n,1)*median(X));
  else
    dev = abs(X - ones(n,1)*mean(X));
  end;

  return;
