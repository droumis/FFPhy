% MINKPDF:  PDF of a Minkowski distribution.
%
%     Usage: y = minkpdf(x,k)
%
%           x = vector of values for which pdf values are desired.
%           k = Minkowski coefficient.
%           ------------------------------------------------------
%           y = coordinates of pdf.
%

% RE Strauss & EG Dyreson, 5/25/01

function y = minkpdf(x,k)
  xismat = 0;
  if (ismatrix(x))
    [r,c] = size(x);
    x = x(:);
    xismat = 1;
  end;

  n = 2000;
  xxmax = range(x)+2;
  if (k>=1)
    xxmax = max([xxmax,15]);
  else
    xxmax = max([xxmax, 100]);
  end;
  xx = linspace(1,xxmax,n)';              % Create pdf
  xx = xx-mean(xx);

  yy = exp(-(abs(xx).^k)/2);   

  dx = xx(2)-xx(1);
  A = sum(yy*dx);
  yy = yy ./ A;

  y = zeros(size(x));
  for i = 1:length(x)
    d = xx - x(i);
    [m,indx] = max(d(d<0));

    x1 = yy(indx);                        % Linearly interpolate
    x2 = yy(indx+1);
    y1 = d(indx);
    y2 = d(indx+1);

    y(i) = x2 - y2*(x2-x1)/(y2-y1);
  end;

  if (xismat)
    y = reshape(y,r,c);
  end;

  return;
