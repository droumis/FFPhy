% MAJAXIS:  Performs major-axis regression on two vectors x & y, 
%           returning the slopes and intercepts of the major and minor axes.
%
%     Usage: b = majaxis(x,y)
%
%           X,Y = vectors of identical size.
%           ---------------------------------------------------------------
%           b =   [2 x 2] row vector of regression coefficients; the first 
%                   column gives the major and minor slopes, the second the 
%                   major and minor intercepts.
%

% RE Strauss, 11/21/97
%   11/9/01 - produce NaN's if cov(x,y)==0.
%   11/20/01 - changed abs(cov) to abs(covxy) in if statement.

function b = majaxis(x,y)
  lenx = length(x);
  leny = length(y);
  if (lenx ~= leny)
    error('  MAJAXIS: Input vectors must be of same size');
  end;

  meanx = mean(x);
  meany = mean(y);
  varx = var(x);
  vary = var(y);
  c = cov(x,y);
  covxy = c(1,2);

  slope = zeros(2,1);
  intcpt = zeros(2,1);

  if (abs(covxy)>eps)
    slope(1) = (vary - varx + sqrt((vary-varx)^2+4*covxy^2)) / (2*covxy);
    intcpt(1) = meany - slope(1)*meanx;
  else
    slope(1) = NaN;
    intcpt(1) = NaN;
  end;

  slope(2) = -1/slope(1);
  intcpt(2) = meany - slope(2)*meanx;

  b = [slope, intcpt];
  return;

