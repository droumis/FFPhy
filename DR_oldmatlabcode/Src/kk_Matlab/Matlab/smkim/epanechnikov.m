function y = epanechnikov(x)
%EPANECHNIKOV Epanechnikov kernel.
%
%   EPANECHNIKOV(X) accepts real X and returns
%     for X >= 1 or X =< 1 : 0
%     for -1 < X < 1 : 0.75*(1 - ABS(X).^2);
%
%Written by smk, 2009 August 19.
%

if ~isnumeric(x) || ~isreal(x)
  error('TRICUBE accepts only real inputs');
end

y = 0.75.*(abs(x) < 1).*(1 - abs(x).^2);

