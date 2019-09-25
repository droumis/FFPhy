function y = tricube(x)
%TRICUBE Tri-cube kernel.
%
%   TRICUBE(X) accepts real X and returns
%     for X >= 1 or X =< 1 : 0
%     for -1 < X < 1 : (1 - ABS(X).^3).^3
%
%Written by smk, 2009 August 19.
%

if ~isnumeric(x) || ~isreal(x)
  error('TRICUBE accepts only real inputs');
end

y = (abs(x) < 1).*(1 - abs(x).^3).^3;

