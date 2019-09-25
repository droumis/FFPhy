function y = bisquare(x)
%BISQUARE Bi-square kernel.
%
%   BISQUARE(X) accepts real X and returns
%     for X >= 1 or X =< 1 : 0
%     for -1 < X < 1 : (1 - X.^2).^2
%
%Written by smk, 2009 August 19.
%

if ~isnumeric(x) || ~isreal(x)
  error('BISQUARE accepts only real inputs');
end

y = (abs(x) < 1).*(1 - abs(x).^2).^2;

