function s=ma(y,wind)
%MA	MA(Y,WIND) is the moving average smoother of Y with
%	window width WIND.  Default for WIND is 5.

y = y(:);
if nargin < 2, wind = 5; end;

[m n] = size(y);
if m < wind, disp('Window wider than sample'); return; end;

s = cumsum( [sum(y(1:wind)); y(wind+1:m)-y(1:m-wind)] )./wind;
