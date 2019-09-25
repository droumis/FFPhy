function zi = CSplines2D(x,y,z,xi,yi,t)

% cardinal splines interpolation
%  xi,    control points, must be sorted, can be 
%            column vector (1-d) or two column vectors (2-d)
%  zi     measured response
%  x      evaluate fct here
%  t      tension parameter for cardinal splines, default= 0.5

if nargin < 6; t=0.5; end
M=[ -t, 2-t, t-2, t; 2*t, t-3, 3-2*t, -t; -t, 0, t, 0;0, 1, 0, 0];
xN=length(x);
yN=length(y);

j=max(find(x<=xi));
k=max(find(y<=yi));

u=(xi-x(j))/(x(j+1)-x(j));
v=(yi-y(k))/(y(k+1)-y(k));

K=[z(j-1:j+2,k-1:k+2)];
zi=([u^3 u^2 u 1]*M*K)*M'*[v^3 v^2 v 1]';


