function zi = CSplines1D(cpx,z,x,t)

% cardinal splines interpolation
%  cpx,    control points, must be sorted, can be 
%            column vector (1-d) or two column vectors (2-d)
%  zi     measured response
%  x      evaluate fct here
%  t      tension parameter for cardinal splines, default= 0.5

if nargin < 4; t=0.5; end
M=[ -t, 2-t, t-2, t; 2*t, t-3, 3-2*t, -t; -t, 0, t, 0;0, 1, 0, 0];

[ncpx tmp]=size(cpx);
if ncpx==1 && tmp > 1
    cpx= cpx';
    ncpx= tmp;
end;

xN= length(x);
zi=zeros(xN,1);
for k=1:xN
    j=max(find(cpx<=x(k)));
    K=z(j-1:j+2);
    u=(x(k)-cpx(j))/(cpx(j+1)-cpx(j));

    zi(k)=[u^3 u^2 u 1]*M*K;
end


