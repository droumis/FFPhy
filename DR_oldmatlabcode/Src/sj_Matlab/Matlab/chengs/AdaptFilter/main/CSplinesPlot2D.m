function CSplinesPlot2D(x,y,z,t)

% cardinal splines interpolation
%  xi,    control points, must be sorted, can be 
%            column vector (1-d) or two column vectors (2-d)
%  zi     measured response
%  x      evaluate fct here
%  t      tension parameter for cardinal splines, default= 0.5

uN=11;

if nargin < 4; t=0.5; end
u=[0:1/(uN-1):1]';
Ones=ones(uN,1);
M=[ -t, 2-t, t-2, t; 2*t, t-3, 3-2*t, -t; -t, 0, t, 0;0, 1, 0, 0];
xN=length(x);
yN=length(y);

zi=zeros(uN*[xN-3,yN-3]);
for j=2:xN-2
    for k=2:yN-2
        K=[z(j-1:j+2,k-1:k+2)];
        zi((j-2)*uN+1:(j-1)*uN,(k-2)*uN+1:(k-1)*uN)= ...
            ([u.^3 u.^2 u Ones]*M*K)*M'*[u.^3 u.^2 u Ones]';
    end
end

%mesh(zi);
imagesc((zi')
colorbar
aixs xy
