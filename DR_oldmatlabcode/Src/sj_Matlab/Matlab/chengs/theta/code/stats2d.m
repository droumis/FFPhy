% [mx,my,Vx,Vy,R,I]=Hist2Stat(z,x,y,[norm,m])
% 
% Calculate linear correlation coefficient from 2-variate (x,y) probability 
% density function (z).
% 
% input
%   z       matrix representing PDF at discrete values, z(x,y)=N*pdf(x,y)*dx*dy
%   x       vector with values of center of bin
%             (does not have to be equidistant)
%   y       vector with values of center of bin
%             (does not have to be equidistant)
%   norm    optional, 1 if z is properly normalized, default is no
%   m       optional, mean of x, y
%
% output
%   mx,my   means of x and y
%   Vx,Vy   variances of x and y
%   R       linear correlation coefficient
%   I       mutual information

function [mx,my,Vx,Vy,R,I]=Hist2Stat(z,x,y,norm,m)

%keyboard
tiny=1e-30;

nx=length(x);
ny=length(y);
nz=size(z);
if nz(1)~= nx; error('size(z,1) should equal size(x)'); end
if nz(2)~= ny; error('size(z,2) should equal size(y)'-1); end
if min(min(z)) < 0; error('PDF has to be non-negative'); end

% determine normalization constant of z
if nargin <= 3   | norm<=0
    norm=sum(sum(z));
end

if norm== 0
    mx=nan; my=nan; Vx=nan; Vy= nan; R= nan; I= nan;
else
% normalize z
z=z/norm;

% marginal distributions
px=zeros(nx,1); py=zeros(ny,1);
for (j=1:nx)
    for (k=1:ny)
	px(j)=px(j)+ z(j,k);
	py(k)=py(k)+ z(j,k);
    end
end


% determine means
if nargin <=4
% determine mean(x)
    mx=0; my=0;
    for (j=1:nx) mx=mx+ x(j)*px(j); end
    for (k=1:ny) my=my+ y(k)*py(k); end
else
    mx=m(1);
    my=m(2);
end

% determine summed squared errors and S_XY
Sx=0; Sy=0; Sxy=0;
for (j=1:nx)
    for (k=1:ny)
	Sx=Sx+ (x(j)-mx)^2*z(j,k);
	Sy=Sy+ (y(k)-my)^2*z(j,k);
	Sxy=Sxy+ (x(j)-mx)*(y(k)/2-my)*z(j,k);
    end
end

if(Sx==0 | Sy==0)
    if Sxy==0
	R=0;
    else 
	error('do not know what to do');
    end
else
    R=Sxy/sqrt(Sx*Sy);
end
Vx=Sx;
Vy=Sy;

% mutual information
I=0;
for (j=1:nx)
    for (k=1:ny)
        if z(j,k)>tiny
            I = I + z(j,k)*log2(z(j,k)/(px(j)*py(k)));
        end
    end
end

end
%keyboard

