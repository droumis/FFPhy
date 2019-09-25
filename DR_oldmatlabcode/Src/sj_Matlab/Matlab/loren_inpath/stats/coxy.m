function r=coxy(x,y,covr,n)
%function r=coxy(x,y,'covr',n)
%   Calculates the covariance or correlation of x to y (+ve lags only)
%   For -ve lags use coxy(y,x). x and y are automatically detrended
%   x, y - series
%   'covr' = 'cov' covariance(default)   or 'cor' correlation
%   n	- number of lags (default = 20)

y=y(:);x=x(:);% col vectors
z=[y x];
z=dtrend(z);
if nargin<=3 , n=20;end
R=covf(z,n);
if nargin==2 , covr='cov'; end
    if covr=='cov'
	r=R(2,:)';
    elseif covr=='cor'
	sccf=sqrt(R(4,1).*R(1,1)); % sqrt(var(y)*var(x))
	r=R(2,:)'./sccf;
    end 	
