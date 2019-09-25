function x=excise(x,miss_code);
%EXCISE	Remove rows with missing values.
%	EXCISE(X) removes rows of X containing any NaNs.
%	EXCISE(X,CODE) uses the value CODE instead of NaN.

%	GKS  3 May 92, 4 June 93.

if nargin==1,
   miss=isnan(x);
else
   miss=(x==miss_code);
end;

[m n]=size(x);
if n==1,
   x(miss) = [];
else
   x(any(miss'),:) = [];
end;
