function [b,s2,mu]=regr(y,x,noconstant);
%REGR    [b,s2,mu]=REGR(y,x,'noconstant').  Least squares regression of y on
%        x.  The form REGR(y,x) includes an intercept, while
%        REGR(y,x,'noconstant') does not.  On output, b holds regression
%        coefficients, s2 the residual mean square and MU the fitted values.

% GKS  June 92.  Last revision 25 Jan 94.

% Ensure data in columns
y = y(:);
[my ny] = size(y);
[mx nx] = size(x);
if mx ~= my,
   if nx == my,
      x = x';
   else
      disp('Regr:  Dimensions of y and x don''t match.');
      return;
   end;
end;

% Remove missing values
x = excise([y x]);
y = x(:,1);

% Include constant term
if nargin < 3,
   x(:,1) = ones(my,ny);
else
   x(:,1) = [];
end;

b=x\y;
mu=x*b;
s2=(y-mu)'*(y-mu)./(length(y)-length(b));
se=sqrt(diag(inv(x'*x))).*sqrt(s2);

disp('     Coef      S.E.');
disp([b se]);

disp('    Res mean square');
disp(s2);
