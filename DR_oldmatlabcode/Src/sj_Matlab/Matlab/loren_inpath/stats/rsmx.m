function xo=rsmx(x)
%function xo=rsmx(x)
%constructs full matrix xo of main effects x, cross terms x1.x2, etc
% and quadratics x1^2, etc. given the effects matrix x
% input : x  main effects matrix  (n experiments * m main effects)
% output: xo  main effects, 2*cross products, quadratics arranged as:
% x1 x2... B(:), where B(:) are the lower triangular cols of B :
% B11
% B12 B22
% B13 B23   B33   so xo will have the form:
%
% xo=[ x1 x2 x3 x1^2 x12 x13 x2^2 x23 x3^2] for 1st row, etc.


[n,m]=size(x);
%xo=zeros(n,(m.^2 +m)/2+m); % set up size of xo
for i=1:n
b1=x(i,:)'*x(i,:); % calc outer product
b2=tril(b1);% lower triangle
B=b2(:)';%form row vector
xo(i,:)=[x(i,:),B];
end
% remove cols of xo that are ALL zero
i=find((xo)==0); % find col #'s of cols which = 0
xo(:,i)=[]; % delete these cols
