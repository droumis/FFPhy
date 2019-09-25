function y=modelrsm(p,x)
% To calc general Quadratic Form y=b0 + x*b +x*B*x'
% Includes all main effects , 2 cross products x(i)*x(j) also x(i)*x(i) (quads)
% used in Response Surface Methodology
% If no Quadratic terms {length(p) is < (w^2+w)/2}  then diag(B) is zeroed
% p=vector(b0,b,B)===> need to reshape p from long vector
% B was further packed since it is symmetric and zeros removed
% p=[b0 b1 b2 b3 B11 B12 B13 B22 B23 B33] for 3 x's (w=3)
% ALL CROSS PARAM'S E.G. B12 (NOT B11) ==> USE 2*B12 IN y=b0 +b1.x +b12.x1.x2
% ie b12=2*B12

% A Jutan May 98

[l,w]=size(x);% order of B =w*w
lp=length(p);
b0=p(1);
b=p(2:w+1);
pp=p(2+w:lp); % select only Bij terms in p
% Check for zero diag in B i.e. no Square terms ***PURE FACTORIAL DATA****
if length(pp) ~=(w^2 +w)/2
   % diags on B=0 reconstruct one order lower B then add diag of zeros
   a1=symunpck(pp,w-1);
   a2=tril(a1); % take lower triangle
   a3=[zeros(1,w-1);a2]; %add top row of zeros
   a4=[a3 , zeros(w,1)]; %add RHS row of zeros
   B=a4+a4';% form symmetric B with zeros on diag
   else
   B=symunpck(pp,w); % recreates symmetric matrix of B from vector form in pp
end
for i=1:l
y(i)=b0 + x(i,:)*b + x(i,:)*B*x(i,:)';
end
y=y';	% make sure col vector returned to leasqr
return
