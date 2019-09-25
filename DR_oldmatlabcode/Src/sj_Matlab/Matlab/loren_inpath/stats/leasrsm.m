function  [b0,b,B]=leasrsm(x,type)
%function  [b0,b,B]=leasrsm(x{,type})
% x=data file [x1 x2 ..xn y]
% Uses Curvefit FROM OPT PACKAGE
% type= 'q' for full BOX composite model with quad terms,else Factorial
% model with only cross terms is fit {}=optional
% This program sets up the data for a call to modelrsm.m
% for Least Squares calc of parameters for General Quadratic
% Model ===> used in Response Surfaces	y=b0 + x*b +x*B*x'
% Uses Subroutines: CURVEFIT ,MODELRSM,SYMUNPCK
%	    A.Jutan. UWO '98
[l,w]=size(x);
y=x(:,w);
xx=x;xx(:,w)=[];  % remove last col store x's in xx
w=w-1; % Number of x's
% INITIAL GUESS'S for params in b0 + x*b +x*B*x'
b0=1.0;
b=2.0*ones(w,1);
B=ones(w,w); % guess ones matrix no zeros
%**********************************************************
% FOR PURE FACTORIAL DATA I.E. NO BOX COMPOSITE (STAR) PTS
%**********************************************************
% SET DIAGONAL OF GUESSED B TO ZEROS
if( nargin == 1)
 B=B -diag(diag(B));
end;
%***********************************************************
%
% PACK GUESS INTO ONE LONG VECTOR
%
B=tril(B); % take lower half only since B symmetric
p1=B(:); % form long vector--includes zeros
i=find(p1==0);% find the zero positions
p1(i)=[]; % eliminate the zero's ==> reduces p vector
% pack into single column vector: length (n^2 +n)/2
% pin=[b0 b1 b2 b3 B11 B12 B13 B22 B23 B33]' for 3 x's (n=3)
% pin=[b0 b1 b2 b3  B12 B13  B23 ] for 3 x s no SQUARE terms
pin=[b0;b(:);p1(:)]; % set up as long vector incl b0 for init param guess
% CHECK FOR SUFFICIENT DATA TO EST PARAMS
lpin=length(pin);
if (l<lpin),error('NOT ENOUGH DATA TO ESTIMATE ALL PARAMETERS'),end;
%
% default iter=20;	% max iterations
%dp=fractional changes for numerical derivatives   default=+.001*ones(pin);
%[f,p,var,iter,cor,std]=leasqr1(t,y,pin,'model'{,tol,iter,wt,dp,'dfdp'});
%[f,p,var,iter,cor,std]=leasqr1(xx,y,pin,'modelrsm');
[p,options,error,jac]=curvefit('modelrsm',pin,xx,y);
ymodel=y+error;
[std,varresid,r2,cor,vcv,varinf]=regdata(p,ymodel,y,jac);
% reshape parameter matrices
disp(' parameter matrices and 95% Confidence Limits +/- CL95')
b0=p(1)
stdb0=std(1)	 ;CL95_b0=2*stdb0
b=p(2:w+1)'
stdb=std(2:w+1)' ;CL95_b=2*stdb
lp=length(p);
pp=p(2+w:lp); % select only Bij terms in p
stdpp=std(2+w:lp); % select only Bij terms in std
lpp=length(pp);
% check for diag B=0 ie NO SQUARE TERMS
if lpp~=(w^2+w)/2 ,w1=w-1;else,w1=w;,end; % Unpack B without diag(quad) PARAMS
B=symunpck(pp,w1);% recreates symmetric matrix of B from vector form in pp
B=tril(B);    % take lower half
stdB=tril( symunpck(stdpp,w1));
% Add Back Zero Diagonals for no quadratic terms
      if lpp~=(w^2+w)/2
      B=[zeros(1,w);B,zeros(w1,1)]; % advanced fiddle see p.2-39 MATLAB
      stdB=[zeros(1,w);stdB,zeros(w1,1)];
      end
B  % show B matrix and 95% CL
CL95_B=2* stdB
disp('p=[b0 b1 b2 b3 B11 B12 B13 B22 B23 B33] for 3 x s (n=3)')
disp('p=[b0 b1 b2 b3  0  B12 B13  0  B23  0 ] for 3 x s no SQUARE terms')
disp('b12=2*B12 in y= b0 + b1.x1 +b2.x2 +b12.x1.x2 + b11.x1^2+b22.x2^2')
disp(' non cross terms b11 b22 equal to B11 B22,etc')
