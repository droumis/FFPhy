function area = gquad6(fun,xlow,xhigh,mparts)
%GQUAD6	Six point Gauss quadrature.
%	AREA = GQUAD6(FUN,XLOW,XHIGH,MPARTS) determines the area under
%	an externally defined function FUN(X) between limits XLOW and
%	XHIGH.  The whole interval is divided into MPARTS subintervals
%	and the integration over each subinterval is done with a six
%	point Gauss formula.  See also GQUAD, GRULE.

%	Howard B. Wilson, U. of Alabama, Spring 1990
%	Revised GKS 5 June 92

%  The weight factors are
wf = [ 1.71324492379170d-01;   3.60761573048139d-01;...
        4.67913934572691d-01]; wf=[wf;wf([3 2 1])];

%  The base points are
bp = [-9.32469514203152d-01;  -6.61209386466265d-01;...
      -2.38619186083197d-01]; bp=[bp;-bp([3 2 1])];

d = (xhigh - xlow)/mparts;  d2 = d/2;  nquad = length(bp);
x = (d2*bp)*ones(1,mparts) + (d*ones(nquad,1))*(1:mparts);
x = x(:) + (xlow-d2); fv=feval(fun,x); wv = wf*ones(1,mparts);

area=d2*(wv(:)'*fv(:));
