function y=maxrsm1(x,p1,p2,p3,p4)
%function y=maxrsm1(x{,p1,p2,p3,p4})
% Subroutine model for SMAXRSM1.M /OPTRSM  for optimising RSM models
% y is of the form y=b0 + x.b +x'.B.x includes quad or no quad terms cases.
% IN: x = guess for optimum factor level to optimise y
%     global :
%     b0_,b_,B_ =  Parameter matrices from RSM Model out of LEASRSM.M
%     opt_ = 'max' or 'min' depends	on problem
% OUT: y = calc from above model
% Note input x is subject to a contraint  as determined by lohi_ matrix
% as input to CONSTR.M

 if nargin ==5 , b0_=p1;b_=p2;B_=p3;opt_=p4;end % set variables if not GLOBAL
 xc=x'; % make row
y=b0_ + xc*b_' +xc*B_*xc'; % calc y for rsm model
if opt_=='max',y=-y;end % switch sign for maximizing default is min
return
