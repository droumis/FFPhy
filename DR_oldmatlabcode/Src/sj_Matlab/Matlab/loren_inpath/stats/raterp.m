function y=raterp(a,b,x)
% y=raterp(a,b,x) interpolates using coefficients
% from function ratcof.
% by H.B. Wilson, Oct 1988
a=fliplr(a(:)'); b=fliplr(b(:)');
y=polyval(a,x)./(1+x.*polyval(b,x));
