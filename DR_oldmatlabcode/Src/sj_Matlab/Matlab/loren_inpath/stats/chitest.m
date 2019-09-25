function Q=chitest(y,theta,sig)
%Q = chitest(y,theta,sig)
%Accepts a univariate time series y, a theta matrix th for the proposed model
%(see ARX, ARMAX) and a significance level sig.  The chi^2 test is performed on
%the first 20 autocorrelations (excluding the zero lag).  Note that the 
%significance level is for the inverse chi^2 cumulative distribution function
%and therefor find the value that exceeds sig*100% of the samples from a chi^2
%distribution with the degree of freedom of the autocorrelation


format compact
%Find the difference between data and model
e = resid(y,theta);
n = length(e);

%Find the autocorrelation of the residual
%for lags 0 to 20

ac = coxy(e,e,'cor',21);

%Disgard autocorrelation at lag 0
%as it is always 1
ac=ac(2:21);

%Degrees of Freedom for the inverse Chi^2

disp('Degrees of Freedom for the Test:')
DOF = 20 - theta(1,4) - theta(1,5)
disp('Chitest exceeds 95% of the samples from a chi^2') 
disp('distribution with DOF degrees of freedom')
%Find the value of the Chi^2 cdf form zero
%to sig at degrees of freedom
chitest = chi2inv(sig, DOF)
disp('You would observe values of Q greater than')
disp('Chitest only (1-sig)*100% of the time by chance')
Q = n.*sum(ac.^2)

if Q <= chitest
disp('Passed at significance level')
else
disp('Failed at significance level')
end