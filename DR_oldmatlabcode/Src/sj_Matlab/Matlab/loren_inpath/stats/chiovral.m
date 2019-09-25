function [S,chitest] = chiovral(z,th)
% [S,chitest]=chiovral(z,th) is the syntax for using this function.
% This function requires the data matrix, z=[y u], and the complete bj 
% model, th, to determine if the model fit is adequate. It uses the Chi-squared
% test with a 95% significance level on cross correlations between 
% the residuals and the input signal. 
% The values which are compared are S and chitest.
% S is the value of the sum of the squared cross correlations between
% the residuals and the input series, multiplied by the number of residuals.  
% chitest is the value of inverse of the chi squared cumulative
% distribution function.  
% 

% Splits the data matrix into output and input series vectors
y=z(:,1);
u=z(:,2);

% Determine the RESIDUAL white noise series
e=resid(z,th);
n=length(e);


% Cross Correlation between input and residuals
r=coxy(u,e,'cor',20);


% Picks out the degree of AR and MA from the model (to be used in
% calculating the degree of freedom for the chi2inv).
na=th(1,4);
nb=th(1,5);

% Calculation of the statistics to be compared using the 
% Chi Squared test. 
degfree=20-na-nb;
S=n.*sum(r.^2);
chitest=chi2inv(.95,degfree);

% Comparison section.
if S<chitest 
 disp(' ')
 disp('S < chitest :)')
 disp('PASSED Test')
else
 disp(' ')
 disp('S > chitest')
 disp('Failed test :(')
end

