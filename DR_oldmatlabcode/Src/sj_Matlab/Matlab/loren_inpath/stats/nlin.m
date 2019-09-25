function [ypred,beta] = nlin(x,y,model,beta0)
%function [ypred,beta] = nlin(x,y,'model',beta0)
% NonLinear fitting routine ~ uses Marquardt and nlinfit from Stats Toolbox
%   nlin(x,y,'model',beta0) finds the coefficients of the nonlinear 
%   function described in MODEL. MODEL is a user supplied function having 
%   the form y = f(beta,X). That is MODEL returns the predicted values of y
%   given initial parameter estimates, beta, and the independent Matrix 
%   variable, X.   
%   [ypred,beta] = nlin(x,y,'model',beta0) returns the estimates ypred
%   and the fitted coefficients beta. 
%   95% Confidence Limits,Covariance and Correlation for Parameters are calc
%	A.Jutan 5-30-96

[beta,resid,J] = nlinfit(x,y,model,beta0);
[ci,varb,corrb,varinf] = nlparci(beta,resid,J);
[ypred, delta] = nlpredci(model,x,beta,resid,J);
disp('Parameter Estimates beta and 95 Percent Confidence Limits')
disp([beta,ci])
disp('Correlation matrix')
disp(corrb)
disp('Covariance matrix')
disp(varb)
disp('Variance Inflation: if>10 ==> multicollinearity in X')
disp(varinf)
disp('')

%Calculate R^2 (Ref Draper & Smith p.46)
      r=corrcoef(y,ypred);
      r2=r(1,2).^2;
      disp('Variance of Residuals  ' )
      disp(  var(resid) )
      disp( 'Correlation Coefficient R^2')
      disp(r2)

%disp('     ydata    ypred')
%disp([y  ypred])

% plot data plus predictions plus 95% errorbars
errorbar(ypred,delta);
hold on
plot(y,'go'); % add data point
hold off
title (' ydata o vs ypred - with 95pc errorbars')
xlabel('Observation number')
ylabel('Data o vs Prediction- and Errorbars')

