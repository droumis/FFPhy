% OrthogPolyRegrP: Produces predicted values for orthogonal polynomial regression.
%
%     Usage: y = orthogpolyregrp(x,b,alpha,beta)
%
%         x =     column vector of values of the independent variable at which
%                   the function is to be evaluated..
%         b =     column vector of regression coefficients.
%         alpha = alpha coefficients used to calculate polynomials.
%         beta =  beta coefficients used to calculate polynomials.
%         --------------------------------------------------------------------
%         y =     predicted values corresponding to x.
%

% RE Strauss, 2/14/02
%   3/7/03 - pass particular values at which regression function is to be evaluated.

function y = orthogpolyregrp(x,b,alpha,beta)
  x = x(:);
  b = b(:);
  k = length(x);
  degree = length(b)-1;
  P = zeros(k,degree+1);

  P(:,1) = ones(k,1);                 % Term 0
  P(:,2) = (x-alpha(2)).*P(:,1);      % Term 1

  for term = 2:degree                 % Remaining terms
    j = term+1;
    P(:,j) = (x-alpha(j)).*P(:,j-1) - beta(j-1)*P(:,j-2);
  end;
  
  y = P*b;                            % Predicted values at x's

  return;
  
