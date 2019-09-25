% OrthogPolyRegrF: Objective function for orthogonal polynomial regression
%
%     Usage: [b,alpha,beta,P] = orthogpolyregrf(x,y,degree)
%
%         x =       [N x 1] vector of unique independent-variable values.
%         y =       [N x 1] vector of sums of dependent-variable values at 
%                     each level of x.
%         degree =  degree of polynomial (=k).
%         ------------------------------------------------------------------------------------
%         b =       [k+1 x 1] vector of regression coefficients.
%         alpha =   [k+1 x 1] vector of alpha coefficients used to calculate polynomials.
%         beta =    [k+1 x 1] vector of beta coefficients used to calculate polynomials.
%         P =       [N x k+1] matrix of observations (rows) for each orthogonal 
%                     polynomial term (cols).
%

% Dutka, AF and FJ Ewens. 1971. A method for improving the accuracy of polynomial
%   regression analysis. J Quality Technology 3:149-155.

% RE Strauss, 2/14/02
%   3/7/03 - revised sequence of output arguments;
%            provide divide-by-zero protection.

function [b,alpha,beta,P] = orthogpolyregrf(x,y,degree)
  k = length(x);
  
  P = zeros(k,degree+1);
  b = zeros(degree+1,1);
  ss = zeros(degree+1,1);
  alpha = zeros(degree+1,1);
  beta = zeros(degree+1,1);
  
  term = 0;                           % Constant term
  j = term+1;
  P(:,j) = ones(k,1);
  b(j) = sum(y.*P(:,j))./sum(P(:,j).^2);
  
  term = 1;                           % Linear term
  j = term+1;
  alpha(j) = sum(x.*P(:,j-1).^2)./sum(P(:,j-1).^2);   
  P(:,j) = (x-alpha(j)).*P(:,j-1);
  b(j) = sum(y.*P(:,j))./sum(P(:,j).^2);
  beta(j) = sum(x.*P(:,j-1).*P(:,j)) ./ sum(P(:,j-1).^2);

  for term = 2:degree                 % Remaining terms
    j = term+1;
    sumP = sum(P(:,j-1).^2);
    if (sumP>eps)
      alpha(j) = sum(x.*P(:,j-1).^2) ./ sumP;
    else
      alpha(j) = NaN;
    end;
    P(:,j) = (x-alpha(j)).*P(:,j-1) - beta(j-1)*P(:,j-2);
    if (sumP > eps)
      beta(j) = sum(x.*P(:,j-1).*P(:,j)) ./ sumP;
    else
      beta(j) = NaN;
    end;
    sumP = sum(P(:,j).^2);
    if (sumP > eps)
      b(j) = sum(y.*P(:,j))./sumP;     
    else
      b(j) = NaN;
    end;
  end;

  return;