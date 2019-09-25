% OrthogPolyRegr: Orthogonal polynomial regression.
%
%     Usage: [b,R2a,F,pF,ss,df,ms,alpha,beta] = orthogpolyregr(x,y,degree,{doplot})
%
%       x,y =     vectors of independent and dependent variables.
%       degree =  degree of polynomial (=k).
%       doplot =  optional boolean flag indicating, if true, that a plot of data
%                   and fitted polynomial is to be produced [default=0].
%       -------------------------------------------------------------------------------
%       b =       [k+1 x 1]  vector of regression coefficients.
%       R2a =     matching vector of cumulative adjusted coefficients of determination.
%       F =       [k x 1] vector of incremental, independent F statistics.
%       pF =      corresponding vector of corresponding probabilities.
%       ss =      vector of sums-of-squares:
%                   ssr:  [k x 1] vector of incremental, independent regression sums-of-squares;
%                   sse:  error sum of squares for full model;
%                   ssto: total sum of squares.
%       df =      corresponding vector of degrees of freedom:
%                   dfr:  [k x 1] vector of regression df;
%                   dfe:  error df;
%                   dfto: total df.
%       ms =      corresponding vector of mean squares:
%                   msr:  [k x 1] vector of regression mean squares;
%                   mse:  error mean square.
%       alpha =   [k+1 x 1] vector of alpha coefficients used to calculate polynomials.
%       beta =    [k+1 x 1] vector of beta coefficients used to calculate polynomials.
%

% Dutka, AF and FJ Ewens. 1971. A method for improving the accuracy of polynomial
%   regression analysis. J Quality Technology 3:149-155.

% RE Strauss, 2/14/02
%   3/7/03 -  revised sequence of output arguments for orthogpolyregrf().
%   3/17/03 - combined ss, df, and ms stats into single vectors.

function [b,R2a,F,pF,ss,df,ms,alpha,beta] = orthogpolyregr(x,y,degree,doplot)
  if (nargin < 4) doplot = []; end;
  
  if (isempty(doplot))
    doplot = 0;
  end;

  x = x(:);
  y = y(:);
  n = length(x);
  if (n~=length(y))
    error('  ORTHOGPOLYREGR: data vectors are incompatible');
  end;
  
  [b,alpha,beta,P] = orthogpolyregrf(x,y,degree);
  
  ssr = zeros(degree+1,1);
  sse = zeros(degree+1,1);
  R2a = zeros(degree+1,1);
  pR2a = zeros(degree+1,1); 
  Fiof = zeros(degree+1,1);
  pFiof = zeros(degree+1,1);
  
  for i = 1:degree+1                      % Sums of squares
    [e,sse(i),ssr(i),ssto,R2a(i),pR2a(i),pred] = orthogpolyregre(y,b(1:i),P(:,1:i));
  end;
  sse = sse(end);
  ssr = ssr(2:end)-ssr(1:end-1);
  
  dfr = ones(degree,1);
  dfto = n-1;
  dfe = dfto - sum(dfr);
  
  mse = sse/dfe;
  msr = ssr./dfr;
  F = msr/mse;
  pF = 1-fcdf(F,dfr,dfe);
  
  ss = [ssr;sse;ssto];
  df = [dfr;dfe;dfto];
  ms = [msr;mse];
  
  if (doplot)
    plot(x,y,'ko');
    xp = linspace(min(x),max(x))';
    yp = orthogpolyregrp(xp,b,alpha,beta);
    hold on;
    plot(xp,yp,'k');
    hold off;
    putbnd([x;xp],[y;yp]);
  end;
  
  return;