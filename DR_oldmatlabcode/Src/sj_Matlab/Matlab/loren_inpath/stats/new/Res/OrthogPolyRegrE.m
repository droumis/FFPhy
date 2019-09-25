% OrthogPolyRegrE: Returns residuals and sums of squares for orthogonal polynomials.
%
%     Usage: [e,sse,ssr,ssto,R2a,pr] = orthogpolyregre(y,b,P)
%
%         y = column vector of dependent variable (original, not just unique
%                   values of x).
%         b =     column vector of regression coefficients.
%         P =     values of observations (rows) for each orthogonal polynomial term 
%                   (cols).
%         -------------------------------------------------------------------------
%         e =     column vector of residuals.
%         sse =   residual sum of squares.
%         ssr =   regression sum of squares.
%         ssto =  total sum of squares.
%         R2a =   adjusted squared multiple correlation.
%         pr =    significance level of R2a.
%

% RE Strauss, 2/8/02

function [e,sse,ssr,ssto,R2a,pr,p] = orthogpolyregre(y,b,P)
  p = P*b;                            % Predicted values at unique x's
  e = y - p;                          % Residuals
  sse = e'*e;                         % SSE
  d = y - b(1);                       % Deviations from y-mean
  ssto = d'*d;                        % SSTO
  ssr = ssto-sse;                     % SSR
  [R2a,F,pr,df] = r2adj(ssr,ssto,length(y),length(b));

  return;
  