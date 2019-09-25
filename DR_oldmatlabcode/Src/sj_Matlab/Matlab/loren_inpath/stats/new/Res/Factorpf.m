% FACTORPF: Objective function for FACTORP.  Given a column vector of loadings for
% the primary factor and a boolean matrix of cells to be sequestered, calculates
% the sum of squared deviations of observed covariances from predicted
% (reproduced) covariances.
%

function sse = factorpf(f,c,skip)
  diagc = diag(c);
  fsq = f.^2;
  neg_indx = find(fsq > diagc);         % Negative diagonal residuals

  if (length(neg_indx)>0)               % Penalize loadings producing negative
    totneg = sum(fsq(neg_indx)-diagc(neg_indx)); %   diagonal elements
    sse = 1e6 * totneg;                 
  else
    rc = f*f';                          % Reproduced covars
    sse = sum((c(~skip)-rc(~skip)).^2); % SSE of non-squestered cells
  end;

  return;
