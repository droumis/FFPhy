% POISSON:  Determines the Poisson parameter, and tests for significant 
%           spatial clustering, given a matrix of observed counts 
%           per quadrate.
%
%     Usage: [lambda,oP,eP,X2,df,pr] = poisson(counts)
%
%           counts =  matrix of quadrate counts.
%           --------------------------------------
%           lambda =  Poisson parameter.
%           oP =      observed pdf.
%           eP =      expected pdf.
%           X2 =      chi-squared statistic.
%           df =      degrees of freedom.
%           pr =      chi-squared probability.
%

% RE Strauss, 7/1/99

function [lambda,oP,eP,X2,df,pr] = poisson(counts)
  counts = counts(:);                   % Create vector of counts
  if (~isintegr(counts))                % Check for non-integers
    error('POISSON: counts must be integers');
  end;

  n = length(counts);
  maxcount = max(counts);               % Create observed pdf

  X = [0:maxcount];
  oP = zeros(maxcount+1,1);             % Observed pdf
  [i,f] = uniquef(counts,1);
  oP(i+1) = f;

  [lambda,eP,X2,df,pr] = poisfit(oP,X);   % Estimate lambda

  return;
