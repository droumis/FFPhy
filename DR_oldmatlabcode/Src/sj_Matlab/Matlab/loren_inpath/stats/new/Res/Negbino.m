% NEGBINO:  Determines the negative-binomial parameters, and tests for 
%           significant spatial clustering, given a matrix of observed counts 
%           per quadrate.  Uses the parameterization of Anscombe (1949).
%
%     Usage: [k,mu,oP,eP,X2,df,pr] = negbino(counts)
%
%           counts =  matrix of quadrate counts.
%           --------------------------------------
%           k =       aggregation parameter.
%           mu =      mean-count parameter.
%           oP =      observed pdf.
%           eP =      expected pdf.
%           X2 =      chi-squared statistic.
%           df =      degrees of freedom.
%           pr =      chi-squared probability.
%

% RE Strauss, 7/1/99
%   1/4/00 -  changed fminu() to fmins().

function [k,mu,oP,eP,X2,df,pr] = negbino(counts)
  counts = counts(:);                   % Create vector of counts
  if (~isintegr(counts))                % Check for non-integers
    error('NEGBINO: counts must be integers');
  end;

  n = length(counts);
  maxcount = max(counts);               % Create observed pdf
  mu = mean(counts);                    % Estimate of parameter mu

  X = [0:maxcount]';
  oP = zeros(maxcount+1,1);
  [i,f] = uniquef(counts,1);
  oP(i+1) = f;

%  k = fmins('negbinof',1,[],[],mu,oP);   % ML estimate of k
% Problem: lets k become huge

  s2 = var(counts);
  if ((s2-mu)>eps)                        % Variance must exceed the mean
    k = mu.*mu./(s2-mu);                  % MM estimate of k
  else
    k = Inf;
  end;

  pdf = zeros(size(oP));                  % Expected pdf
  pdf(1) = (k/(mu+k)).^k;
  for x = 1:maxcount
    pdf(x+1) = ((k+x-1)/x)*(mu/(mu+k))*pdf(x);
  end;
  eP = n*pdf;

  df = length(oP)-2;
  [X2,df,pr] = goodfit(oP,eP,df);

  i = find(~finite(X2));
  if (~isempty(i))
    X2(i) = zeros(length(i),1);
  end;

  return;
