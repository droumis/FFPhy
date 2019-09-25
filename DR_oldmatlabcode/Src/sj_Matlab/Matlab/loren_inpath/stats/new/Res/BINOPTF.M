% BINOPTF: Objective function for BINOPT.  Finds the G-statistic
%          (Williams-corrected) for the agreement between an observed
%          discrete distribution and a binomial distribution with
%          parameters B(N,p).  Combines cells with small freqs
%          (exp<5, obs<1) into single cell.
%
%     Usage: [G,df,e,f] = binoptf(p,N,x,f)
%

% RE Strauss, 3/23/95
%   9/20/99 - update handling of null input arguments.
%   1/25/99 - renamed and changed usage of unique().

function [G,df] = binoptf(p,N,x,f)
  len = length(x);
  df = len-1-2;                     % Intrinsic hypothesis, 2 params
  n = sum(f);                       % Total sample size
  e = n * binopdf(x,N,p);           % Expected binomial counts

  % Check for and combine cells w/ small freqs
  omit = uniquef([find(e<5), find(f<1)],1);
  if (~isempty(omit))
    ecomb = sum(e([omit]));         % Combine small cells
    fcomb = sum(f([omit]));

    indx = ones(1,len);             % Find complement of omission set
    indx(omit) = indx(omit)-1;
    incl = find(indx);

    e = [e(incl), ecomb];           % Append combined cell to
    f = [f(incl), fcomb];           %   non-small ones
    df = length(e)-1-2;             % Adjust df
  end;

  q = 1 + ((length(e)^2-1)/(6*n*max(df,1)));
  G = 2 * sum(f .* log(f./e)) / q;  % G-statistic

  return;
