% KS1DF: Objective function for KSTEST1D, called by BOOTSTRP.
%
%     Usage: solution = ks1df(X,not_used1,not_used2,tdist,stat,rnge)
%
%         Xf =        [n x 1] vector of data scores.
%         not_used1 = variable passed by BOOTSTRP (grps=[]).
%         not_used2 = variable passed by BOOTSTRP (initial_solution).
%         tdist =     row vector representing the discrete theoretical 
%                       cumulative distribution.
%         stat =      flag indicating the statistic to be calculated:
%                       0 = MSE statistic,
%                       1 = KS statistic,
%         rnge =      2-element vector specifying min & max of original data.
%         ----------------------------------------------------------------------
%         solution =  test-statistic value.
%

% RE Strauss, 9/25/96

function solution = ks1df(Xf,not_used1,not_used2,tdist,stat,rnge)
  MSE = 0;  KS = 1;

  X = sort(X);                          % X = sorted data sample
  x = [0:1/(length(X)-1):1]';           % Cumulative increments

  lenx = length(X);
  t = zeros(lenx,1);

  for i = 1:lenx
    j = max(1,max(find(tdist(:,2) <= X(i))));
    t(i) = tdist(j,1);
  end;

  diff = x - t;                         % Difference between distribs
  if (stat==KS)                           % Calc KS statistic
    solution = max(abs(diff));
  elseif (stat==MSE)                      % Calc MSE statistic
    solution = mean(diff.*diff);
  elseif (stat==X2)                       % Calc X2 statistic
    solution = sum(diff.^2./t);
  end;

  return;


