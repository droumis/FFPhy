% KS1CF: Objective function for KSTEST1C, called by BOOTSTRP.
%
%     Syntax: solution = ks1cf(X,nu1,nu2,nullflag,pp,stat,rnge)
%
%         X =         [n x 1] vector of data scores.
%         nu1,nu2 =   variables passed by BOOTSTRP but unused.
%         nullflag =  flag indicating that null distribution is being generated.
%         pp    =     row vector representing the piecewise-polynomial form of
%                       the integral of the continuous reference probability distribution.
%         stat =      flag indicating the statistic to be calculated:
%                       0 = MSE statistic
%                       1 = KS statistic
%         rnge =      2-element vector specifying min & max of original data.
%         ----------------------------------------------------------------------
%         solution =  test-statistic value.
%

% RE Strauss, 9/25/96

function solution = ks1cf(X,nu1,nu2,nullflag,pp,stat,rnge)
  MSE = 0; KS = 1; X2 = 2;              % Statistic indicators

  y = [0:1/(length(X)-1):1]';           % Cumulative increments

  if (~nullflag)                        % Compare bootstrapped sample & reference distrib
    X = sort(X);                          % X = sorted bootstrapped sample
    t = spintgrl(pp,[rnge(1);X;rnge(2)]); % Evaluate integral at data points
    t = t - t(1);                         % Anchor to 0 at left
    t = t / t(length(t));                 % Anchor to 1 at right
    t(1) = [];                            % Delete original min & max
    t(length(t)) = [];
  end;

  if (nullflag)                         % Compare data & random sample from reference distrib
    r = rand(length(X),1);                % Random sample of same size as data
    r = sort(r*(rnge(2)-rnge(1)) - rnge(1)); % Scale to same possible range as data

    t = spintgrl(pp,[rnge(1);r;rnge(2)]); % Evaluate integral at random points
    t = t - t(1);                         % Anchor to 0 at left
    t = t / t(length(t));                 % Anchor to 1 at right
    t(1) = [];                            % Delete original min & max
    t(length(t)) = [];
  end;

  diff = y - t;                         % Difference between distribs
  if (stat==KS)                           % Calc KS statistic
    solution = max(abs(diff));
  elseif (stat==MSE)                      % Calc MSE statistic
    solution = mean(diff.*diff);
  end;

  return;


