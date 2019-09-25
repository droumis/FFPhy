% BINOPT: Fits a binomial distribution B(N,p) to a discrete
%         distribution (X=0,1,2,...,N) by minimizing the G-statistic.
%
%     Usage: [N,p,G,df,prob] = binopt(x,freq)
%
%         x =    vector of integral values of X (of N) represented by
%                  distribution
%         freq = vector of corresponding absolute frequencies
%         N,p =  parameters of best-fitting binomial distribution
%         G =    corresponding G-statistic value
%         df =   corresponding degrees of freedom
%         prob = probability of the observed G-statistic
%

% RE Strauss, 3/25/95

function [N,p,G,df,prob] = binopt(x,freq)
  lenx = length(x);
  xbar = meanwt(x,freq);
  options = [];                 % No foption() changes

  N = x(lenx);             % Initial N at max observed X
  p = fmin('binoptf',0,1,options,N,x,freq);   % Optimize p
  [G,df] = binoptf(p,N,x,freq);               % G-statistic
  prob = 1 - chi2cdf(G,df);     % Corresponding probability

  G_save = G + 1;
  while (G < G_save)
    N_save = N;                 % Save current values
    p_save = p;
    G_save = G;
    df_save = df;
    prob_save = prob;

    N = N + 1;                  % Next N, next values
    p = fmin('binoptf',0,1,options,N,x,freq); % Optimize p
    [G,df] = binoptf(p,N,x,freq);             % G-statistic
    prob = 1 - chi2cdf(G,df);   % Corresponding probability
  end;

  N = N_save;
  p = p_save;
  G = G_save;
  df = df_save;
  prob = prob_save;

  return;
