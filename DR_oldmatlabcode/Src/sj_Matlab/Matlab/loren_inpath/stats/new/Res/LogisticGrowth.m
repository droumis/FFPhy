% LogisticGrowth: logistic birth-death population model.
%
%     Usage:  N = poplgrowth(N0,r,{sr},K,{sK},tmax,{harvest})
%
%         N0 = initial population size.
%         r = intrinsic growth rate.
%         sr = optional standard error of r for stochastic variation
%               [default = 0].
%         K = asymptotic population size.
%         sK = standard error of K for stochastic variation
%               [default = 0].
%         tmax = maximum number of generations.
%         harvest = optional proportion of population harvested per 
%             time unit [default = 0].
%         -----------------------------------------------------------
%         N = [tmax+1 x 1] vector of population sizes.
%

function N = poplgrowth(N0,r,sr,K,sK,tmax,harvest)
  if (nargin < 3)
    sr = [];
  end;
  if (nargin < 5)
    sK = [];
  end;
  if (nargin < 7)
    harvest = [];
  end;
  
  if (isempty(sr))
    sr = 0;
  end;
  if (isempty(sK))
    sK = 0;
  end;
  if (isempty(harvest))
    harvest = 0;
  end;

  N = zeros(tmax+1,1);
  N(1) = N0;
  
  for t = 2:tmax+1
    rr = r + sr*randn;
    rK = K + sK*randn;
    dN = rr*N(t-1)*((rK-N(t-1))/rK);
    N(t) = N(t-1) + dN;
    N(t) = N(t) - harvest*N(t);
    N(t) = max([0,floor(N(t))]);
  end;
  
  plot(0:tmax,N);
  putxlab('Time');
  putylab('Population size');
  
  return;
  