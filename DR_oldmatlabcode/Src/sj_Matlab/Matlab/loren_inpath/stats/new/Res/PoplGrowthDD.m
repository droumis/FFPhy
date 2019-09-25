% PoplGrowthDD: simple birth-death population model.
%
%     Usage:  N = poplgrowthdd(N0,b0,b,d0,d,tmax)
%
%         N0 = initial population size.
%         b0 = per capita birth rate at population size zero.
%         b = change in birth rate per unit increase in N.
%         d0 = per capita death rate at population size zero.
%         d = change in death rate per unit increase in N.
%         tmax = maximum number of generations.
%         ---------------------------------------------------
%         N = [tmax+1 x 1] vector of population sizes.
%

function N = poplgrowth(N0,b0,b,d0,d,tmax)
  N = zeros(tmax+1,1);
  N(1) = N0;
  
  for t = 2:tmax+1
    B = b0 - b*N(t-1);
    D = d0 + d*N(t-1);
    N(t) = N(t-1) + B*N(t-1) - D*N(t-1);
  end;
  
  plot(0:tmax,N);
  putxlab('Time');
  putylab('Population size');
  
  return;
  