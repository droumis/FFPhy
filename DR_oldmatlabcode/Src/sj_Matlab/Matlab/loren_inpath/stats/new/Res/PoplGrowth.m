% PoplGrowth: simple birth-death population model.
%
%     Usage:  N = poplgrowth(N0,B,D,tmax)
%
%         N0 = initial population size.
%         B = per capita birth rate.
%         D = per capita death rate.
%         tmax = maximum number of generations.
%         ------------------------------------------
%         N = [tmax+1 x 1] vector of population sizes.
%

function N = poplgrowth(N0,B,D,tmax)
  N = zeros(tmax+1,1);
  N(1) = N0;
  
  for t = 2:tmax+1
    N(t) = N(t-1) + B*N(t-1) - D*N(t-1);
    N(t) = max([0,floor(N(t))]);
  end;
  
  plot(0:tmax,N);
  putxlab('Time');
  putylab('Population size');
  
  return;
  