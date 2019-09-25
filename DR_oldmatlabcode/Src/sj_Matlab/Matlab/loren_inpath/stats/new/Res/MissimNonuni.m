% MISSIMNONUNI: Simulate non-uniform distribution of missing values in matrix.

N = 30;
P = 10;                     % Must be even values
propmiss = 0.10;
iter = 200;

suites = 2*ones(1,P/2);     % 2 vars/suite
nsuites = length(suites);

tofileappend = 0;

for in = 1:length(N)
  n = N(in);
  
  for ip = 1:length(P)
    p = P(ip);
    
    for ipm = 1:length(propmiss)
      pmiss = propmiss(ipm);
      
      for nsuitemiss = 1:nsuites
        tofile([n,p,pmiss,nsuitemiss],'Missim.txt',-4,tofileappend);
        [mse,h,mac,e,deltavar,reldevrand] = missim(n,suites,pmiss,iter,0,[],nsuitemiss);
        tofileappend = 1;
      end;
    end;
  end;
end;
