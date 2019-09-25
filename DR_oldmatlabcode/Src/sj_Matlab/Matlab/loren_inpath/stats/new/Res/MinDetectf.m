% MINDETECTF: Estimates probability of detecting a specified relative frequency 
%             for a specified sample size.
%
%     Usage: pfind = mindetectf(p0,N,iter)
%
%         p0 =    relative frequency of occurrence of event.
%         N =     sample size.
%         iter =  number of iterations.
%         --------------------------------------------------
%         pfind = probability of detection of event.
%

% RE Strauss, 5/4/01

function pfind = mindetectf(p0,N,iter)
  r = binornd(N,p0,iter,1);
  pfind = sum(r>0)/iter;

  return;
