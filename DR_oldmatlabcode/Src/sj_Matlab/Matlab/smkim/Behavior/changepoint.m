function [cest] = changepoint(I)

%Change-point test for binary observations
%Siegel and Castellan, p65, 1988

%takes a binary sequence of zeros and ones and
%estimates where the distribution changes

N  = length(I);   %need N>25 for this test to be good             

sj = cumsum(I);
m  = sum(I);
n  = N - m;

for j = 1:N-1

 v      = N/m/n*(sj(j) - j*m/N);
 dmn(j) = abs(v);

end

 crit          = 1.36*sqrt(N/m/n);
 [maxdmn, jj]  = max(dmn);

 if(maxdmn > crit)
  cest = jj;
 else
  cest = NaN;
 end



