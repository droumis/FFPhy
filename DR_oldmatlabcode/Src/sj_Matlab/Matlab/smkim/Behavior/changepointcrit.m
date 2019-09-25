
function [cest ] = changepointcrit(I, pb, pvalue)

%Anne Smith, April, 2003
%Finds the position in a binary sequence where the first sub-sequence
%of consecutive ones of length crit occurs. crit is computed
%by calling findj.m

%First, call the routine findj.m to get the number in a row correct responses
%required in this case

  crit = findj(pb, length(I), pvalue)

%Second, move along the binary sequence I and find the first point where
%get crit in a row

realcnt = 0;

cest=NaN;
seq = I;
n   = length(I);
j   = crit;

cnt = 0;
for pp = 1:n-j+1
 
  s1 = sum(seq(pp:pp+j-1));

  if(s1>=j)
   cnt = cnt+1;
   cest = pp;
   break
  end
end
if(cnt>0)
realcnt = realcnt + 1;
end




