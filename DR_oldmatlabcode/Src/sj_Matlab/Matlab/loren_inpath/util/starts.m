function [y, yi] = starts(x, len, del)
% starts - find sequences of entries in a vector
%
% [y, yi] = starts(x) returns a vector of entries and indices in x.  
% Each index marks the beginning of a sequence of numbers, each one 1 
% greater than the previous. 
% 
% [y, yi] = starts(x, min, del) returns a vector of entries and
% indices of sequences in x, that have at least min adjacent elements
% separated by delta. (x, x+delta, x+2*delta ... x + (min-1)*delta ...)
%
% y = starts(find(x)) finds the beginnings of sequences of non zero 
% elements in x. 

if nargin < 2  len = 1; end
if nargin < 3  del = 1; end

y = [];
yi = [];

c = 1;			% candidate index
n = 1;			% number in sequence

for i = 2:length(x)
  if x(i) == (x(i-1) + del)
    n = n+1;
  else
    if n >= len
      y  = [y, x(c)];
      yi = [yi, c];
    end
    c = i;
    n = 1;
  end
end

if n >= len
  y = [y, x(c)];
  yi = [yi, c];
end


