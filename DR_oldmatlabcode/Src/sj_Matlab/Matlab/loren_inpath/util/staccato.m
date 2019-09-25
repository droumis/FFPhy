function d = staccato (sz, times)
% staccato - create column vector(s) of zeros with occasional ones
%	     
% d = staccato (size, ones)
% 
% size  = overall size of the vector
% ones  = vector of indices at which to place the ones
%
% If size is a 2 element vector, d is a matrix and the second element of size 
% is the number of (identical) columns in d.

d = zeros (sz(1), 1);
d(times) = ones (max(size(times)), 1);

if max(size(sz)) > 1
  d = d * ones(1, sz(2));
end

