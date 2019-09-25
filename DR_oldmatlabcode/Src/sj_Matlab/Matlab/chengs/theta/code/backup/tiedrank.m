function [r, tieadj] = tiedrank(x)
%TIEDRANK Compute the ranks of a sample, adjusting for ties
%   [R, TIEADJ] = TIEDRANK(X) computes the ranks of the values in
%   the vector X.  If any X values are tied, TIEDRANK computes
%   their average rank.  The return value TIEADJ is an adjustment
%   for ties required by the nonparametric tests SIGNRANK and
%   RANKSUM.

%   Tom Lane, 4-16-99
%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.5 $  $Date: 2002/01/17 21:32:04 $


[rowx, colx] = size(x);

if min(rowx, colx) ~= 1
   error('TIEDRANK requires a vector.');
end 
if rowx == 1
   x = x';
end
   
[sx, rowidx] = sort(x);
ranks = 1:length(x);

% Adjust for ties
tieloc = find(~diff(sx));
tieadj = 0;
while (length(tieloc) > 0)
   tiestart = tieloc(1);
   ntied = 1 + sum(sx(tiestart) == sx(tiestart+1:end));
   tieadj = tieadj + ntied*(ntied-1)*(ntied+1)/2;
   ranks(tiestart:tiestart+ntied-1) = tiestart + (ntied-1)/2;
   tieloc(1:ntied-1) = [];
end

r(rowidx) = ranks;
