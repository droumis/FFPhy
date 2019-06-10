%FIND_NEARBY_MEX Find elements in target vector that are nearby each element of a reference vector.
%
%   IDX = FIND_NEARBY_MEX(X,Y,LOWER_BOUND,UPPER_BOUND) takes two column vectors
%   of real non-NaN double elements X and Y, and two double real scalars
%   LOWER_BOUND and UPPER_BOUND, and returns a column cell array IDX containing
%   index vectors. IDX is the same size as X, such that
%
%   IDX{i} = find((Y >= X(i) - LOWER_BOUND) & (Y < X(i) + UPPER_BOUND));
%   all(Y(IDX{i}) - X(i) < UPPER_BOUND)
%   all(Y(IDX{i}) - X(i) >= LOWER_BOUND)
%
%   This MEX function depends on standard C libraries.
%
%Written by SMK, 2009 May 29.
%
