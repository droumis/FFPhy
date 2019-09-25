%FIND_NEAREST_MEX Find element in target vector closest to specified offset relative to each element of reference vector.
%
%   [IDX, DELTA] = FIND_NEAREST_MEX(X,Y,OFFSET) takes two column vectors of real
%   non-NaN double elements X and Y, and a double real scalar OFFSET, and
%   returns column vectors IDX and DELTA of the same size as X, such that
%
%   DELTA = X + OFFSET - Y(IDX)
%   abs(DELTA(i)) = min(abs(X(i) + OFFSET - Y)) for all i = 1:numel(X)
%
%   If the minimum is not unique, then IDX(i) is chosen to be as close to i as
%   possible (which is useful behavior when X is identical to Y).
%
%   This MEX function depends on standard C libraries.
%
%Written by SMK, 2009 November 25.
%
