function h = mplot(matrix, varargin)

% MPLOT  Plot rows of a data matrix.
%    MPLOT(MATRIX) makes a line plot of the rows of an (M x N) matrix using
%    only one line object.  For large M, these plots can be much faster than
%    PLOT(MATRIX') since only one object is created and that object does not
%    require front-to-back sorting.  The tradeoff is that all data rows are
%    plotted with a single linestyle (color, linewidth, etc).
% 
%    MPLOT(MATRIX, ...) passes additional arguments through to PLOT.
%
%    H = MPLOT(MATRIX) returns a handle to the line object.
%
% Last Modified: sbm, 8/12/03

[M,N] = size(matrix);

matrix = padmatrix(matrix, [0 1 0 0], NaN)';
inds = repmat([1:N NaN]', [M,1]);

h = plot(inds, matrix(:), varargin{:});

if (nargout < 1)
    clear h;
end