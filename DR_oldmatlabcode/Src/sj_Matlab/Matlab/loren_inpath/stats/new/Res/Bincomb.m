% BINCOMB:  Returns, as rows of a binary matrix, all combinations of 0s and 1s 
%           [the binary values of decimal 0 through (2^ncols - 1)].
%
%     Usage: binvals = bincomb(ncols)
%
%         ncols = number of columns of output matrix.
%         ----------------------------------------------------------
%         binvals = [2^ncols x ncols] matrix of binary combinations,  
%                     from 000... to 111...
%

% RE Strauss, 7/20/00

function binvals = bincomb(ncols)
  nrows = 2^ncols;
  binvals = zeros(nrows,ncols);

  g = 1;
  for p = ncols:-1:1
    for i = 1:(2*g):nrows
      binvals(i:(i+2*g-1),p) = [zeros(g,1); ones(g,1)];
    end;
    g = 2*g;
  end;

  return;
