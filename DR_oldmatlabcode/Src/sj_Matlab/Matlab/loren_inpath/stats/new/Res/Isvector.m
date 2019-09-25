% ISVECTOR: Returns 1 if the input matrix is a vector, 0 if not.  Optionally 
%           returns the vector's length (=total number of cells, if not a 
%           vector), and a boolean flag indicating whether the vector is a 
%           colummn.
%
%     Usage: [isvect,ncells,iscol] = isvector(X)
%
%           X =       [r x c] matrix.
%           --------------------------------------------------------------------
%           isvect =  boolean flag indicating whether (=1) or not (=0) the input 
%                       matrix is a vector.
%           ncells =  number of cells.
%           iscol =   boolean flag indicating whether (=1) or not (=0) the 
%                       vector is a column vector.
%

% RE Strauss, 11/20/99
%   5/4/00 - corrected wrong decision for scalar input.

function [isvect,ncells,iscol] = isvector(X)
  [r,c] = size(X);

  if ([r c]==[1 1])                     % Is a scalar
    isvect = 0;
    ncells = 1;
    iscol = 0;
  elseif (min([r,c])==1)                % Is a vector
    isvect = 1;
    ncells = max([r,c]);
    iscol = 0;
    if (r>1)
      iscol = 1;
    end;
  else                                  % Is not a vector
    isvect = 0;
    ncells = r.*c;
    iscol = 0;
  end;

  return;