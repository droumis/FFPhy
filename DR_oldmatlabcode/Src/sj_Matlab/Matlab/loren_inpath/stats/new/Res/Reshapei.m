% RESHAPEI: Complements the reshape() function by mapping the row/col indices 
%           (subscripts) from the original matrix to the reshaped matrix.
%
%     Usage: ind2 = reshapei(dim1,dim2,ind1)
%
%         dim1 = 2-element vector containing the dimensions of the input matrix. 
%                  If the input matrix is a vector, dim1 can be a scalar 
%                  indicating its length.
%         dim2 = 2-element vector containing the dimensions of the output matrix.
%                  If the output matrix is a vector, dim1 can be a scalar 
%                  indicating its length.
%         ind1 = 2-col matrix whose rows contain sets of index (subscript) 
%                  pairs i,j for the input matrix; can be a vector if the input 
%                  matrix is a vector.
%         -----------------------------------------------------------------------
%         ind2 = 2-col matrix whose rows contain sets of index (subscript) 
%                  pairs i,j for the output matrix, or a column vector if the
%                  output matrix is a vector.
%

% RE Strauss, 11/23/99

function ind2 = reshapei(dim1,dim2,ind1)
  if (prod(dim1) ~= prod(dim2))
    error('  RESHAPEI: inconsistent dimensions for input and output matrices.');
  end;
  if (any(max(ind1)>max(dim1)))
    error('  RESHAPEI: subscripts out of range.');
  end;

  if (length(dim1)>1 & min(dim1)>1)             % If input matrix is not a vector,
    ind1 = (ind1(:,1)-1).*dim1(1) + ind1(:,2);  %   convert indices to vector indices
    dim1 = length(ind1);
  end;

  if (length(dim2)==1 | min(dim2)==1)           % If output matrix is a vector,
    ind2 = ind1;                                %   copy input indices
  else                                          % else
    ind2 = zeros(dim1,2);                       %   allocate output matrix
    ind2(:,1) = ceil(ind1/dim2(2));             %   row subscripts
    ind2(:,2) = ind1 - (ind2(:,1)-1)*dim2(2);   %   col subscripts
  end;

  return;

