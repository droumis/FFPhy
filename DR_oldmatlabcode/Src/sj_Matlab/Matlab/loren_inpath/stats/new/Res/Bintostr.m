% BINTOSTR: Converts binary matrix to character equivalent.
%
%     Usage: s = bintostr(b)
%
%           b = binary matrix.
%           ---------------------------------------------------
%           s = character matrix of same size as binary matrix.
%

% RE Strauss, 1/22/00

function s = bintostr(b)
  [r,c] = size(b);

  s = char(b);                            % Allocate output matrix

  for i = 1:r
    for j = 1:c
      if (b(i,j))
        s(i,j) = '1';
      else
        s(i,j) = '0';
      end;
    end;
  end;

  return;
