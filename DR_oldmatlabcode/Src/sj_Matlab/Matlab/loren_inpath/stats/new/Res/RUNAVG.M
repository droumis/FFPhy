% RUNAVG: Sweeps an input matrix by column, producing a matrix of column
%         vectors in which each element has been replaced by its
%         m-point running average.  Terminal elements are replaced by the
%         nearest running-average values.
%            'm' should be an odd number representing the number of elements
%         within the window centered on the current vector element.
%
%     Usage: vectout = runavg(vectin,m)
%
%         vectin =  input matrix, of arbitrary size.
%         m =       width of averaging window [default = 3].
%         --------------------------------------------------
%         vectout = smoothed output matrix.
%

% RE Strauss, 6/27/93
%   7/3/02 - replace terminal values by nearest running average;
%            use iseven() function.

function vectout = runavg(vectin,m)
  if (nargin < 2) m = []; end;

  if (isempty(m))                 % Default window-width
    m = 3;
  end;

  [r,c] = size(vectin);
  if (r==1),                      % If input is row vector,
    vectin = vectin';               % Transpose
    rowvector = 1;                  % Set flag
    [r,c] = size(vectin);           % Resize
  else
    rowvector = 0;
  end;

  if (m > r)                      % If m too large, exit and return input
    vectout = vectin;             %    matrix (in case called from function)
    error('  RUNAVG: Averaging window too large for input matrix');
  end;

  vectout = [];                   % Begin with null output matrix
  for col = 1:c,                  % For all columns,
    if (iseven(m))                  % If m is even,
      m = m+1;                      %   increment to next odd value
    end;
    k = floor(m/2);                 % Number of neighboring elements
                                    %   within window
    outvector = vectin(:,col);      % Copy input column to output vector
    for i = (k+1):(r-k),            % Sweep input vector
      avg = 0;
      for j = (i-k):(i+k),            % Calculate running average for
        avg = avg + vectin(j,col)/m;  %   current element
      end;
      outvector(i) = avg;
    end;
    outvector(1:k,col) = outvector(k+1,col);  % Overwrite terminal elements
    outvector((r-k+1):end,col) = outvector(r-k,col);
    vectout = [vectout outvector];  % Append smoothed vector to output
  end;

  if (rowvector),                   % If input was row vector,
    vectout = vectout';               % Re-transpose
  end;

  return;

