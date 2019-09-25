% SUBMATROWS: Reduce a set of corresponding matrices to identical subsets of 
%             rows.
%
%     Usage: [s1,...,s9] = submatrows(i,m1,...,m9)
%
%         i = vector of indices of rows to be included in subsets.
%         m1,...,m9 = list of matrices to be subsetted.
%
%         s1,...,s9 = corresponding list of subsetted matrices
%

% RE Strauss, 6/3/01
%   6/14/01 - corrected problem with saving rows of matrix.
%   6/16/01 - added error messages.

function [s1,s2,s3,s4,s5,s6,s7,s8,s9] = submatrows(i,m1,m2,m3,m4,m5,m6,m7,m8,m9)
  if (nargin > 1)
    if (isvector(m1))
      matlen = length(m1);
    else
      matlen = size(m1,1);
    end;
  end;
  imax = max(i);

  if (nargin>1) 
    if (isvector(m1))
      mlen = length(m1);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s1 = m1(i); 
    else 
      mlen = size(m1,1);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s1 = m1(i,:); 
    end; 
  end;

  if (nargin>2) 
    if (isvector(m2))
      mlen = length(m2);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s2 = m2(i); 
    else 
      mlen = size(m2,1);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s2 = m2(i,:); 
    end; 
  end;

  if (nargin>3) 
    if (isvector(m3))
      mlen = length(m3);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s3 = m3(i); 
    else 
      mlen = size(m3,1);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s3 = m3(i,:); 
    end; 
  end;

  if (nargin>4) 
    if (isvector(m4))
      mlen = length(m4);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s4 = m4(i); 
    else 
      mlen = size(m4,1);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s4 = m4(i,:); 
    end; 
  end;

  if (nargin>5) 
    if (isvector(m5))
      mlen = length(m5);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s5 = m5(i); 
    else 
      mlen = size(m5,1);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s5 = m5(i,:); 
    end; 
  end;

  if (nargin>6) 
    if (isvector(m6))
      mlen = length(m6);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s6 = m6(i); 
    else 
      mlen = size(m6,1);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s6 = m6(i,:); 
    end; 
  end;

  if (nargin>7) 
    if (isvector(m7))
      mlen = length(m7);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s7 = m7(i); 
    else 
      mlen = size(m7,1);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s7 = m7(i,:); 
    end; 
  end;

  if (nargin>8) 
    if (isvector(m8))
      mlen = length(m8);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s8 = m8(i); 
    else 
      mlen = size(m8,1);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s8 = m8(i,:); 
    end; 
  end;

  if (nargin>9) 
    if (isvector(m9))
      mlen = length(m9);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s9 = m9(i); 
    else 
      mlen = size(m9,1);
      if (mlen ~= matlen)
        error('  SUBMATROWS: all input matrices not of compatible size.');
      end;
      if (imax > mlen)
        error('  SUBMATROWS: subscript out of bounds.');
      end;
      s9 = m9(i,:); 
    end; 
  end;

  return;