% TOFILE : Writes matrix to a specified file.  
%
%     Syntax: tofile(A,'filename',{digits},{append})
%
%           A =        single matrix to be exported.
%           filename = DOS pathname, in single quotes.  The matrix is 
%                        appended to the file, which will be created if it 
%                        does not exist.  File will be written to the current 
%                        directly if a directory path is not specified.
%           digits =   optional number of significant digits (if positive), or 
%                        decimal positions (if negative) [default=+8].
%           append =   optional flag indicating that the matrix is to be 
%                        appended to the file, rather than overwritten 
%                        [default=0].
%

% RE Strauss, 6/23/96
%   9/20/99 - update handling of null input arguments.
%   1/20/00 - allow for writing character matrices.
%   4/20/00 - allow for specification of decimal positions rather than 
%               significant digits.
%   5/16/00 - fixed problem with setting decimal positions in last column.
%   1/28/01 - fixed problem with matrices having cols of all missing data.
%   2/2/01 -  if appending to file, write initial blank line.
%   7/13/01 - added error message for invalid fid.
%   8/29/01 - allow for empty matrix.

function tofile(A,fname,digits,append)
  if (nargin < 3) digits = []; end;
  if (nargin < 4) append = []; end;

  if (isempty(append))
    append = 0;
  end;
  if (isempty(digits))
    digits = 8;
  end;

  fix_dp = 0;
  if (digits < 0)
    digits = -digits;
    fix_dp = 1;
  end;

  [r,c] = size(A);

  if (append)
    fid = fopen(fname,'a');
    fprintf(fid,'\n');
  else
    fid = fopen(fname,'w');
  end;
  if (fid < 0)
    error('  TOFILE: invalid file name (must be in single quotes)');
  end;
  
  if (isempty(A))
    fprintf(fid,'\n');
    fclose(fid);
    return;
  end;

  if (ischar(A))
    for i = 1:r
      fprintf(fid,['%',tostr(c),'s'],A(i,:));
      fprintf(fid,'\n');
    end;
  
    fclose(fid);
    return;
  end;

  if (r > 1)
    maxA = max(abs(A));               % Min & max values per column
    minA = min(A);
  else  % if row vector
    maxA = abs(A);
    minA = A;
  end;

  i = find(abs(maxA) < eps);
  if (~isempty(i))
    maxA(i) = eps * ones(1,length(i));
  end;

  f = max(floor(log10(maxA)+eps)+1,1);    % Number of digits before decimal point
  if (fix_dp)                             % Number of digits after  decimal point
    d = digits*ones(1,length(f));           % Fix dp
  else
    d = max([digits-f;zeros(1,length(f))]); % Fix significant digits
    i = find(maxA<1);                   % If all values < 1, allow full number of digits
    d(i) = d(i)+1;                      %   despite leading zero
  end;

  i = find(minA<0);                   % If any values <0,
  f(i) = f(i)+1;                      %   allow for sign

  w = zeros(1,c);
  editstr = [];

  for j = 1:c                         % For each column,
    B = A(:,j);

    if (isintegr(B))                    % If vector of all integers
      if (all(finite(A(:,j))))          %   truncate fractional portion
        w(j) = f(j);
        editstr = [editstr,' %',int2str(w(j)),'.0f'];
      else
        w(j) = max([3,f(j)]);
        editstr = [editstr,' %',int2str(w(j)),'.0f'];
      end;
      d(j) = 0;
    else                                % If real numbers, full number of digits
      epsilon = 10.^(-d(j)-1);
      for k = 1:d(j)
        BB = B * 10^k;
        e = epsilon * 10^k;
        if (isintegr(BB,e))
          d(j) = k;
          break;
        end;
      end;

      w(j) = f(j)+d(j)+1;
      editstr = [editstr,' %',int2str(w(j)),'.',int2str(d(j)),'f'];
    end;
  end;
  editstr = [editstr,'\n'];

  for i = 1:r
    if (all(finite(A(i,:))))
      fprintf(fid,editstr,A(i,:));
    else
      for j = 1:c
        if (finite(A(i,j)))
          if (d(j) == 0)
            estr = [' %',int2str(w(j)),'.0f'];
          else
            estr = [' %',int2str(w(j)),'.',int2str(d(j)),'f'];
          end;
          fprintf(fid,estr,A(i,j));
        else
          wj = w(j)-3;
          if (wj==0)
            estr = ' NaN';
          else
            estr = [' '];
            for k = 1:floor(wj/2)
              estr = [estr, ' '];
            end;
            estr = [estr,'NaN'];
            for k = 1:ceil(wj/2)
              estr = [estr, ' '];
            end;
          end;
          fprintf(fid,estr);
        end;
      end;
      fprintf(fid,'\n');
    end;
  end;

  fclose(fid);

  return;
