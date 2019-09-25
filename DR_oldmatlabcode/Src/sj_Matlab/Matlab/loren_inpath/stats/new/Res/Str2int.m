% Str2int:  Performs the inverse operation from the Matlab int2str function, converting
%           the rows of a string matrix into integers.  The only allowable characters
%           are leading blanks, the digits 0-9, and the + and - signs; other characters 
%           are ignored.
%
%     Usage: intvals = str2int(charvals)
%
%         charvals = [r x c] string matrix.
%         -------------------------------------------------
%         intvals =  [r x 1] vector of integer equivalents.
%

% RE Strauss, 3/8/02

function intvals = str2int(charvals)
  intvals = [];
  if (~ischar(charvals))
    error('  STR2INT: input must be a character matrix.');
  end;
  
  [r,c] = size(charvals);
  intvals = zeros(r,1);
  
  digits = '0123456789';
  
  for ir = 1:r                          % Cycle through rows of string matrix
    c = charvals(ir,:);                   % Isolate current string
    iblank = find(isin(c,' '));             % Find blanks
    idigit = find(isin(c,digits));          % Find digits
    isign =  find(isin(c,'-+'));            % Find sign
    
    if (isempty(iblank))
      iblank = 0;
    end;    
    if (isempty(isign))
      isign = 0;
    end;
    if (isempty(idigit))
      idigit = 0;
    end;
    
    err = 0;                                % Check for error conditions
    if (any(~isin(c,[' -+',digits])))         % Invalid character
      err = 1;
    end;
    if (max(isign)>=min(idigit))              % Sign to right of first digit
      err = 1;
    end;
    if (isin(' ',c(min(idigit):max(idigit)))) % Blank in midst of digits
      err = 1;
    end;
    if (length(isign)>1)                      % Multiple signs
      err = 1;
    end;
    
    if (~err)
      s = '+';
      if (isign>0)
        s = c(isign);                         % Stash sign
      end;
      c = c(idigit);                          % Reduce to digits
      nc = length(c);
      n = 0;
      for i = nc:-1:1
        [b,j] = isin(c(i),digits);
        n = n + (j-1)*10^(nc-i);
      end;
      if (s=='-')
        n = -n;
      end;
      intvals(ir) = n;
    else                          
      intvals(ir) = NaN;
    end;
  end;

  return;
  