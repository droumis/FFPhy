% REPLACE: Specifies replacement of numeric elements in a matrix.  Target items 
%            may be exact or within specified bounds.
%
%     Usage: outmat = replace(inmat,what,by)
%
%         inmat = [r x c] input matrix.
%         what  = vector (length n) or [n x 2] matrix of values to be replaced; 
%                   if a matrix, the first and second columns represent lower 
%                   and upper bounds on values to be replaced.
%         by =    corresponding vector (length n) of replacement values.
%         -------------------------------------------------------------------
%         outmat = matching [r x c] output matrix.
%

% RE Strauss, 12/21/99

function outmat = replace(inmat,what,by)
  outmat = inmat;

  rng = 0;
  [r,c] = size(what);
  if (min([r,c])>1)
    if (c>2)
      error('  REPLACE: "what" matrix must have one or two columns');
    else
      rng = 1;
      n = max([r,c]);
    end;
  else
    n = length(what);
  end;

  by = by(:);
  if (length(by) ~= n)
    error('  REPLACE: "what" and "by" matrices incompatible');
  end;

  for ir = 1:n
    if (rng)
      [i,j] = find(inmat>=what(ir,1) & inmat<=what(ir,2));
    else
      [i,j] = find(inmat==what(ir));
    end;
    if (~isempty(i))
      for k = 1:length(i)
        outmat(i(k),j(k)) = by(ir);
      end;
    end;
  end;

  return;
