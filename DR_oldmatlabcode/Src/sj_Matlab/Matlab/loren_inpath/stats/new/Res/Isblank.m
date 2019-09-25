% ISBLANK: Determines whether or not a character matrix consists of all blanks.
%
%     Usage: b = isblank(C)
%
%         C = character matrix of arbitrary size.
%         ---------------------------------------------------------------------
%         b = boolean value indicating whether C does (=1) or does not (=0) 
%               consist of all blanks.
%

% RE Strauss, 7/11/00

function b = isblank(C)
  if (~ischar(C))
    error('  ISBLANK: input must be a character matrix.');
  end;

  C = C(:)';

  b = 0;
  if (C == blanks(length(C)))
    b = 1;
  end;

  return;
