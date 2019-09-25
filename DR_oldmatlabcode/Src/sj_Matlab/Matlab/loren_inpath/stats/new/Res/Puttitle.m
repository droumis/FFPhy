% PUTTITLE: Adds title text to the current figure, with default font 
%           characteristics.
%
%     Syntax: puttitle(textstring,{fontsize})
%
%           textstring - character string to be printed above figure.
%           fontsize -   optional font size [default = 16].
%

% RE Strauss, 12/8/96
%   11/5/98 - added optional font size specification.

function puttitle(textstring,fontsize)
  if (nargin < 2) fontsize = []; end;

  if (isempty(fontsize))
    fontsize = 16;
  end;

  title(textstring,'FontName','helvetica','FontSize',fontsize);
  return;

