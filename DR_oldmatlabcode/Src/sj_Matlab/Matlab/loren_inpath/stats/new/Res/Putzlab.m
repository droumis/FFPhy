% PUTZLAB: Adds Ylabel text to the current figure, with default font 
%          characteristics.  If an optional percentage value is provided as a second 
%           argument, it is printed (in parentheses) following the text string.
%
%     Syntax: putzlab(textstring,{percent},{fontsize})
%
%         textstring - character string to be printed.
%         percent -    optional "percentage account for" value.
%         fontsize -   optional font size [default = 14].
%

% RE Strauss, 11/10/02, modified from putylab().

function putzlab(textstring,percent,fontsize)
  if (nargin < 2) percent = []; end;
  if (nargin < 3) fontsize = []; end;

  if (~isempty(percent))
    p = sprintf(' (%3.1f',percent);
    textstring = [textstring p '%)'];
  end;

  if (isempty(fontsize))
    fontsize = 14;
  end;

  ZLabel(textstring,'FontName','helvetica','FontSize',fontsize);
  return;

