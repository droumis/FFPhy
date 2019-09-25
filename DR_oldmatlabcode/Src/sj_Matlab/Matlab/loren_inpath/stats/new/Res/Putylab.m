% PUTYLAB: Adds Ylabel text to the current figure, with default font 
%          characteristics.  If an optional percentage value is provided as a second 
%           argument, it is printed (in parentheses) following the text string.
%
%     Syntax: putylab(textstring,{percent},{fontsize})
%
%         textstring - character string to be printed.
%         percent -    optional "percentage account for" value.
%         fontsize -   optional font size [default = 14].
%

% RE Strauss
%   3/23/98 - add optional percent string
%   11/5/98 - add optional font size

function putylab(textstring,percent,fontsize)
  if (nargin < 2) percent = []; end;
  if (nargin < 3) fontsize = []; end;

  if (~isempty(percent))
    p = sprintf(' (%3.1f',percent);
    textstring = [textstring p '%)'];
  end;

  if (isempty(fontsize))
    fontsize = 14;
  end;

  YLabel(textstring,'FontName','helvetica','FontSize',fontsize);
  return;

