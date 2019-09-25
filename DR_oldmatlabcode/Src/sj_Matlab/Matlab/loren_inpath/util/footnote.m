function hout = footnote(str, posdesc)
% footnote - put a footnote at the bottom of the page.
%
%	FOOTNOTE('text') adds text to the bottom of the figure below
%	all subplots. Use this function after all subplot commands.

% Copyright John Pezaris 4/97, pz@caltech.edu.

if (nargin < 2)
  posdesc = 'right';
end

if strcmp(posdesc, 'left')
  pos = 0.1;
elseif strcmp(posdesc, 'right')
  pos = 0.9;
elseif strcmp(posdesc, 'center')
  pos = 0.5;
end
  
footnote_xpos = 0.4;
footnote_ypos = 0.4;
footnote_pos  = [pos 0 0.1 0.1];

% Fontsize for footnote
fs = get(gcf, 'defaultaxesfontsize') - 4;

old_axes = gca;		% save this for future use

% look for old footnote
footnote_axes = findobj(gcf, 'Tag', 'footnote'); 

% if there, get position info (unless position was specified), 
% and then kill it
if (~isempty(footnote_axes))
  if (nargin < 2)
     footnote_pos = get(footnote_axes, 'pos');
  end
  delete(footnote_axes);
end

new_plot = get(gcf, 'nextplot');
set(gcf, 'nextplot', 'add');
footnote_axes = axes('pos', footnote_pos, 'visible', 'off', 'Tag', 'footnote');
text_handle = text(footnote_xpos, footnote_ypos, str);
set(text_handle, 'horizontalalignment', posdesc, 'fontsize', fs);
set(gcf, 'nextplot', new_plot);

axes(old_axes);		% restore previously selected axis

if nargout > 0
   hout = ht;
end
