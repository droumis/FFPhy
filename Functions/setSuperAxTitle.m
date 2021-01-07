

%{
set super figure axis title i.e. above and independent of all the subplots

arg1: str: title string

varg-vertPos: float: % to top
varg-horzPos: float: % to right
varg-titleVargs: cellarray: varargs input to set title text object
%}

function setSuperAxTitle(titlestring, varargin)

vertPos = .99;
horzPos = .5;
titleVargs = {'FontWeight','bold','FontName','Arial','horizontalAlignment', 'center', ...
    'FontSize',14};
if ~isempty(varargin)
    assign(varargin{:})
end

sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', gcf);
iStitle = text(horzPos, vertPos, {titlestring}, 'Parent', sprax, 'Units', 'normalized');
set(iStitle,titleVargs{:});
h = get(gcf,'Children');

% put the super axes at the bottom of the stack to gain access to
% underlying subplots for direct matlab-figure interaction
set(gcf,'Children',flip(h));
