

%{
set super figure axis title
%}

function sprtit = setSuperAxTitle(animal, figname, varargin)
addtitle = '';
if ~isempty(varargin)
    assign(varargin{:})
end
sprtit = sprintf('%s %s %s', figname,  animal, addtitle);
sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', gcf);
iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment', 'center', ...
    'FontSize',12);
h = get(gcf,'Children');
set(gcf,'Children',flip(h));
