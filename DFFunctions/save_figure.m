function save_figure(figdirectory, filename, varargin) 
set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
subdir = '';
if ~isempty(varargin)
    assign(varargin{:});
end
if ~isempty(subdir)
    figdir = sprintf('%s/%s/',figdirectory, subdir);
else
    figdir = sprintf('%s/',figdirectory);
end
if ~isdir(figdir)
    mkdir(figdir);
end
strsave = strrep(filename,' ', '_'); %put '_' character in place of whitespace for filename
print(sprintf('%s/%s',figdir, strsave),'-dpng', '-r0')
fprintf('saved %s/%s.png\n', figdir, strsave)
end