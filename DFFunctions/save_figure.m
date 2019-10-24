
% save matlab figs to common academic and data analysis formats
% save_figure(figdirectory, filename, varargin)
% varg - subdir: string (default = '')
% varg - savefigas: cellarray of file formats (default={'png'}) [png, mfig, eps, svg, pdf]

% 2019 Demetris Roumis

function save_figure(figdirectory, filename, varargin)

set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
subdir = '';
savefigas = {'png'};
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
strsave = strrep(filename,' ', '_'); % '_' in place of whitespace for filename
if ~isa(savefigas,'cell')
    savefigas = {savefigas};
end
    for s = 1:length(savefigas)
        switch savefigas{s}
            case 'png'
                print(sprintf('%s/%s',figdir, strsave),'-dpng', '-r0')
                fprintf('saved %s/%s.png\n', figdir, strsave)
            case 'mfig'
                savefig(gcf, sprintf('%s/%s',figdir, strsave),'compact') % compact: 2014b or higher
                fprintf('saved %s/%s.fig\n', figdir, strsave)
            case 'eps'
                % using a loose bounding box and a TIFF preview
                print(sprintf('%s/%s',figdir, strsave),'-depsc','-tiff','-loose')
                fprintf('saved %s/%s.eps\n', figdir, strsave)
            case 'svg'
                print(sprintf('%s/%s',figdir, strsave),'-dsvg')
                fprintf('saved %s/%s.svg\n', figdir, strsave)
            case 'pdf'
                print(sprintf('%s/%s',figdir, strsave),'-dpdf', '-bestfit')
                fprintf('saved %s/%s.pdf\n', figdir, strsave)
        end
    end
end