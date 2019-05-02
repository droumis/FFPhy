function save_figure(figdirectory, filenamesave, sprtit) 

figdir = sprintf('%s/%s/',figdirectory, filenamesave);

if ~isdir(figdir);
    mkdir(figdir);
end
if ~isdir(figdir);
    mkdir(figdir);
end
sprtitsave = strrep(sprtit,' ', '_'); %put '_' character in place of whitespace for filename
set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
currfigfile = sprintf('%s/%s',figdir, sprtitsave);
print(currfigfile,'-dpng', '-r0')
fprintf('plot %s saved \n', sprtit)
end