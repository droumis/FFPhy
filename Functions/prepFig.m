function [ifig1, sfrows, sfcols] = prepFig(savefigs, pausefigs, nt, fig, varargin)
backcolor = 'white';
figspecs = '';
nDays = [];
if ~isempty(varargin)
    assign(varargin{:});
end
%if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
if savefigs && ~pausefigs
    ifig1 = figure('Visible','off','units','normalized','position',fig.position);
else
    ifig1 = figure('units','normalized','position',fig.position);
end
set(ifig1,'color',backcolor)

if strcmp(figspecs,'AreasByDays')
    sfrows = size(nt.areas,1);
    sfcols = nDays;
else
try
    sfrows = floor(sqrt(nt.nNTrodes));
    sfcols = ceil(sqrt(nt.nNTrodes));
catch
    sfrows = 1;
    sfcols = length(nt.areas{nt.ian});
end
end

end