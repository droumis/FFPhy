function selectCells(region, celltype, cleanup)
%function selectCells(selectOut, cleanup)
%
%  Select cells (by selectOut) from cells in selectIn.
%  region= {'CA1', 'EC', 'Ecto', 'Sub'}
%  celltype= {'PE', 'FS'}
%   cleanup= 1: remove empty cell entry 
%   

% parameters
if nargin<3; cleanup= 1; end
selectOut= [region celltype];

global fmaux linbehav linmeanvel behavdata nbinx select task tetloc spikes

loadFile(fmaux.select);
loadVar(fmaux.datadir, 'tetloc', 0, fmaux.prefix, 1);
loadVar(fmaux.datadir, 'task', 0, fmaux.prefix, 1);
if ~isfield(select, 'x');
    select.x= {};
end

[d,e,t,c]= startCellList;
while ~isempty(d)
    n= fmaux.currentCell;
    loadVar(fmaux.datadir, 'spikes', d, fmaux.prefix, 1);
%    disp([tetloc{d}{t}.region '-' spikes{d}{e}{t}{c}.celltype]);

    if  n>length(select.a) | isempty(select.a{n}) | ...
        isempty(tetloc{d}) | isempty(tetloc{d}{t}) | ...
        ~strcmp(tetloc{d}{t}.region, region) | ...  % only CA1 cells count
        ~strcmp(spikes{d}{e}{t}{c}.celltype, celltype) % only excitatory cells
        % reject traj
        select.a{n}= {};
        select.x{n}= {};
    end
    [d,e,t,c]= getNextCell;
end

if(cleanup)
    disp('**** clean up ****')
    select= cleanSelect(select);
end

save([fmaux.data2dir '/select-' selectOut], 'select');

