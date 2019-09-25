function [ncells]= resetCellList
% function [ncells]= resetCellList

global fmaux select
load(fmaux.select);

if isempty(select)
	fmaux.nCells= 0;
else
	fmaux.nCells=size(select.cellnum,1);
end

if isfield(fmaux,'startCell') & ~isempty(fmaux.startCell)
    fmaux.currentCell= fmaux.startCell-1;
else
    fmaux.currentCell= 0;
end

ncells= fmaux.nCells-fmaux.currentCell;
