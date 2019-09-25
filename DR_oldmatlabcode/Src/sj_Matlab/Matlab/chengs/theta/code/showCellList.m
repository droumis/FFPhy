function cl= showCellList(selectid)

if nargin<1
    global fmaux
    load(fmaux.select);
    cl= select.cellnum;
else
    load(['select-' selectid]);
    cl= select.cellnum;
    select
end
