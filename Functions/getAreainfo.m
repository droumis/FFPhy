
function [out] = getAreainfo(areas, ian, varargin);
colorSet = 'DR1';
if ~isempty(varargin)
    assign(varargin{:})
end
strAreas = areas;
strSubAreas = repmat({'all'}, size(strAreas,1),size(strAreas,2));
icolors = colorPicker(colorSet, strAreas, strSubAreas);
iout = whos;
outvars = {iout(:).name};
eval(['outvals = {', strjoin(outvars, ','), '};']);
out = cell2struct(outvals', outvars,1);
end