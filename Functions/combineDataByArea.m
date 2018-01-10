function [areas2procStr,areaDataMean, areas2procNum]  = combineDataByArea(index, animals, ian, dataInput, varargin)
ntindcolumn = 3;
if ~isempty(varargin)
    assign(varargin{:});
end
[nt] = getNTinfo(index, animals, ian, 'ntindcolumn', ntindcolumn);
[~, unqindsStr] = unique(cell2mat(nt.strAreas), 'rows', 'stable');
[~, unqinds] = unique(cell2mat(nt.numsumAreas), 'rows', 'stable');
if any(unqindsStr ~= unqinds)
    error('int16 sum of area strings needs to be unqique..')
end
areas2procWrefStr = nt.strAreas(unqinds, :);
areas2procWrefNum = nt.numsumAreas(unqinds, :);
% areas2proc2 = cellstr(areas2proc);
useareas2procindsWref = cellfun(@(x) isempty(x), strfind(areas2procWrefStr, 'ref'), 'un', 1);
useareas2procinds = logical(prod(useareas2procindsWref, 2));
areas2procStr = areas2procWrefStr(useareas2procinds, :);
areas2procNum = cell2mat(areas2procWrefNum(useareas2procinds, :));
% areas2proc = areas2procWref(logical(arrayfun(@(x) prod([useareas2procindsWref(x,:)]), [1:size(useareas2procindsWref,1)],'un',1)'),:);
% areas2proc = areas2proc(logical(sum(bsxfun(@times, useareas2procinds,fliplr(useareas2procinds)),2)), :);
% areas2proc = areas2proc(find(strfind(areas2proc, 'ref'))); %ignore the reference ntrode
% areas{ian} = areas2procStr;
% areasNums{ian} = areas2procNum;
% a = cellstr(nt.strAreas);
% a = cellfun(@(x) strjoin(x), nt.strAreas, 'un', 0);
nt.numsumAreas = cell2mat(nt.numsumAreas);
for iar = 1:size(areas2procNum,1)
    areaID = areas2procNum(iar,:);
    areainds = ismember(nt.numsumAreas,areaID, 'rows');
    areaoutputtmp = dataInput(areainds);
    for idT = 1:length(areaoutputtmp{1})
        idTdata = arrayfun(@(x) areaoutputtmp{x}{idT},[1:length(areaoutputtmp)], 'un', 0);
        areaDataMean{iar}{idT} = mean(cat(4,idTdata{:}), 4);
    end
    %     b = reshape(cat(2,areaoutputtmp{:}),size(areaoutputtmp,1),size(areaoutputtmp{1},2));
    %     c = cat(4,b{:});
    %     areaDataMean{iar} = mean(cat(4,areaoutputtmp{:}), 4);
    %     [areaDataMean,  = getAreaMean(ixpc, areainds, ian, iar, datatype);
    %     eval(sprintf,'ixpc.area%smean
end
end