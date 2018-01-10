
function [out] = getNTinfo(index, animals, ian, varargin)
colorSet = 'DR1';
ntindcolumn = 3;
if ~isempty(varargin)
    assign(varargin{:})
end
animalinfo = animaldef(lower(animals{ian}));
animalID = animalinfo{1,3}; %use anim prefix for name
FFanimdir =  sprintf('%s',animalinfo{1,2});
%% ---- loadtetinfostruct ----
load([FFanimdir, animalID, 'tetinfo']);
tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
%         for iday = days
%             ntrodesIndices = ixpc.index{ianimal}{iday};
ntrodesIndices = index;
nNTrodes = size(unique(ntrodesIndices(:,ntindcolumn),'rows','stable'),1);
ntrodes = unique(ntrodesIndices(:,ntindcolumn), 'rows', 'stable');
if nNTrodes ~= length(ntrodes);
    error('data doesnt match indices')
end
%% ---- reorder the LFP traces by suparea|area|subarea tags (in that priority) ----
if length(ntindcolumn) > 1
    for intr = 1:length(ntindcolumn);
        icol = ntindcolumn(intr);
        [~, tagIndMap(:,icol)] = ismember(ntrodes(:,icol),tetinfoAll.index(:,3), 'rows');
        ntrodeTags(:,icol) = tetinfoAll.values(tagIndMap(:,icol));
    end
    
else
    [~, tagIndMap] = ismember(ntrodes,tetinfoAll.index(:,3), 'rows');
    ntrodeTags = tetinfoAll.values(tagIndMap);
end

try
    numsumSupAreas = cellfun(@(x) sum(uint16(x.suparea)), ntrodeTags, 'UniformOutput', false);
    numsumAreas = cellfun(@(x) sum(uint16(x.area)), ntrodeTags, 'UniformOutput', false);
    numsumSubAreas = cellfun(@(x) sum(uint16(x.subarea)), ntrodeTags, 'UniformOutput', false);
    strSupAreas = cellfun(@(x) x.suparea, ntrodeTags, 'UniformOutput', false);
    strAreas = cellfun(@(x) x.area, ntrodeTags, 'UniformOutput', false);
    strSubAreas = cellfun(@(x) x.subarea, ntrodeTags, 'UniformOutput', false);
catch
    error('all ntrodes need to have a suparea, subarea, and area tag, even if blank')
end
icolors = colorPicker(colorSet, strAreas, strSubAreas);
numsumallareatags = cell2mat([numsumSupAreas numsumAreas numsumSubAreas]);
[numsumallSort, numsumSortInds] = sortrows(numsumallareatags);%,[-1 -2 -3]); % -Col to 'descend'
icolors = icolors(numsumSortInds,:);
iout = whos;
outvars = {iout(:).name};
eval(['outvals = {', strjoin(outvars, ','), '};']);
% outvars = {'ntrodes', 'strSupAreas', 'strAreas', 'strSubAreas', 'numsumSupAreas', 'numsumAreas', 'numsumSubAreas', 'icolors', 'ntrodeTags'};
% outvals = {ntrodes, strSupAreas, strAreas, strSubAreas, numsumSupAreas, numsumAreas, numsumSubAreas, icolors, ntrodeTags};
out = cell2struct(outvals', outvars,1);

end