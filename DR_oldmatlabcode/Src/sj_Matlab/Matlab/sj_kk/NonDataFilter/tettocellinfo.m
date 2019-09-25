function tettocellinfo(animdirect,fileprefix,tetrodes,field)

% Transfers a given fields value from tetrode to cell, kk 5.3.13

    load([animdirect, fileprefix, 'tetinfo']);
    load([animdirect, fileprefix,'cellinfo']);
    o = cellfetch(cellinfo,'numspikes');
    targetcells = o.index(find(ismember(o.index(:,3),tetrodes)),:);
    
for i = 1:size(targetcells,1)
cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.(field)= ...
    tetinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}.(field);
end

save([animdirect, fileprefix,'cellinfo'], 'cellinfo')