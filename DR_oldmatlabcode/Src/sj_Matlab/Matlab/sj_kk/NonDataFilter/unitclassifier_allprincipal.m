function unitclassifier_allprincipal(animdirect,fileprefix,days,tetrodes)

% This function assigns all units in a cellinfostruct to be principal cells from the selected
% tetrodes.

load([animdirect, fileprefix,'cellinfo']);
o = cellfetch(cellinfo,'numspikes');
targetcells = o.index(find(ismember(o.index(:,1),days) & ismember(o.index(:,3),tetrodes)),:);
for i = 1:size(targetcells,1)
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.type = 'principal';
end

save([animdirect, fileprefix,'cellinfo'], 'cellinfo')

end


 