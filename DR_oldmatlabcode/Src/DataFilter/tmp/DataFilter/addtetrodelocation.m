function addtetrodelocation(animdirect,fileprefix,tetrodes,location)

load([animdirect, fileprefix,'cellinfo']);
o = cellfetch(cellinfo,'numspikes');
targetcells = o.index(find(ismember(o.index(:,3),tetrodes)),:);
for i = 1:size(targetcells,1)
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.area = location;
end
save([animdirect, fileprefix,'cellinfo'], 'cellinfo');
 