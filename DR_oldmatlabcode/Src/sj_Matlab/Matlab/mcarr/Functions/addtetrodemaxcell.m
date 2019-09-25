function addtetrodemaxcell(animdirect,fileprefix,tetrodes)
%This function marks the tetrode that has the most cells of the list of
%tetrodes.

load([animdirect, fileprefix, 'tetinfo']);

a = cellfetch(tetinfo, 'numcells');
targettets = a.index(find(ismember(a.index(:,3),tetrodes)),:);

dind = unique(targettets(:,1));
for d = dind(1):dind(end)
    eind = unique(targettets(targettets(:,1)==d,2));
    for e = eind(1):eind(end)
        tetrodes = targettets(targettets(:,1)==d & targettets(:,2)==e,3);
        ind = find(a.index(:,1)==d & a.index(:,2)==e & ismember(a.index(:,3),tetrodes));
        if max(cell2mat(a.values(ind))) > 0
            maxcells = find(cell2mat(a.values(ind))==max(cell2mat(a.values(ind)))); 
            tetindex = a.index(ind(maxcells(end)),:);
            tetinfo{tetindex(1)}{tetindex(2)}{tetindex(3)}.maxcell = 1;
        end
    end
end

save([animdirect, fileprefix,'tetinfo'], 'tetinfo');

end