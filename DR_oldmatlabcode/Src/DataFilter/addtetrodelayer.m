function addtetrodelayer(animdirect,fileprefix,tetrodes,layer)
% This function adds tetrode location to both cellinfo and tetinfo. If
% tetinfo is not defined, this function will add location only to cellinfo.
% Note that this assumes each tetrode was in one and only one location
% throughout the recordings.


load([animdirect, fileprefix,'cellinfo']);
o = cellfetch(cellinfo,'numspikes');
targetcells = o.index(find(ismember(o.index(:,3),tetrodes)),:);
for i = 1:size(targetcells,1)
    cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.layer = layer;
end

save([animdirect, fileprefix,'cellinfo'], 'cellinfo')

try
    load([animdirect, fileprefix, 'tetinfo']);
    a = cellfetch(tetinfo, 'depth');
    targettets = a.index(find(ismember(a.index(:,3),tetrodes)),:);
    for i =1:size(targettets,1)
        tetinfo{targettets(i,1)}{targettets(i,2)}{targettets(i,3)}.layer = layer;
    end
    save([animdirect, fileprefix,'tetinfo'], 'tetinfo');
catch
    warning('tetinfo structure does not exist')
end


 