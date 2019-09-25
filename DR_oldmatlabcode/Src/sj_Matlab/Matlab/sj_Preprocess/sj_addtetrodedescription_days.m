function sj_addtetrodedescription_days(animdirect,fileprefix,tetrodes,descrip, days)
% Similar to sj_addtetrodelocation except instead of location ($area), adds
% a description ($descrip)

% Add only for given days

% This function adds tetrode location to both cellinfo and tetinfo. If
% tetinfo is not defined, this function will add location only to cellinfo.
% Note that this assumes each tetrode was in one and only one location
% throughout the recordings.

try
    load([animdirect, fileprefix,'cellinfo']);
    o = cellfetch(cellinfo,'numspikes');
    targetcells = o.index(find(ismember(o.index(:,1),days) & ismember(o.index(:,3),tetrodes)),:);
    for i = 1:size(targetcells,1)
        cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.descrip = descrip;
        save([animdirect, fileprefix,'cellinfo'], 'cellinfo')
    end
catch
    warning('cellinfo structure does not exist')
end



try
    load([animdirect, fileprefix, 'tetinfo']);
    a = cellfetch(tetinfo, 'depth');
    targettets = a.index(find(ismember(a.index(:,1),days) & ismember(a.index(:,3),tetrodes)),:);
    for i =1:size(targettets,1)
        tetinfo{targettets(i,1)}{targettets(i,2)}{targettets(i,3)}.descrip = descrip;
    end
    save([animdirect, fileprefix,'tetinfo'], 'tetinfo');
catch
    warning('tetinfo structure does not exist')
end


 