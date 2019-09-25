function addtetrodelocation_depth(animdirect,fileprefix,day,tetrode,depth)

% adds depth to tetinfostruct for a tetrode to all epochs of a day, kk 3.18.13

try
    load([animdirect, fileprefix, 'tetinfo']);
    a = cellfetch(tetinfo, 'depth');
    targetfields = a.index(find(ismember(a.index(:,1),day) & ismember(a.index(:,3),tetrode)),:);
    for i =1:size(targetfields,1)
        tetinfo{targetfields(i,1)}{targetfields(i,2)}{targetfields(i,3)}.depth = {depth};
    end
    save([animdirect, fileprefix,'tetinfo'], 'tetinfo');
catch
    warning('tetinfo structure does not exist')
end


 