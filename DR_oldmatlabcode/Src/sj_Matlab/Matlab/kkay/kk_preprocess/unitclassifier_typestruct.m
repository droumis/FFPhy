function unitclassifier_typestruct(animdirect,fileprefix,days,tetrodes,typestruct)

% This function assigns all units in a cellinfostruct to be principal cells from the selected
% tetrodes.

load([animdirect, fileprefix,'cellinfo']);
o = cellfetch(cellinfo,'numspikes');

targetcells = o.index(find(ismember(o.index(:,1),days) & ismember(o.index(:,3),tetrodes)),:);

for k=1:length(typestruct)     % iterate through cell types
        
    for i = 1:size(targetcells,1)
        if rowfind(targetcells(i,[1 3 4]),typestruct(k).daytetcell)
            cellinfo{targetcells(i,1)}{targetcells(i,2)}{targetcells(i,3)}{targetcells(i,4)}.type = typestruct(k).type;
            s=sprintf('changed .type of %d %d %d to %s',targetcells(i,[1 3 4]),typestruct(k).type);
            disp(s)
        end
    end
    
end

save([animdirect, fileprefix,'cellinfo'], 'cellinfo')

end


 