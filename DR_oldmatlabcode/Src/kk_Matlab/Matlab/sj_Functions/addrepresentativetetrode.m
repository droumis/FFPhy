function addrepresentativetetrode(animdirect,fileprefix,tetrodes)
%This function marks the tetrode that is listed as the most representative
%of that layer. Use histology and drive location to pick most
%representative tetrode for given MEC layer.

load([animdirect, fileprefix, 'tetinfo']);

a = cellfetch(tetinfo, 'depth');
targettets = a.index(find(ismember(a.index(:,3),tetrodes)),:);

for i =1:size(targettets,1)
    tetinfo{targettets(i,1)}{targettets(i,2)}{targettets(i,3)}.representative = 1;
end

save([animdirect, fileprefix,'tetinfo'], 'tetinfo');

end

 