function addtetrodehemisphere(animdirect,fileprefix,tetrodes,location)
%This function defines a hemisphere for each tetrode. This information is
%added to both tetinfo and cell info.

%Options are 'right' and 'left'

load([animdirect, fileprefix, 'tetinfo']);
load([animdirect, fileprefix, 'cellinfo']);

a = cellfetch(cellinfo, 'depth');
targettets = a.index(ismember(a.index(:,3),tetrodes),:);

for i =1:size(targettets,1)
    cellinfo{targettets(i,1)}{targettets(i,2)}{targettets(i,3)}{targettets(i,4)}.hemisphere = location;
end
% 
% a = cellfetch(tetinfo,'depth');
% targettets = a.index(find(ismember(a.index(:,3),tetrodes)),:);
% 
% for i = 1:size(targettets,1)
%     tetinfo{targettets(i,1)}{targettets(i,2)}{targettets(i,3)}.hemisphere = location;
% end

%save([animdirect, fileprefix,'tetinfo'], 'tetinfo');
save([animdirect, fileprefix,'cellinfo'],'cellinfo');

end

 