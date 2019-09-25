function h = plotrasters(adir, aname, clist, times);
% h = plotrasters(animaldir, animalname, celllist, times);
%     plots, in order, the spikes from the cell list between the given times
%

spikes = loaddatastruct(adir, aname, 'spikes');
scoord = [];
for c = 1:size(clist, 1)
    tmps = spikes{clist(c,1)}{clist(c,2)}{clist(c,3)}{clist(c,4)}.data(:,1);
    sind = find((tmps >= times(1)) & (tmps <= times(2)));
    if (~isempty(sind))
	scoord = [scoord ; tmps(sind) ones(size(sind))*c];
    end
end
scoord(:,2) = scoord(:,2) - 0.4;
plotraster(scoord(:,1), scoord(:,2), .8, []);
set(gca, 'XLim', times);
set(gca, 'YLim', [0 ceil(max(scoord(:,2)))+1]);


