i = [1 1 1];

for event = 1:length(decodefilter(i(1)).output{i(2)}(i(3)).eventtraj)
% figure out the order of the first spikes from each cell

matches = rowfind(trainingfilter(i(1)).output{i(2)}(i(3)).index(:,[1 3 4]),decodefilter(i(1)).output{1}(i(3)).index(:,[1 3 4]));
traininddata = [];
eventcellsactive = [];
eventcellpeaks = [];
trainingdata = [];
spikedata = [];
decodedata = [];
indexlist = [];
for trainingcell = 1:length(matches)
    if (matches(trainingcell) > 0) %we have a match
        indexlist = [indexlist; trainingfilter(i(1)).output{i(2)}(i(3)).index(trainingcell,:)];
        trainingdata = [trainingdata; trainingfilter(i(1)).output{i(2)}(i(3)).rates(trainingcell,:)];      
        [tmppeak, tmppeakind] = max(trainingfilter(i(1)).output{i(2)}(i(3)).rates(trainingcell,:));
        tmppeakdist = trainingfilter(i(1)).output{i(2)}(i(3)).dist(tmppeakind);
	tmpspiketimes = decodefilter(i(1)).output{1}(i(3)).eventdata(event).spiketimes(find(decodefilter(i(1)).output{1}(i(3)).eventdata(event).cellindex == matches(trainingcell)));
        %save all the info for the active cells
        if ~isempty(tmpspiketimes)
	    eventcellsactive = [eventcellsactive matches(trainingcell)];
	    eventcellpeaks = [eventcellpeaks tmppeakdist];
	end
    end
end
eventcells = [eventcellsactive' eventcellpeaks'];
eventcells = sortrows(eventcells,2);
eventcellsactive = eventcells(:,1);
eventcellpeaks = eventcells(:,2);


rastertime = mean(decodefilter(i(1)).output{1}(i(3)).eventdata(event).spiketimes);
eventtime = decodefilter(i(1)).output{i(2)}(i(3)).eventtime(event,:);
clist = decodefilter(i(1)).output{1}(i(3)).index(eventcellsactive,:);
%plotfilterfields2d(trainingfilter, i, clist, [2 8]);
adir = decodefilter(i(1)).animal{2};
aname = decodefilter(i(1)).animal{3};
plotrasters(adir, aname, clist, [eventtime(1) eventtime(2)]);
pause
end