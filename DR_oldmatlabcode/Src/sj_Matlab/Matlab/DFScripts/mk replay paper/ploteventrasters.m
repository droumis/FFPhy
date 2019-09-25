% plot the fields for a decoding event
%event = 38;

% Fig 2, first raster
event = 55
i = [2 1 2]

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
	tmpspiketimes = decodefilter(i(1)).output{i(2)}(i(3)).eventdata(event).spiketimes(find(decodefilter(i(1)).output{i(2)}(i(3)).eventdata(event).cellindex == matches(trainingcell)));
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

rastertime = mean(decodefilter(i(1)).output{i(2)}(i(3)).eventdata(event).spiketimes);

clist = decodefilter(i(1)).output{i(2)}(i(3)).index(eventcellsactive,:);
%plotfilterfields2d(trainingfilter, i, clist, [2 8]);
adir = decodefilter(i(1)).animal{2};
aname = decodefilter(i(1)).animal{3};
plotrasters(adir, aname, clist, [3913.5 3914.5]);


% Fig 2, second raster
figure
event = 38 

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
	tmpspiketimes = decodefilter(i(1)).output{i(2)}(i(3)).eventdata(event).spiketimes(find(decodefilter(i(1)).output{i(2)}(i(3)).eventdata(event).cellindex == matches(trainingcell)));
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

rastertime = mean(decodefilter(i(1)).output{i(2)}(i(3)).eventdata(event).spiketimes);

clist = decodefilter(i(1)).output{i(2)}(i(3)).index(eventcellsactive,:);
%plotfilterfields2d(trainingfilter, i, clist, [2 8]);
adir = decodefilter(i(1)).animal{2};
aname = decodefilter(i(1)).animal{3};
plotrasters(adir, aname, clist, [3775.8 3776.5]);


% Fig 4, first raster
i = [3 1 7]
event = 11
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
	tmpspiketimes = decodefilter(i(1)).output{i(2)}(i(3)).eventdata(event).spiketimes(find(decodefilter(i(1)).output{i(2)}(i(3)).eventdata(event).cellindex == matches(trainingcell)));
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

rastertime = mean(decodefilter(i(1)).output{i(2)}(i(3)).eventdata(event).spiketimes);

tmpclist = decodefilter(i(1)).output{i(2)}(i(3)).index(eventcellsactive,:);
% we have to switch some around to get the order right for this outside to
% outside decode
clist = [tmpclist [8 7 9 6 5 4 10 11 1 2 3]'];
clist = sortrows(clist, 5);
clist = clist(:,1:4)
%plotfilterfields2d(trainingfilter, i, clist, [2 8]);
adir = decodefilter(i(1)).animal{2};
aname = decodefilter(i(1)).animal{3};
plotrasters(adir, aname, clist, [7412.6 7413.6]);


% Fig 4, second replay event
i = [2 2 2]
event = 106

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

clist = decodefilter(i(1)).output{1}(i(3)).index(eventcellsactive,:);
%plotfilterfields2d(trainingfilter, i, clist, [2 8]);
adir = decodefilter(i(1)).animal{2};
aname = decodefilter(i(1)).animal{3};
plotrasters(adir, aname, clist, [7789 7790]);

