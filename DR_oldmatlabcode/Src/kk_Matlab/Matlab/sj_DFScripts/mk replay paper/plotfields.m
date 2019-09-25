% plot the fields for a decoding event
%event = 38;
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

clist = decodefilter(i(1)).output{i(2)}(i(3)).index(eventcellsactive,:);
plotfilterfields2d(trainingfilter, i, clist, [2 8]);




% figure 4 replay #1
event = 11
i = [3 1 7]

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

clist = decodefilter(i(1)).output{i(2)}(i(3)).index(eventcellsactive,:);
plotfilterfields2d(trainingfilter, i, clist, [2 6]);


% figure 4 replay #2

event = 106
i = [2 2 2]

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

clist = decodefilter(i(1)).output{1}(i(3)).index(eventcellsactive,:);
plotfilterfields2d(trainingfilter, i, clist, [5 3]);


% plot all of the fields from Frank dataset 8
anum = 1;
epochnum = 7;
i1 = placefieldfilter(anum).output{1}(epochnum).index;   
i2 = placefieldfilter(anum).output{2}(epochnum).index;   
t1 = trainingfilter(anum).output{1}(epochnum).index;   
t2 = trainingfilter(anum).output{2}(epochnum).index;   
corresp = rowfind(i1(:,[1 3 4]), i2(:,[1 3 4]));
valid = find(corresp);
clist1 = i1(valid,:);
clist2 = i2(corresp(valid),:);
% find the cells that have place fields in either t1 or t2
valid1 = rowfind(clist1, t1);
valid2 = rowfind(clist2, t2);
clist1 = clist1(find(valid1 | valid2), :);
clist2 = clist2(find(valid1 | valid2), :);
c = [24 19 25 14 28:33];
clist1 = clist1(c,:);
clist2 = clist2(c,:);
i = [anum 1 epochnum];
orient tall
plotfilterfields2d(placefieldfilter, i, clist1, [22 length(clist1)]);
figure
orient tall
i = [anum 2 epochnum];
plotfilterfields2d(placefieldfilter, i, clist2, [22 length(clist2)]);

