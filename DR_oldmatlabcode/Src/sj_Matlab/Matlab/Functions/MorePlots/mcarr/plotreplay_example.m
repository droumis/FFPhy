% Plot example of replay, CA1 and CA3 raw eeg, and filtered CA1 and CA3
% ripples

%Load decodefilter
load '/data13/mcarr/RipplePaper/decodefilter.mat'
load '/data13/mcarr/RipplePaper/trainingfilter.mat'
load '/data13/mcarr/Bon/boncellinfo.mat'

%Find largest replay event
animal = 1;
event_index = 54;
%Determine which cells fire
matches = rowfind(trainingfilter(animal).output{1}(1).index(:,[1 3 4]),decodefilter(animal).output{1}(1).index(:,[1 3 4])); %find the matching cell indices
startevent = decodefilter(animal).output{1}(1).eventtime(event_index,1);
endevent = decodefilter(animal).output{1}(1).eventtime(event_index,2);

binsize = 0.015;
timebins = startevent:binsize:endevent;

trainingdata = []; spikedata = []; indexlist = []; activecount = 0;
activerates = []; activespiketimes = [];
for trainingcell = 1:length(matches)
    if (matches(trainingcell) > 0) %we have a match
        trainingdata = [trainingdata; trainingfilter(animal).output{1}(1).rates(trainingcell,:)];
        tmpspiketimes = decodefilter(animal).output{1}(1).eventdata(event_index).spiketimes(find(decodefilter(animal).output{1}(1).eventdata(event_index).cellindex == matches(trainingcell)));
        %save all the info for the active cells
        if ~isempty(tmpspiketimes)
            activecount = activecount+1;
        	activespiketimes{activecount} = tmpspiketimes;
            activerates = [activerates; trainingfilter(animal).output{1}(1).rates(trainingcell,:)];
            indexlist = [indexlist; trainingfilter(animal).output{1}(1).index(trainingcell,:)];
        end
    end
end
trainingdata = trainingdata*binsize; %transform rates to expected number of spikes
activerates = activerates*binsize;

distvector = trainingfilter(animal).output{1}(1).dist;
totalsamples = 10000;
spikedata = [];
cellsactive = [];
exponentmatrix = exp(-activerates);
%histogram the spikes for each active cell
for i = 1:length(activespiketimes)
    spikebins = lookup(activespiketimes{i},timebins);  
    spikecount = zeros(1,length(timebins));
    for j = 1:length(spikebins)
        spikecount(spikebins(j)) = spikecount(spikebins(j))+1;
    end
    spikedata = [spikedata; spikecount];
    cellsactive = [cellsactive; (spikecount > 0)];
end

%the decoded data contains spatial probabilities, and is x by t
decodedata = zeros(size(activerates,2),size(spikedata,2));
naninds = find(isnan(activerates(1,:)));
for t = 1:size(spikedata,2) %time  
    Tspikecount = repmat(spikedata(:,t),1,size(activerates,2));
    %calculate P(spikecount|x) for this timebin across all cells and all x 
    spatialprob = prod(((activerates.^Tspikecount)./factorial(Tspikecount)).*exp(-activerates),1)'; 
    spatialprob(naninds) = 0;
    spatialprob = spatialprob/sum(spatialprob);  %normalize across space
    decodedata(:,t) = spatialprob;  
end

trajmapping = [1 1 2 2]; xdata = {[],[]}; ydata = {[],[]};
probdata = {}; probdata = cell(2,1);

%combine the outbound and inbound trajectories
for traj = 1:4;    
    trajindex = find(trainingfilter(animal).output{1}(1).traj == traj);
    xdata{trajmapping(traj)} = trainingfilter(animal).output{1}(1).dist(trajindex);
    ydata{trajmapping(traj)} = stack(ydata{trajmapping(traj)}, decodedata(trajindex,1)');    
    if isempty(probdata{trajmapping(traj)})
        probdata{trajmapping(traj)} = decodedata(trajindex,:);
    else
        probdata{trajmapping(traj)} = probdata{trajmapping(traj)} + decodedata(trajindex,:);
    end
end

%Define color for replay plot
color = zeros(19,3);
color(1:12,1) = [0.2 0.4 0.6 0.8 1 1 1 1 1 1 1 0.2];
color(6:15,2) = [0.2 0.4 0.6 0.8 1 1 1 1 0.8 0.2];
color(11:19,3) = [0.2 0.8 1 1 1 1 0.8 0.6 0.4];
color = color(end:-1:1,:);

%Plot replay plot
traj1 = mean(mean(probdata{1}))>=mean(mean(probdata{2}));
if traj1; traj = 1; else traj = 2; end
figure(1); hold on
for i = size(probdata{traj},2):-1:1
    plot(xdata{traj}+i,i*0.005+probdata{traj}(:,i),'color',color(i,:))
end

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_replay_example.pdf', m, d, y);
print('-dpdf', savestring)

%% Plot rasters, color coded for CA1-CA3, left and right hemispheres
ca1_left = evaluatefilter(cellinfo{3}{2},'isequal($hemisphere,''left'') & isequal($area,''CA1'')');
ca3_left = evaluatefilter(cellinfo{3}{2},'isequal($hemisphere,''left'') & isequal($area,''CA3'')');
ca1_right = evaluatefilter(cellinfo{3}{2},'isequal($hemisphere,''right'') & isequal($area,''CA1'')');
ca3_right = evaluatefilter(cellinfo{3}{2},'isequal($hemisphere,''right'') & isequal($area,''CA3'')');
cell_identity = zeros(size(cellsactive));

%Try sorting
ind = zeros(size(cellsactive,1),1);
for i = 1:size(cellsactive,1)
    ind(i) = find(cellsactive(i,:),1);
end
[y ind] = sort(ind,1,'ascend'); clear y
cellsactive = cellsactive(ind,:);
indexlist = indexlist(ind,:);

for c = 1:size(indexlist,1);
    if rowfind(indexlist(c,[3 4]),ca1_left)
        cell_identity(c,find(cellsactive(c,:))) = -1;
    elseif rowfind(indexlist(c,[3 4]),ca3_left)
        cell_identity(c,find(cellsactive(c,:))) = -3;
    elseif rowfind(indexlist(c,[3 4]),ca1_right)
        cell_identity(c,find(cellsactive(c,:))) = 1;
    elseif rowfind(indexlist(c,[3 4]),ca3_right)
        cell_identity(c,find(cellsactive(c,:))) = 3;
    end
end


for i = 1:size(cellsactive,1)
    if any(cell_identity(i,:)==-1)
        plotraster(activespiketimes{ind(i)},i,0.8,2,'color','b')
    end
    if any(cell_identity(i,:)==1)
        plotraster(activespiketimes{ind(i)},i,0.8,2,'color','c')
    end
    if any(cell_identity(i,:)==-3)
        plotraster(activespiketimes{ind(i)},i,0.8,2,'color','r')
    end
    if any(cell_identity(i,:)==3)
        plotraster(activespiketimes{ind(i)},i,0.8,2,'color','m')
    end   
end
set(gca,'xlim',[timebins(1) timebins(end)],'xtick',timebins(1:3:end),...
    'xticklabel',timebins(1:3:end)-startevent,'ytick',1.4:1:size(cellsactive,1)+0.4,...
    'yticklabel',1:size(cellsactive,1))
box off
xlabel('Time since ripple detection (sec)')
ylabel('Cell count')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_replay_example_raster.pdf', m, d, y);
print('-dpdf', savestring)

%% Plot EEGs for this event

%CA1 left
load '/data13/mcarr/Bon/EEGnonreference/boneeg03-2-13.mat'
etimes = geteegtimes(eeg{3}{2}{13});
eind = lookup([startevent endevent],etimes);
ca1_left = eeg{3}{2}{13}.data(eind(1) - 500:eind(2)+500);
eegtime = etimes(eind(1)-500:eind(2)+500);
starttime = etimes(eind(1));

load '/data13/mcarr/Bon/EEG/bonripple03-2-13.mat'
ca1_left_ripple = ripple{3}{2}{13}.data(eind(1) - 500:eind(2)+500,1);

load '/data13/mcarr/Bon/EEGnonreference/bonlowgamma03-2-13.mat'
ca1_left_gamma = lowgamma{3}{2}{13}.data(eind(1) - 500:eind(2)+500,1);

%CA1 right
load '/data13/mcarr/Bon/EEGnonreference/boneeg03-2-03.mat'
etimes = geteegtimes(eeg{3}{2}{3});
eind = lookup(eegtime,etimes);
ca1_right = eeg{3}{2}{3}.data(eind);

load '/data13/mcarr/Bon/EEG/bonripple03-2-03.mat'
ca1_right_ripple = ripple{3}{2}{3}.data(eind);

load '/data13/mcarr/Bon/EEGnonreference/bonlowgamma03-2-03.mat'
ca1_right_gamma = lowgamma{3}{2}{3}.data(eind);

%CA3 left
load '/data13/mcarr/Bon/EEGnonreference/boneeg03-2-18.mat'
etimes = geteegtimes(eeg{3}{2}{18});
eind = lookup(eegtime,etimes);
ca3_left = eeg{3}{2}{18}.data(eind);

load '/data13/mcarr/Bon/EEG/bonripple03-2-18.mat'
ca3_left_ripple = ripple{3}{2}{18}.data(eind);

load '/data13/mcarr/Bon/EEGnonreference/bonlowgamma03-2-18.mat'
ca3_left_gamma = lowgamma{3}{2}{18}.data(eind);

%CA3 right
load '/data13/mcarr/Bon/EEGnonreference/boneeg03-2-01.mat'
etimes = geteegtimes(eeg{3}{2}{1});
eind = lookup(eegtime,etimes);
ca3_right = eeg{3}{2}{1}.data(eind);

load '/data13/mcarr/Bon/EEG/bonripple03-2-01.mat'
ca3_right_ripple = ripple{3}{2}{1}.data(eind);

load '/data13/mcarr/Bon/EEGnonreference/bonlowgamma03-2-01.mat'
ca3_right_gamma = lowgamma{3}{2}{1}.data(eind);

eegtime = eegtime - starttime;
clear eeg etimes eind ripple lowgamma

%Plot the four eeg traces
figure
plot(eegtime,ca1_left,'b',eegtime,1000+ca1_right,'c',eegtime,2000+ca3_left,'r',eegtime,3000+ca3_right,'m')
set(gca,'xlim',[-0.4 0.7],'ylim',[-800 5000])
xlabel('Time since ripple detection (sec)')
legend([{'CA1 left'},{'CA1 right'},{'CA3 left'},{'CA3 right'}])

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_replay_example_eeg.pdf', m, d, y);
print('-dpdf', savestring)

%Plot the four ripple traces at large time scale and zoomed in
figure
plot(eegtime,ca1_left_ripple,'b',eegtime,200+ca1_right_ripple,'c',...
    eegtime,400+ca3_left_ripple,'r',eegtime,600+ca3_right_ripple,'m')
set(gca,'xlim',[-0.4 0.7],'ylim',[-150 750])
xlabel('Time since ripple detection (sec)')
legend([{'CA1 left'},{'CA1 right'},{'CA3 left'},{'CA3 right'}])
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_replay_example_ripple_trace.pdf', m, d, y);
print('-dpdf', savestring)

figure
timeind = lookup([0 0.2],eegtime);
plot(eegtime(timeind(1):timeind(2)),ca1_left_ripple(timeind(1):timeind(2)),'b',...
    eegtime(timeind(1):timeind(2)),200+ca1_right_ripple(timeind(1):timeind(2)),'c',...
    eegtime(timeind(1):timeind(2)),400+ca3_left_ripple(timeind(1):timeind(2)),'r',...
    eegtime(timeind(1):timeind(2)),600+ca3_right_ripple(timeind(1):timeind(2)),'m')
set(gca,'xlim',[-0.05 0.25],'ylim',[-150 750])
xlabel('Time since ripple detection (sec)')
legend([{'CA1 left'},{'CA1 right'},{'CA3 left'},{'CA3 right'}])
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_replay_example_ripple_trace_zoomedin.pdf', m, d, y);
print('-dpdf', savestring)

%Plot the four gamma traces
figure
plot(eegtime,ca1_left_gamma,'b',eegtime,800+ca1_right_gamma,'c',eegtime,1600+ca3_left_gamma,'r',eegtime,2400+ca3_right_gamma,'m')
set(gca,'xlim',[-0.4 0.7],'ylim',[-500 3000])
xlabel('Time since ripple detection (sec)')
legend([{'CA1 left'},{'CA1 right'},{'CA3 left'},{'CA3 right'}])

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_replay_example_gamma.pdf', m, d, y);
print('-dpdf', savestring)
