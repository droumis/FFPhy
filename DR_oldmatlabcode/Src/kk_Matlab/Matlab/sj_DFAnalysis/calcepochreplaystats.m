function out = calcepochreplaystats(index, trainingindex, trainingfilter, decodefilter)

animal = index(1);
day = index(2);
epoch = index(3);
trajmapping = [1 1 2 2];
binsize = .015; %default temporal bin
out = [];

for eventindex = 1:length(decodefilter(animal).output{day}(epoch).eventdata)
    disp(eventindex)
    trainingdata = [];
    spikedata = [];
    decodedata = [];
    indexlist = [];
    activespiketimes = [];
    activerates = [];

    %pick out all the matching cells from the training data and the
    %decoding data
    %traindata contains linear rates, and is n by x, where n is the
    %number of cells and x is the number of spatial bins
    %spikedata contains spikecounts, and is n by t, where t is the
    %number of temporal bins in the data to be decoded.
    matches = rowfind(trainingfilter(animal).output{trainingindex(1)}(trainingindex(2)).index(:,[1 3 4]),decodefilter(animal).output{day}(epoch).index(:,[1 3 4])); %find the matching cell indices
    startevent = decodefilter(animal).output{day}(epoch).eventtime(eventindex,1);
    endevent = decodefilter(animal).output{day}(epoch).eventtime(eventindex,2);
    if ((endevent-startevent) < 2)
        timebins = startevent:binsize:endevent;
        eventcellsactive = [];
        activecount = 0;
        for trainingcell = 1:length(matches)
            if (matches(trainingcell) > 0) %we have a match
                indexlist = [indexlist; trainingfilter(animal).output{trainingindex(1)}(trainingindex(2)).index(trainingcell,:)];
                trainingdata = [trainingdata; trainingfilter(animal).output{trainingindex(1)}(trainingindex(2)).rates(trainingcell,:)];

                tmpspiketimes = decodefilter(animal).output{day}(epoch).eventdata(eventindex).spiketimes(find(decodefilter(animal).output{day}(epoch).eventdata(eventindex).cellindex == matches(trainingcell)));
                %save all the info for the active cells
                if ~isempty(tmpspiketimes)
                    activecount = activecount+1;
                    activespiketimes{activecount} = tmpspiketimes;
                    activerates = [activerates; trainingfilter(animal).output{trainingindex(1)}(trainingindex(2)).rates(trainingcell,:)];
                end

            end
        end
        trainingdata = trainingdata*binsize; %transform rates to expected number of spikes
        activerates = activerates*binsize;

        if (length(activespiketimes) >= 4)
            %out = [out; calcReplayStats(activespiketimes,activerates,timebins,trainingfilter(animal).output{trainingindex(1)}(trainingindex(2)).dist)];
            %out = [out; [calcReplayStats(activespiketimes,activerates,timebins,trainingfilter(animal).output{trainingindex(1)}(trainingindex(2)).dist) decodefilter(animal).output{day}(epoch).eventimmobiletime(eventindex) length(activespiketimes) ]];  %decodefilter(animal).output{1}(epochnum).std(eventindex)
            %out = [out; [decodefilter(animal).output{day}(epoch).eventimmobiletime(eventindex) length(activespiketimes) ]];  %decodefilter(animal).output{day}(epoch).std(eventindex)
            %out = [out;index(1)];
            out = [out; [calcReplayStats(activespiketimes,activerates,timebins,trainingfilter(animal).output{trainingindex(1)}(trainingindex(2)).dist) length(activespiketimes) eventindex]];
        end
    end
end