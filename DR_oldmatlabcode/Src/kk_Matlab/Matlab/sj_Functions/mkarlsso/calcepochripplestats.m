function out = calcepochripplestats(index, trainingindex, trainingfilter, decodefilter)

animalnum = index(1);
epochnum = index(2);
trajmapping = [1 1 2 2];
binsize = .015; %default temporal bin
out = [];

for eventindex = 1:length(decodefilter(animalnum).output{1}(epochnum).eventdata)
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
    matches = rowfind(trainingfilter(animalnum).output{trainingindex}(epochnum).index(:,[1 3 4]),decodefilter(animalnum).output{1}(epochnum).index(:,[1 3 4])); %find the matching cell indices
    startevent = decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,1);
    endevent = decodefilter(animalnum).output{1}(epochnum).eventtime(eventindex,2);
    if ((endevent-startevent) < 2)
        timebins = startevent:binsize:endevent;
        eventcellsactive = [];
        activecount = 0;
        for trainingcell = 1:length(matches)
            if (matches(trainingcell) > 0) %we have a match
                indexlist = [indexlist; trainingfilter(animalnum).output{trainingindex}(epochnum).index(trainingcell,:)];
                trainingdata = [trainingdata; trainingfilter(animalnum).output{trainingindex}(epochnum).rates(trainingcell,:)];

                tmpspiketimes = decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).spiketimes(find(decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).cellindex == matches(trainingcell)));
                %save all the info for the active cells
                if ~isempty(tmpspiketimes)
                    activecount = activecount+1;
                    activespiketimes{activecount} = tmpspiketimes;
                    activerates = [activerates; trainingfilter(animalnum).output{trainingindex}(epochnum).rates(trainingcell,:)];
                end

            end
        end
        trainingdata = trainingdata*binsize; %transform rates to expected number of spikes
        activerates = activerates*binsize;          
    
        %out = [out; [decodefilter(animalnum).output{1}(epochnum).eventimmobiletime(eventindex) length(decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).cellindex)]];  %decodefilter(animalnum).output{1}(epochnum).std(eventindex)
        out = [out; [decodefilter(animalnum).output{1}(epochnum).eventimmobiletime(eventindex) length(activespiketimes)]];  %decodefilter(animalnum).output{1}(epochnum).std(eventindex)
    end
end