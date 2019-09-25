function out = calcepochactivationstats(index, trainingindex, trainingfilter, decodefilter)

animalnum = index(1);
epochnum = index(2);
out = [];

for eventindex = 1:length(decodefilter(animalnum).output{1}(epochnum).eventdata)
    

    %pick out all the matching cells from the training data and the
    %decoding data
   
    matches = rowfind(trainingfilter(animalnum).output{trainingindex}(epochnum).index(:,[1 3 4]),decodefilter(animalnum).output{1}(epochnum).index(:,[1 3 4])); %find the matching cell indices
    tmppeak = decodefilter(animalnum).output{1}(epochnum).peak(eventindex);
    activecount = 0;
    spikecount = 0;
    for trainingcell = 1:length(matches)
        if (matches(trainingcell) > 0) %we have a match
            
            tmpspiketimes = decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).spiketimes(find(decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).cellindex == matches(trainingcell)));
            
            %save all the info for the active cells
            if ~isempty(tmpspiketimes)
                activecount = activecount+1;
                spikecount = spikecount + length(tmpspiketimes);
                
            end

        end
    end
    percentactive = activecount/size(trainingfilter(animalnum).output{trainingindex}(epochnum).index,1);
    
    

    if (activecount > 4)
        out = [out; [activecount percentactive spikecount tmppeak]];
    end
end