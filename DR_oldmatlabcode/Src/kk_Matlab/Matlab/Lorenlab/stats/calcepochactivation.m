function [counts, totals] = calcepochactivation(index, trainingfilter, decodefilter)

animalnum = index(1);
epochnum = index(2);
counts = [];


%find the total number of active cells in the training data
totals(1) = size(trainingfilter(animalnum).output{1}(epochnum).index,1);
totals(2) = size(trainingfilter(animalnum).output{2}(epochnum).index,1);

%count the total cells activated from the set of cells in each training set
for eventindex = 1:length(decodefilter(animalnum).output{1}(epochnum).eventdata)  
    activecount = [0 0];
    for trainingindex = 1:2
        matches = rowfind(trainingfilter(animalnum).output{trainingindex}(epochnum).index(:,[1 3 4]),decodefilter(animalnum).output{1}(epochnum).index(:,[1 3 4])); %find the matching cell indices       
        for trainingcell = 1:length(matches)
            if (matches(trainingcell) > 0) %we have a match

                tmpspiketimes = decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).spiketimes(find(decodefilter(animalnum).output{1}(epochnum).eventdata(eventindex).cellindex == matches(trainingcell)));
                %save all the info for the active cells
                if ~isempty(tmpspiketimes)
                    activecount(trainingindex) = activecount(trainingindex)+1;
                end
            end
        end
    end
    counts = [counts; activecount];
end