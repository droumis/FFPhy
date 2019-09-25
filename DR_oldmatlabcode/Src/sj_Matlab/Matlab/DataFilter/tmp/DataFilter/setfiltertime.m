function f = setfiltertime(f, timefilter)
% f = setfiltertime(f, timefilter)
% Sets the 'time' field of the filter.  TIMEFILTER is a cell array, where every cell is a call to a specific filtering
% function. Inside every cell is another cell array, where the first
% two cells are always the name of the function and the search terms for
% the function.  Any additional arguments are taken as options to the
% function, in the form 'optionname', optionvalue', ...
% 
% The out put of the filter function must be  a cell structure, 
% for example spikes{}{}{}{}.#### There must be a time field containing a vector of times.  
% All other fields are searcheable, and must be the same length as the time field.
% Example:
% f = setfiltertime(f,{{'linbehavefilter', '(($traj == 1) | ($traj == 3))','includeStates', 2}, {'2Dbehavefilter', '$velocity > 6'}})

for an = 1:length(f)
    if isempty(f(an).animal)
        error(['You must define an animal for the filter before filtering time'])
    end
    if isempty(f(an).epochs)
        error(['You must define the desired epochs for the filter before filtering time'])
    end
    datadir = f(an).animal{2};
    animalprefix = f(an).animal{3};
    
    %find all unique epochs to analyze for the current animal
    totalepochs = [];
    for e = 1:length(f(an).epochs)
        totalepochs = [totalepochs; f(an).epochs{e}];
    end
    totalepochs = unique(totalepochs, 'rows');
    
    %call each filter function using the totalepochs list
    filterresults = [];
    for i = 1:length(timefilter)
        if (length(timefilter{i}) == 2)
            dataStruct = feval(timefilter{i}{1}, datadir, animalprefix, totalepochs);
        elseif (length(timefilter{i}) >2)
            dataStruct = feval(timefilter{i}{1}, datadir, animalprefix, totalepochs, timefilter{i}{3:end});
        else
            error('Each time filter call must have at least two inputs: {function, searchstring, option1, optionvalue1, ...}');
        end
        filterresults{i} = evaluatefilter2(dataStruct, timefilter{i}{2}); %evaluate the filter string for every structure and return results in a cell array
    end
    
    %save the results of each filter sequence as the start and end times
    %for each exclusion period
    for i = 1:length(f(an).epochs)
        if isempty(f(an).epochs{i})
            f(an).excludetime{i} = [];
        end 
        for j = 1:size(f(an).epochs{i},1)
            combinedexclude = [];
            for k = 1:length(timefilter)       
                tmpexcludePeriods = getExcludePeriods(filterresults{k}{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}(:,1),filterresults{k}{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}(:,2));
                combinedexclude = combineExcludePeriods(combinedexclude, tmpexcludePeriods);
            end
            f(an).excludetime{i}{j} = combinedexclude;
        end
    end
end
    
    
