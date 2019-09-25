function f = setfiltertime(f, timefilter)
% f = setfiltertime(f, timefilter)
% Sets the 'time' field of the filter.  TIMEFILTER is a cell array, where every cell is a call to a specific filtering
% function. Inside every cell is another cell array, where 
% the first cell is the name of a function.
% the second cell is either a string 
%	1. specifitying the search terms for that function if the function 
%	returns a cell structure (see below) or 
%	2. 'excludetimelist' which indicates that the function will a cell
%	array where each element c{dataset}{epoch} is a list of exclude times 
%	in an Nx2 double matrix where each row is the start and end time for 
%	the exclude list.
%
% cells 3-n are options to the function 
% 
% If 'excludetimelist' is not used, the output of the filter function must be  a cell structure, 
% for example spikes{}{}{}{}.#### There must be a time field containing a vector of times.  
% All other fields are searcheable, and must be the same length as the time field.
% Example:
% f = setfiltertime(f,{{'linbehavefilter', '(($traj == 1) | ($traj == 3))','includeStates', 2}, {'2Dbehavefilter', '$velocity > 6'}})
%
% if exclude time list is used then, as mentioned above, the function must
% output an Nx2 double matrix with exclude times

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
    loaddays = unique(totalepochs(:,1));
    if ~isempty(totalepochs)
        %call each filter function using the totalepochs list
        filterresults = [];
        for i = 1:length(timefilter) %one filter for each experimental condition
            if ((length(timefilter{i}) >= 2) && (mod(length(timefilter{i}),2) == 0))          
                for j = 1:2:length(timefilter{i}) %for combined filters
                    if (strfind(lower(timefilter{i}{j}),'<function>')) %this is a function call
                        fparam = parseFunctionCallString(timefilter{i}{j});
                        dataStruct = feval(fparam.funcName, datadir, animalprefix, totalepochs, fparam.funcVarargin{:});
                        filterresults{i}{(j+1)/2} = evaluatefilter2(dataStruct, timefilter{i}{j+1}); %evaluate the filter string for every structure and return results in a cell array
                    elseif (strfind(lower(timefilter{i}{j}),'<excludetimelist>')) %the input is already a list for each epoch (a double cell array)
                        filterresults{i}{(j+1)/2} = timefilter{i}{j+1};
                    else %assume the input is the name of a variable to search                                          
                        dataStruct = loaddatastruct(datadir, animalprefix, timefilter{i}{j}, loaddays);
                        filterresults{i}{(j+1)/2} = evaluatefilter2(dataStruct, timefilter{i}{j+1}); %evaluate the filter string for every structure and return results in a cell array
                    end
                end
            else
                error('Each time filter call must have inputs: {datatype1, searchstring1, datatype2, searchstring2, ... }');
            end
                                   
        end
    end
    
    %save the results of each filter sequence as the start and end times
    %for each exclusion period
    for i = 1:length(f(an).epochs)
        if isempty(f(an).epochs{i})
            f(an).excludetime{i} = [];
        end
        for j = 1:size(f(an).epochs{i},1)
            combinedexclude = [];
            
            if (length(filterresults) == 1)
                for k = 1:length(filterresults{1})
                    tmpexcludePeriods = getExcludePeriods(filterresults{1}{k}{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}(:,1),filterresults{1}{k}{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}(:,2));                    
                    
                    combinedexclude = combineExcludePeriods(combinedexclude, tmpexcludePeriods);
                end
            elseif (length(filterresults) == length(f(an).epochs))
                for k = 1:length(filterresults{i})
                    tmpexcludePeriods = getExcludePeriods(filterresults{i}{k}{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}(:,1),filterresults{i}{k}{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}(:,2));                    
                    combinedexclude = combineExcludePeriods(combinedexclude, tmpexcludePeriods);
                end
            else 
                error('The filter must have length 1 or the same length as the number of epochs');
                
            end
            f(an).excludetime{i}{j} = combinedexclude;
        end
    end
    f(an).arguments.excludetime = timefilter;
end

    
