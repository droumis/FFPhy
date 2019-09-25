function f = setfiltertrials(f, filterInput)

for an = 1:length(f)
    if isempty(f(an).animal)
        error(['You must define an animal for the filter before filtering the trials'])
    end
    if isempty(f(an).epochs)
        error(['You must define the desired epochs for the filter before filtering the trials'])
    end
    datadir = f(an).animal{2};
    animalprefix = f(an).animal{3};
    trials = loaddatastruct(datadir,animalprefix,'trials');
    
    %find all unique epochs to analyze for the current animal
    totalepochs = [];
    for e = 1:length(f(an).epochs)
        totalepochs = [totalepochs; f(an).epochs{e}];
    end
    totalepochs = unique(totalepochs, 'rows');
    
    if ~isempty(totalepochs)
        if (iscell(filterInput) && (length(filterInput) == length(f(an).epochs))) %if there are multiple filters in a cell array, create multiple epoch groups
            for j = 1:length(filterInput)
                if isstr(filterInput{j})
                    filterresult = evaluatefilter2(trials,filterInput{j});
                    
                    for k = 1:size(f(an).epochs{j},1)
                        f(an).trials{j}{k} = filterresult{f(an).epochs{j}(k,1)}{f(an).epochs{j}(k,2)};
                        
                    end
                else
                    error('Each cell in filterInput must contain a string');
                end
            end
        elseif (iscell(filterInput) && (length(filterInput) == 1)) %if there is only one filter string, use the same filter for all groups
            filterresult = evaluatefilter2(trials,filterInput{1});
            for j = 1:length(f(an).epochs)
                for k = 1:size(f(an).epochs{j},1)
                    f(an).trials{j}{k} = filterresult{f(an).epochs{j}(k,1)}{f(an).epochs{j}(k,2)};                                    
                end
            end
        else
            error('Trial filter must either be a cell array of length 1 or the same length as the number of groups');
        end
        
    end
    f(an).arguments.trials = filterInput;
    
end

