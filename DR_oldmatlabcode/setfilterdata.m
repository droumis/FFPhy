function f = setfilterdata(f, filterStrings)
% f = setfiltercells(f, filterString)
% For each epoch in the filter F, this function finds the indices to the
% cells that satisfy the given filter condtions in filterString.  The syntax
% for filterString is defined in EVALUATEFILTER.m. The animal and desired epochs
% for the filter need to be predefined. Assumes that each animal's data
% folder contains a file 'cellinfo.mat' that contains a cell structure with
% information about each cell.

if ~iscell(filterStrings)
    error('The cell filter input must be a cell array');
end

for an = 1:length(f)
    if isempty(f(an).animal)
        error(['You must define an animal for the filter before filtering the cells'])
    end
    if isempty(f(an).epochs)
        error(['You must define the desired epochs for the filter before filtering the cells'])
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
    
    for g = 1:length(f(an).epochs)
        
            %call each filter function using the totalepochs list
            filterresults = [];
            if (length(filterStrings) == 1)
                i = 1;
            elseif (length(filterStrings) == length(f(an).epochs))
                i = g;
            else
                error('The number of data filer conditions must match the number of conditions set up with the epoch filter.');
            end
            
            %for i = 1:length(filterStrings) %one filter for each experimental condition
                if ((length(filterStrings{i}) >= 2) && (mod(length(filterStrings{i}),2) == 0))
                    for j = 1:2:length(filterStrings{i}) %for combined filters
                        if (strfind(lower(filterStrings{i}{j}),'<function>')) %this is a function call
                            fparam = parseFunctionCallString(filterStrings{i}{j});
                            dataStruct = feval(fparam.funcName, datadir, animalprefix, totalepochs, fparam.funcVarargin{:});
                            filterresults{i}{(j+1)/2} = evaluatefilter(dataStruct, filterStrings{i}{j+1}); %evaluate the filter string for every structure and return results in a cell array
                        else %assume the input is the name of a variable to search
                            dataStruct = loaddatastruct(datadir, animalprefix, filterStrings{i}{j}, loaddays);
                            tmpresult = evaluatefilter(dataStruct, filterStrings{i}{j+1}); %evaluate the filter string for every structure and return results in a cell array
                            tmpresultsorted = [];
                            for e = 1:size(f(an).epochs{g},1)
                                eind = find(ismember(tmpresult(:,1:2),f(an).epochs{g}(e,1:2),'rows'));
                                tmpresultsorted{e} = tmpresult(eind,3:end);
                            end
                            %f(an).data = setfield(f(an).data,{an},filterStrings{i}{j},{g},{tmpresultsorted});
                            f(an).data = setfield(f(an).data,filterStrings{i}{j},{g},{tmpresultsorted});
                            
                        end
                    end
                else
                    error('Each inputdata filter call must have inputs: {datatype1, searchstring1, datatype2, searchstring2, ... }');
                end
                
            %end
        
    end
                  
    f(an).arguments.data = filterStrings;
end