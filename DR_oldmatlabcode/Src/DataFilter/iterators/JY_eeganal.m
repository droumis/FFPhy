function f = JY_eeganal(f)
% f = singleeeganal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.
%
% Each function call is for one tetrode, and it is assumed that
% the function's first input is the index to the tetrode ([day epoch tetrode]).
% The second input is a list of exclusion periods [starttime endtime].
% The next inputs are the load variables.  Note that eeg load variables,
% specified in iseegvar(), are loaded individually for each tetrode.
% The final inputs are the options.
% out = fname(index, excludeperiods, var1, var2, ..., option1, option2,...).
%
% The output of the call function can either be a 1 by N vector, or a structure.
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.

%iterate through all animals
for an = 1:length(f)
    %find all unique epochs to analyze for the current animal
    animaldir = f(an).animal{2};
    animalprefix = f(an).animal{3};
    totalepochs = [];
    for g = 1:length(f(an).epochs)
        totalepochs = [totalepochs; f(an).epochs{g}];
    end
    totaldays = unique(totalepochs(:,1)); %get all of the days across groups
    tepochs=unique(totalepochs(:,2));
    %load all the variables that the function requires
    foptions = f(an).function.options;
    
    %iterate through the epochs within each data group
    for g = 1:length(f(an).epochs)
        
        for e = 1:size(f(an).epochs{g},1)
            
            for c = 1:size(f(an).eegdata{g}{e},1)
                tmpindex = [f(an).epochs{g}(e,:) f(an).eegdata{g}{e}(c,:)];
                excludeperiods = f(an).excludetime{g}{e};
                loadstring = [];
                for i = 1:length(f(an).function.loadvariables)
                    % differ from eeganal
                    if strcmp(f(an).function.loadvariables{i},'eeg')
                        
                        eval([f(an).function.loadvariables{i},' = loadeegstruct(animaldir, animalprefix, f(an).function.loadvariables{i}, tmpindex(1),tmpindex(2),tmpindex(4));']);
                        loadstring = [loadstring, f(an).function.loadvariables{i},','];
                    elseif strcmp(f(an).function.loadvariables{i},'meaneegspectrograms')
                        
                        eval([f(an).function.loadvariables{i},' = loadeegstruct(animaldir, animalprefix, f(an).function.loadvariables{i}, tmpindex(1),tmpindex(2),tmpindex(4));']);
                        loadstring = [loadstring, f(an).function.loadvariables{i},','];   
                        
                        
                    elseif strcmp(f(an).function.loadvariables{i},'theta')
                        
                        eval([f(an).function.loadvariables{i},' = loadeegstruct(animaldir, animalprefix, f(an).function.loadvariables{i}, tmpindex(1),tmpindex(2),tmpindex(4));']);
                        loadstring = [loadstring, f(an).function.loadvariables{i},','];
 
                    elseif strcmp(f(an).function.loadvariables{i},'ripples')

                        eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, tmpindex(1));']);
                        loadstring = [loadstring, f(an).function.loadvariables{i},','];
                    elseif  strcmp(f(an).function.loadvariables{i},'tetinfo')

                        eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, tmpindex(1));']);
                        loadstring = [loadstring, f(an).function.loadvariables{i},',']; 
                        elseif  strcmp(f(an).function.loadvariables{i},'data')

                        eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, tmpindex(1));']);
                        loadstring = [loadstring, f(an).function.loadvariables{i},',']; 
                    
                                            elseif  strcmp(f(an).function.loadvariables{i},'linpos')

                        eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, tmpindex(1));']);
                        loadstring = [loadstring, f(an).function.loadvariables{i},',']; 
                    
                    
                    
                    elseif strcmp(f(an).function.loadvariables{i},'ripple')

                        eval([f(an).function.loadvariables{i},' = loadeegstruct(animaldir, animalprefix, f(an).function.loadvariables{i}, totaldays,tepochs,[]);']);
                        loadstring = [loadstring, f(an).function.loadvariables{i},','];
           
                     elseif strcmp(f(an).function.loadvariables{i},'gamma')

                        eval([f(an).function.loadvariables{i},' = loadeegstruct(animaldir, animalprefix, f(an).function.loadvariables{i},tmpindex(1),tmpindex(2),tmpindex(4));']);
                        loadstring = [loadstring, f(an).function.loadvariables{i},','];
           


else 

                        eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, tmpindex(1), tmpindex(2));']);
                        loadstring = [loadstring, f(an).function.loadvariables{i},',']; 
                    
                    end
                end
                %run the designated function: fout = fname(tmpindex, var1, var2, ..., option1, option2, ...)
                eval(['fout = ',f(an).function.name,'(tmpindex,excludeperiods,', loadstring, 'foptions{:});']);
                
                %save the function output in the filter variable.  Allows numeric or struct outputs
                if isstruct(fout)
                    if (isempty(f(an).output) | (length(f(an).output) < g))
                        f(an).output{g}(1) = fout;
                    else
                        f(an).output{g}(end+1) = fout;
                    end
                elseif isnumeric(fout)
                    if ((isempty(f(an).output)) | (length(f(an).output) < g))
                        f(an).output{g} = [];
                    end
                    if (size(fout,1) == 1)
                        f(an).output{g} = [f(an).output{g}; fout];
                    else
                        error(['In calling ', f(an).function.name, ': Numeric function outputs must be 1 by N.  Use a structure output for more complicated outputs']);
                    end
                else
                    error(['In calling ', f(an).function.name, ': Function output must be either numeric or a structure']);
                end
            end
        end
    end
end




