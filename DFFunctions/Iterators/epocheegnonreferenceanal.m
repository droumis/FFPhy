function f = epocheegnonreferenceanal(f, varargin)
% f = singleepochanal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.
%
% This iterator should be used when you want to use perform an analysis where
% you select a single eeg channel and behavioral data (pos, linpos, etc.).
%
% Each function call is for one epoch, and it is assumed that
% the function's first input is the index to the tetrode for that day ([day epoch tetrode]).  
% The second input is a list of exclusion periods [starttime endtime].
% The next inputs are the load variables, and the final inputs are the options.
% out = fname(index, excludeperiods, var1, var2, ..., option1, option2,...).
%
% The output of the call function can either be a 1 by N vector, or a structure.
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.

eegvar = [];

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

    %load all the variables that the function requires
    loadstring = [];
    for i = 1:length(f(an).function.loadvariables)
        if (~iseegvar(f(an).function.loadvariables{i}))
            eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, totaldays);']);
        end
        loadstring = [loadstring, f(an).function.loadvariables{i},','];
    end
    foptions = f(an).function.options;

    %iterate through the epochs within each data group
    for g = 1:length(f(an).epochs)

        for e = 1:size(f(an).epochs{g},1)
            
            for c = 1:size(f(an).eegdata{g}{e},1)
                index = [f(an).epochs{g}(e,:) f(an).eegdata{g}{e}(c,:)];
                excludeperiods = f(an).excludetime{g}{e};
                % load the eeg data for this day epoch
                for i = 1:length(f(an).function.loadvariables)
                    if (iseegvar(f(an).function.loadvariables{i}))
                        eval([f(an).function.loadvariables{i},' = loadeegnonreferencestruct(animaldir, animalprefix, f(an).function.loadvariables{i}, index(1), index(2), index(3:end));']);
                    end
                end
            

                %run the designated function: fout = fname(tmpindex, var1, var2, ..., option1, option2, ...)
                eval(['fout = ',f(an).function.name,'(index,excludeperiods,', loadstring, 'foptions{:});']);

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
                    if (size(fout,1) >= 0) %edited by AS changed from (size(fout,1) == 1)
                            f(an).output{g} = stack(f(an).output{g}, fout);
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

