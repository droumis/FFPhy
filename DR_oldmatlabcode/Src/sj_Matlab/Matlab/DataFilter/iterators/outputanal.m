function f = outputanal(f)
% f = outputanal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name and loads the values in f().input which are designated 
% as strings in f().function.loadvariables{:}.  Also the function call 
% appends any options in the f().function.options{} cell array.
%
% This iterator should be used when you want to perform an analysis on
% output from a previous filter and have since moved to f().input. The 
% iterator is agnostic as to the data type, it just uses the output structure
% defined by the previous filter.
%
% Each function call is for one unit of output(whatever that may be), and 
% it is assumed that the function's first input is the index ([day epoch]).  
% The second input is a list of exclusion periods [starttime endtime].
% The next input is whatever is present in input for the index, and the final 
% inputs are the options.
% out = fname(index, excludeperiods, input, ..., option1, option2,...).
%
% The output of the call function can either be a 1 by N vector, or a structure.
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.
% 
% This has only been verified for eegdata, not sure it works on cell data
% mcarr july2009


%iterate through all animals
for an = 1:length(f)
    foptions = f(an).function.options;
    
    %iterate through the epochs within each data group
    for g = 1:length(f(an).epochs)
        
        for e = 1:size(f(an).epochs{g},1)

            index = [f(an).epochs{g}(e,:) f(an).eegdata{g}{e}];
            % load the data for this day epoch
            data = f(an).input{g}(e);
            excludeperiods = f(an).excludetime{g}{e};
           
            %run the designated function: fout = fname(tmpindex, excludeperiods, data, ..., option1, option2, ...)
            eval(['fout = ',f(an).function.name,'(index,excludeperiods,data,foptions{:});']);

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
                if (size(fout,1) >= 0)
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

