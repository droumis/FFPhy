function f = singlecellanal(f)
% f = singlecellanal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.  
% 
% Each function call is for one cell (or one tetrode), and it is assumed that 
% the function's first input is the index to the cell or tetrode 
% ([day epoch tetrode cell] or [day epoch tetrode]).  
% The second input is a list of exclusion periods [starttime endtime].
% The next inputs are the load variables.  
%
% 	Note that a load variable of 'filteroutput' specifies that the current
% 	element of the filter output field is fed into the function.  This
% 	allows for sequential analyses.
%
% The final inputs are the options. 
% out = fname(index, excludeperiods, var1, var2, ..., option1, option2,...).
%
% The output of the call function can either be a 1 by N vector, or a structure.
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.

%iterate through all animals
for an = 1:length(f)
    output = [];
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
    usefilteroutput = 0;
    for i = 1:length(f(an).function.loadvariables)
	% one load variable ('filteroutput') should not be loaded but should
	% instead be taken as the corresponding output element from the current
	% filter
	if (~strcmp(f(an).function.loadvariables{i}, 'filteroutput'))
	    eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, totaldays);']);
	else
	    usefilteroutput = 1;
	end
	loadstring = [loadstring, f(an).function.loadvariables{i},','];
    end
    foptions = f(an).function.options;
    
    %iterate through the epochs within each data group
    for g = 1:length(f(an).epochs)
        
        for e = 1:size(f(an).epochs{g},1)
            
            for c = 1:size(f(an).data{g}{e},1)
                tmpindex = [f(an).epochs{g}(e,:) f(an).data{g}{e}(c,:)];
                excludeperiods = f(an).excludetime{g}{e};
                %run the designated function: fout = fname(tmpindex, var1, var2, ..., option1, option2, ...)
		% check to see if we should use the filter output as an input
		% variable
		if (usefilteroutput)
		    if (isstruct(f(an).output{g}(c)))
			filteroutput = f(an).output{g}(c);
		    else
			filteroutput = f(an).output{g}(c,:);
		    end
		end
                eval(['fout = ',f(an).function.name,'(tmpindex,excludeperiods,', loadstring, 'foptions{:});']);
                
                %save the function output in the filter variable.  Allows numeric or struct outputs
                if isstruct(fout)
                    if (isempty(output))
                        output = fout;
                    else
                        output(end+1) = fout;
                    end
                elseif isnumeric(fout)
                    if (isempty(output))
                        output = [];
                    end
                    if (size(fout,1) == 1)
                        output = [output; fout];
                    else
                        error(['In calling ', f(an).function.name, ': Numeric function outputs must be 1 by N.  Use a structure output for more complicated outputs']);
                    end
                else
                    error(['In calling ', f(an).function.name, ': Function output must be either numeric or a structure']);
                end
            end
        end
	% assign the output
	f(an).output{g} = output;
	output = [];
    end
end
                        
                
            
            
