function f = singlecellsave(f, varargin)
% f = singlecellsave(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variable designated as strings in
% f().function.loadvariables{1}.  Also the function call appends any
% options in the f().function.options{} cell array.  
% 
% Each function call is for one cell (or one tetrode), and it is assumed that 
% the function's first input is the index to the cell or tetrode ([day epoch tetrode
% cell]).  The second input is a list of exclusion periods [starttime endtime].
% The next input is a single variable to be loaded
% out = fname(index, excludeperiods, var1);
%
% This function applies the exclude periods to the loaded variables and returns
% the resulting cell array.
%
%

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
    % check to see that only one variable was to be loaded
    if length(f(an).function.loadvariables) > 1
	error('load', 'singlecellsave can only load a single spike variable');
    end
    eval([f(an).function.loadvariables{1},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{1}, totaldays);']);
    loadstring = [loadstring, f(an).function.loadvariables{1},','];
    foptions = f(an).function.options;
    
    %iterate through the epochs within each data group
    for g = 1:length(f(an).epochs)
        for e = 1:size(f(an).epochs{g},1)
	    d = f(an).epochs{g}(e,1);
	    etmp = f(an).epochs{g}(e,2);
            for c = 1:size(f(an).data{g}{e},1)
		t = f(an).data{g}{e}(c,1);
		cnum = f(an).data{g}{e}(c,2);
                tmpindex = [f(an).epochs{g}(e,:) f(an).data{g}{e}(c,:)];
                excludeperiods = f(an).excludetime{g}{e};
                %run the designated function: fout = fname(tmpindex, var1, var2, ..., option1, option2, ...)
                eval(['fout = ',f(an).function.name,'(tmpindex,excludeperiods,', loadstring, 'foptions{:});']);
                
                %the result of this application will be a new spikestruct for
		%this animal
		if (length(f) == 1)
		    f(an).output{d}{etmp}{t}{cnum} = fout;
		else
		    f(an).output{d}{etmp}{t}{cnum} = fout;
		end
            end
        end
    end
end
