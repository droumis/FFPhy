function f = singletetrodeanal(f, varargin)
% f = multicellanal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.  
% 
% Each function call is for one epoch, and it is assumed that the
% function's first input is a list of indices to the tetrode
% ([day epoch tetrode]). The second input is a list of exclusion
% periods [starttime endtime].  The next inputs are the load
% variables, and the final inputs are the options.  out =
% fname(index, excludeperiods, var1, var2, ..., option1,
% option2,...).
%
% The output of the call function can either be numeric, or a structure.
% If the output if numeric, data is appended along the first dimenion for
% the group.
outputDayEpTetCells = 0; %output into a nested cell structure {day}{ep}{tet}
allowpairs = 1;
if ~isempty(varargin)
    assign(varargin{:})
end

%iterate through all animals
for an = 1:length(f)
   %find all unique epochs to analyze for the current animal
   animaldir = f(an).animal{2};
   animalprefix = f(an).animal{3};
   totalepochs = [];
   totaltetrodes = [];
   for g = 1:length(f(an).epochs)
      totalepochs = [totalepochs; f(an).epochs{g}];
      for e = 1:size(f(an).epochs{g},1)
         totaltetrodes = [totaltetrodes; f(an).data{g}{e}];
      end
   end
   totaldays = unique(totalepochs(:,1)); %get all of the days across groups
   totaltetrodes = unique(totaltetrodes);
   
   %load all the variables that the function requires
   loadstring = [];
   for i = 1:length(f(an).function.loadvariables)
        if (iseegvar(f(an).function.loadvariables{i}))
         eval([f(an).function.loadvariables{i},' = loadeegstruct(animaldir, animalprefix, f(an).function.loadvariables{i}, totaldays);']);
        else
         eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, totaldays);']);
        end
      loadstring = [loadstring, f(an).function.loadvariables{i},','];
   end
   foptions = f(an).function.options;
    
    %iterate through the epochs within each data group
    fprintf(':::::::: tetrode iterator :::::::: \n');
    for g = 1:length(f(an).epochs)
    for e = 1:size(f(an).epochs{g},1)
    for c = 1:size(f(an).data{g}{e},1)
        tmpindex = [f(an).epochs{g}(e,:) f(an).data{g}{e}(c,:)];
        excludeperiods = f(an).excludetime{g}{e};

        day = tmpindex(1);
        ep = tmpindex(2);
        tet = tmpindex(3:end);
        foptions = [foptions {'an', an}]; %add the animal index to the varargins
        if outputDayEpTetCells
            eval(['f.output{' num2str(day) '}{' num2str(ep) '}{' num2str(tet) ...
                '} = ',f.function.name,'(tmpindex,excludeperiods,', loadstring, ...
                'foptions{:});']);
        else
            %run the designated function: fout = fname(tmpindex, var1, var2, ..., option1, option2, ...)
            fprintf(sprintf('%s %d %d %d :: %s \n', f(an).animal{1}, day, ...
                ep, tet, f(an).function.name));
            eval(['fout = ',f(an).function.name,'(tmpindex,excludeperiods,', ...
                loadstring, 'foptions{:});']);
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
                if (size(fout,1) >= 1) %edited by AS changed from (size(fout,1) == 1)
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
end
