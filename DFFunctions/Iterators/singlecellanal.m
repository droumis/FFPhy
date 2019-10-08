function f = singlecellanal(f, varargin)
% f = singlecellanal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.
%
% Each function call is for one cell (or one tetrode), and it is assumed that
% the function's first input is the index to the cell or tetrode ([day epoch tetrode
% cell]).  The second input is a list of exclusion periods [starttime endtime].
% The next inputs are the load variables, and the final inputs are the options.
% out = fname(index, excludeperiods, var1, var2, ..., option1, option2,...).
%
% The output of the call function can either be a 1 by N vector, or a structure.
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.
returnStruct = 1;
if ~isempty(varargin)
    assign(varargin{:})
end

for an = 1:length(f)
    %find all unique epochs to analyze for the current animal
    animaldir = f(an).animal{2};
    animal = f(an).animal{3};
    totalepochs = [];
    for g = 1:length(f(an).epochs)
        totalepochs = [totalepochs; f(an).epochs{g}];
    end
    totaldays = unique(totalepochs(:,1)); %get all of the days across groups
    
    %load all the variables that the function requires
    loadstring = [];
    foptions = f(an).function.options;
    for i = 1:length(f(an).function.loadvariables)
%         try
        eval([f(an).function.loadvariables{i}, ...
' = loaddatastruct(animaldir, animal, f(an).function.loadvariables{i}, totaldays);']);
%         catch
%         eval([f(an).function.loadvariables{i}, ...
% ' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i});']);
%         end
        loadstring = [loadstring, f(an).function.loadvariables{i},','];
    end
    %iterate through the epochs within each data group
    fprintf(':::::::: single cell iterator :::::::: \n');
    for g = 1:length(f(an).epochs) % what is this intended for?
        for e = 1:size(f(an).epochs{g},1)
            day = f(an).epochs{g}(e,1);
            ep = f(an).epochs{g}(e,2);
            excludeperiods = f(an).excludetime{g}{e};
            numcells = size(f(an).data{g}{e},1);
            fout = cell(1,numcells);
            % just pass in this epoch's worth of data to each worker
            d = [];
            foptions = f(an).function.options;
            for i = 1:length(f(an).function.loadvariables)
                try
                    d{day}{ep} = eval([f(an).function.loadvariables{i} '{day}{ep}']);
                    eval(['foptions = [foptions {f(an).function.loadvariables{i},d}];']);
                catch
%                     continue
                    error('cant load f(an).function.loadvariables{i}\n')
                end
            end
            for c = 1:numcells % can use parfor
                tet = f(an).data{g}{e}(c,1);
                clust = f(an).data{g}{e}(c,2);
                cindex = [day ep f(an).data{g}{e}(c,:)];
                %         foptions = [foptions {'an', an}]; %add the animal index to the varargins
                %         if outputDayEpTetCells
                %             eval(['f(an).output{' num2str(ep) '}{' num2str(tet) ...
                %                 '} = ' f.function.name,'(tmpindex,excludeperiods,', loadstring, ...
                %                 'foptions{:});']);
                %         else
                %run the designated function: fout = fname(tmpindex, var1, var2, ..., option1, option2, ...)
                fprintf(sprintf('%s %d %d %d %d:: %s \n', f(an).animal{1}, day, ...
                    ep, tet, clust, f(an).function.name));
                        % run the specified filter function on this set of animal/epoch/ntrodes
                foptions = [foptions {'animal',animal}];
                fout{c} = feval(f(an).function.name, cindex, excludeperiods, foptions{:});
                %save the function output in the filter variable.  Allows numeric or struct outputs
%                 if isstruct(fout)
%                     if (isempty(f(an).output) | (length(f(an).output) < g))
%                         f(an).output{g}(1) = fout;
%                     else
%                         f(an).output{g}(end+1) = fout;
%                     end
%                 elseif isnumeric(fout)
%                     if ((isempty(f(an).output)) | (length(f(an).output) < g))
%                         f(an).output{g} = [];
%                     end
%                     if (size(fout,1) >= 1) %edited by AS changed from (size(fout,1) == 1)
%                         f(an).output{g} = stack(f(an).output{g}, fout);
%                     else
%                         error(['In calling ', f(an).function.name, ': Numeric function outputs must be 1 by N.  Use a structure output for more complicated outputs']);
%                     end
%                 else
%                     error(['In calling ', f(an).function.name, ': Function output must be either numeric or a structure']);
%                 end
            end
            f(an).output{g}{e} = fout;
        end
        f(an).output{g} = f(an).output{g}{:};
        % i need to have an option to save as struct array, numeric array, cell array
%         if returnStruct
%             f(an).output{g} = [f(an).output{g}{:}];
%         end
    end
end

end

