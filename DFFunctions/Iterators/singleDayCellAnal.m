function f = singleDayCellAnal(f,varargin)
% f = singleDayCellAnal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.
%
% Each function call is for one cluster, across epochs in a day, and it is assumed that
% the function's first input is the index to the cell ([day ntrode cell]).  
% The second input is a list of exclusion periods [starttime endtime].
% The next inputs are the load variables, and the final inputs are the options.
% out = fname(index, excludeperiods, var1, var2, ..., option1, option2,...).
%
% The output of the call function can either be a 1 by N vector, or a structure.
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.
% $DR19

%iterate through all animals
for a = 1:length(f)
    %find all unique days
    animaldir = f(a).animal{2};
    animal = f(a).animal{3};
    foptions = f(a).function.options;
    foptions = [foptions {'animal', animal}];
    totalepochs = [];
    for g = 1:length(f(a).epochs)
        totalepochs = [totalepochs; f(a).epochs{g}];
    end
    unqdays = unique(totalepochs(:,1)); %get all of the days across groups
    
    %load all the variables that the function requires except the eeg
    loadstring = [];
    for i = 1:length(f(a).function.loadvariables)
        eval([f(a).function.loadvariables{i}, ...
    ' = loaddatastruct(animaldir, animal, f(a).function.loadvariables{i}, unqdays);']);
        loadstring = [loadstring, f(a).function.loadvariables{i},','];
    end
    % iterate through the days within each data group
    g = 1;% what is this intended for?
    fprintf(':::::::: single Day Cell iterator :::::::: \n');
    tmp = [];
    for d = 1:length(unqdays)
        day = unqdays(d);
        deIdx = find(f(a).epochs{g}(:,1)==unqdays(d));
        excludeperiods = [];
        % add this cell's full day of data to foptions
        % fuck i need to fix this.. looks like cells are not defined for
        % whole day automatically.. 
        tmp = [];
        for i = 1:length(f(a).function.loadvariables)
            for e = 1:length(deIdx)
                ep = f(a).epochs{g}(deIdx(e),2);
                excludeperiods = [excludeperiods; f(a).excludetime{g}{e}];
                numcells = size(f(a).data{g}{deIdx(e)},1);
                % add loaded vars to varargin
                tmp{day}{ep} = eval([f(a).function.loadvariables{i} '{day}{ep}']);
            end
            eval(['foptions = [foptions {f(a).function.loadvariables{i}, tmp}];']);
        end
        
        fout = cell(numcells,1);
        for c = 1:numcells % can use parfor
            tet = f(a).data{g}{deIdx(e)}(c,1);
            clust = f(a).data{g}{deIdx(e)}(c,2);
            eps = f(a).epochs{g}(deIdx,2)';
            cindex = [day f(a).data{g}{deIdx(e)}(c,:) eps];
            fprintf(sprintf('%s %s day%d nt%d cl%d :: eps %d %d\n', f(a).function.name, ...
                f(a).animal{1}, day, tet, clust, eps));
            % run the specified filter function on this set of animal/epoch/ntrodes
            fout{c} = feval(f(a).function.name, cindex, excludeperiods, foptions{:});
        end
        f(a).output{g}{d} = fout;
    end
    f(a).output{g} = f(a).output{g}{:};
end
end