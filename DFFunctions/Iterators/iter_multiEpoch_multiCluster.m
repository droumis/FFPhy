function f = iter_multiEpoch_multiCluster(f, varargin)
% f = iter_multiEpoch_multiCluster(f)
% Iterator for a filter object. Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariable{:}. Also the function call appends any options
% in the f().function.options{} cell array.

% iterate througn all animals
for a = 1:length(f)
    f(a).output = [];
    % find all unique days
    animaldir = f(a).animal{2};
    animal = f(a).animal{3};
    unqEpochs = [];
    for g = 1:length(f(a).epochs)
        unqEpochs = [unqEpochs; f(a).epochs{g}];
    end
    unqDays = unique(unqEpochs(:,1)); % get unique days across filter groups
    
    % load all the variables that the function requires except the eeg
    loadstring = [];
    for i = 1:length(f(a).function.loadvariables)
        % i hate that this still uses eval
        eval([f(a).function.loadvariables{i}, ...
            '= loaddatastruct(animaldir, animal, f(a).function.loadvariables{i}, unqDays);']);
    end
    
    % Iterate through the days within each data group
    g = 1; % this was intended for multiple epoch filter groups, but isn't currently used?
    fprintf(':::::::: multi epoch multi cell iterator :::::::: \n');
    
    for d = 1:length(unqDays)
        day = unqDays(d);
        foptions = f(a).function.options; %reset per day of data
        foptions = [foptions {'animal', animal}]; % add animal to options
        deIdx = find(f(a).epochs{g}(:,1)==unqDays(d)); % day epoch index
        % collect day's timefilter
        excludeperiods = [];
        for e = 1:length(deIdx)
            ep = f(a).epochs{g}(deIdx(e),2);
            excludeperiods = [excludeperiods; f(a).excludetime{g}(e)];
        end
        % collect day's data
        tmp = [];
        for i = 1:length(f(a).function.loadvariables)
            for e = 1:length(deIdx)
                ep = f(a).epochs{g}(deIdx(e),2);
                % add loaded vars to varargin
                tmp{day}{ep} = eval([f(a).function.loadvariables{i} '{day}{ep}']);
            end
            eval(['foptions = [foptions {f(a).function.loadvariables{i}, tmp}];']);
        end
        % get unq cells
        ntClust = [];
        for e = 1:length(deIdx)
            ntClust = [ntClust; f(a).data{g}{deIdx(e)}];
        end
        unqNtCell = unique(ntClust, 'rows');
        numCells = size(unqNtCell,1);
        % evaluate function per cell pair
        fout = cell(numCells,1);
        for c = 1:numCells % can use parfor (parallelizes cells within this day)
            ntA = unqNtCell(c,1);
            clustA = unqNtCell(c,2);
            ntB = unqNtCell(c,3);
            clustB = unqNtCell(c,4);
            eps = f(a).epochs{g}(deIdx,2)';
            cindex = [day ntA clustA ntB clustB eps];
            fprintf('--%s %s day %d nt %d cl %d :: eps %s\n', f(a).function.name, ...
                f(a).animal{1}, day, ntA, clustA, ntB, clustB,...
                strjoin(num2cell(num2str(eps(:)))',' '));
            % run the specified filter function on this day cell
            fout{c,1} = feval(f(a).function.name, cindex, excludeperiods, foptions{:});
        end
        f(a).output{g}{d,1} = [fout{:}];
    end
    f(a).output{g} = [f(a).output{g}{:}];
end
end
