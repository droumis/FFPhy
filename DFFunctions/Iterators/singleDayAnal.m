function f = singleDayAnal(f,varargin)
% f = singleDayAnal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.
%
%                 Cerulean City
%                        .|
%                        | |
%                        |'|            ._____
%                ___    |  |            |.   |' .---"|
%        _    .-'   '-. |  |     .--'|  ||   | _|    |
%     .-'|  _.|  |    ||   '-__  |   |  |    ||      |
%     |' | |.    |    ||       | |   |  |    ||      |
%  ___|  '-'     '    ""       '-'   '-.'    '`      |____

%{
% Notes:
%   - city.alien
% Each function call is for all epochs in a day, and it is assumed that
% the function's first input is the index to the day, then all epochs ([day ep1 ep2..]).  
% The second input is a list of exclusion periods [starttime endtime].
% The next inputs are the load variables, and the final inputs are the options.
% out = fname(index, excludeperiods, var1, var2, ..., option1, option2,...).
%
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.

% @DKR
%}

%iterate through all animals
for a = 1:length(f)
    f(a).output = [];
    %find all unique days
    animaldir = f(a).animal{2};
    animal = f(a).animal{3};
    unqEpochs = [];
    for g = 1:length(f(a).epochs)
        unqEpochs = [unqEpochs; f(a).epochs{g}];
    end
    unqDays = unique(unqEpochs(:,1)); %get all of the days across groups
    
    %load all the variables that the function requires except the eeg
    loadstring = [];
    for i = 1:length(f(a).function.loadvariables)
        eval([f(a).function.loadvariables{i}, ...
    ' = loaddatastruct(animaldir, animal, f(a).function.loadvariables{i}, unqDays);']);
        loadstring = [loadstring, f(a).function.loadvariables{i},','];
    end
    % iterate through the days within each data group
    g = 1; % this was intended for multiple epoch filter groups, but isn't currently used?
    fprintf(':::::::: single Day iterator (Cerulean City):::::::: \n');
    tmp = [];
    fout = cell(numel(unqDays),1);
    for d = 1:numel(unqDays) % can use parfor
        day = unqDays(d);
        foptions = f(a).function.options; % reset per day of data
        foptions = [foptions {'animal', animal}];
        deIdx = find(f(a).epochs{g}(:,1)==unqDays(d));
        % collect day's timefilter
        excludeperiods = [];
        for e = 1:length(deIdx)
            ep = f(a).epochs{g}(deIdx(e),2);
            excludeperiods = [excludeperiods; f(a).excludetime{g}{e}];
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
        eps = f(a).epochs{g}(deIdx,2)';
        index = [day eps];
        fprintf('--%s %s day %d :: eps %s\n', f(a).function.name, ...
            f(a).animal{1}, day, strjoin(num2cell(num2str(eps(:)))',' '));
        % run the specified filter function on this day
        fout{d,1} = feval(f(a).function.name, index, excludeperiods, foptions{:});
    end
    f(a).output{g} = [fout{:}];
end
end