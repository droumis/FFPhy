function f = multitetrodeanal(f, varargin)
% [f] = multicellanal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.
%
%                     Fecund Forest
%        ^  ^  ^   ^      ___I_      ^  ^   ^  ^  ^   ^  ^
%       /|\/|\/|\ /|\    /\-_--\    /|\/|\ /|\/|\/|\ /|\/|\
%       /|\/|\/|\ /|\   /  \_-__\   /|\/|\ /|\/|\/|\ /|\/|\
%       /|\/|\/|\ /|\   |[]| [] |   /|\/|\ /|\/|\/|\ /|\/|\
% funcs:
% - bear: dfa_eventTrigLFP

% Each function call is for one epoch, and it is assumed that
% the function's first input is a list of indices to the cell or tetrode ([day epoch tetrode
% cell]).  The second input is a list of exclusion periods [starttime endtime].
% The next inputs are the load variables, and the final inputs are the options.
% out = fname(index, excludeperiods, var1, var2, ..., option1, option2,...).
%
% The output of the call function can either be numeric, or a structure.
% If the output if numeric, data is appended along the first dimenion for
% the group.
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.

%iterate through all animals
for an = 1:length(f)

    %find all unique days/epochs
    animaldir = f(an).animal{2};
    animal = f(an).animal{3};
    unqdayepochs = [];
    for iepochs = 1:length(f(an).epochs)
        unqdayepochs = [unqdayepochs; f(an).epochs{1}];
    end
    unqdays = unique(f(an).epochs{1}(:,1)); %get all of the days across groups
    
    %load all the variables that the function requires except the eeg
    loadstring = [];
    for i = 1:length(f(an).function.loadvariables)
        eval([f(an).function.loadvariables{i}, ...
' = loaddatastruct(animaldir, animal, f(an).function.loadvariables{i}, unqdays);']);
        loadstring = [loadstring, f(an).function.loadvariables{i},','];
    end
    
    % run func per epoch, all ntrodes
    numepochs = length(f(an).epochs{1}(:,1));
    fout = cell(1,numepochs);
    for idayep = 1:numepochs % cannot use parfor bc eval is used within loop
        day = f(an).epochs{1}(idayep,1);
        epoch = f(an).epochs{1}(idayep,2);
        ntrodes = f(an).eegdata{1}{idayep};
        numntrodes = size(ntrodes,1);
        indices = [repmat(f(an).epochs{1}(idayep,:),[numntrodes 1]) ntrodes];
        excludeperiods = f(an).excludetime{1}{idayep};
        if isempty(indices)
            fprintf(sprintf('no data for %s d%d e%d.. skipping\n',animal,day,epoch));
            continue
        end
        
        % for this day/epoch load the eeg data for all ntrodes
        for i = 1:length(f(an).function.loadvariables)
            if (iseegvar(f(an).function.loadvariables{i}))
                iv = f(an).function.loadvariables{i};
                eval([f(an).function.loadvariables{i}, ...
                    '=loadeegstruct(animaldir,animal,iv, day, epoch, ntrodes);']);
            end
        end
        
        % add the eeg data vars to varargin cell array
        foptions = f(an).function.options;
        d = [];
        for i = 1:length(f(an).function.loadvariables)
            try
                d{day}{epoch} = eval([f(an).function.loadvariables{i} '{day}{epoch}']);
                eval(['foptions = [foptions {f(an).function.loadvariables{i},d}];']);
                continue
            end
        end
        
        % run the specified filter function on this set of animal/epoch/ntrodes
        fprintf(sprintf('%s %d %d:: %s \n', animal, day, epoch, f(an).function.name));
        foptions = [foptions {'animal',animal}];
        fout{idayep} = feval(f(an).function.name, indices, excludeperiods, foptions{:});
%         
%         eval(['fout = ',f(an).function.name,'(indices, excludeperiods,' loadstring, ...
%             'foptions{:});']);
        
%         if isempty(fout.data)
%             continue
%         end
        %save the function output in the filter variable.  Allows numeric or struct outputs
%         if isstruct(fout)
%             f(an).output{day}(epoch) = fout;
%             
% %             elseif isnumeric(fout)
% %                 if (length(f(an).output) < day)
% %                     f(an).output{day} = [];
% %                 end
% %                 f(an).output{day} = stack(f(an).output{day}, fout);
%         else
%             error(['In calling ', f(an).function.name, ': Function output must be either numeric or a structure']);
%         end
    end
    f(an).output = fout;
end
end