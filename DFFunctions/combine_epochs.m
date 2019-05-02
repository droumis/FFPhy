
function ppF = combine_epochs(F,Fp, saveCombinedEpochs, paths, varargin)
% combine the data for a given day tet cell across epochs of the same type
% switch flag for different filter output types
% paths is an object from the constructor make_paths.m
ppF = struct;
if ~isempty(varargin)
    assign(varargin{:})
end

switch Fp.filtfunction
    case Fp.filtfunction

for iAn = 1:length(Fp.animals)

    % load anim/tet info 
    animalinfo = animaldef(lower(Fp.animals{iAn}));
    animalID = animalinfo{1,3}; %use anim prefix for name
    FFanimdir =  sprintf('%s',animalinfo{1,2});
    load([FFanimdir, animalID, 'tetinfo']);
    tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);

    % find the unique day/tet outputs
    matInds = cell2mat({F.output{1}.index}');
    [daytetcells, ~, daytetInds2 ] = unique(matInds(:,[1 3 4]), 'rows', 'stable');

    % ---------- for each cell, get all epochs--------------------------------
    % for each allepoch-unique cell
    for ic = 1:size(daytetcells,1);
        cellID = daytetcells(ic, 3);
        % indices into (F)ilter output for this day tet cell
        icellFoutInds = find(ic == daytetInds2);
        epochInds = matInds(icellFoutInds,:);
        
        numeps = size(epochInds,1);
        iInd = epochInds(1,:);
        day = iInd(1);
        epochF = iInd(2);
        ntrode = iInd(3);

        load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'task', day));
        % ---------- get unique epoch environment types  --------------------------------
        eptypes = [];
        for iepoch = 1:numeps;
            epoch = epochInds(iepoch,2);
            eptypes{iepoch,1} = task{day}{epoch}.environment;
        end
        [envtypes,~,IndC] = unique(eptypes, 'stable');
        ppF = initialize_riptrigspiking_output(ppF, ic, numeps, daytetcells(ic,:));
        
        for ienv = 1:size(envtypes,1) %for each environment type
            ienvTypeInds = find(ienv == IndC);
            ienvFInds = icellFoutInds(ienvTypeInds);
            
            % the business end
            for iepienv = 1:length(ienvFInds)
                epID = epochInds(ienvTypeInds(iepienv),2);
                iout = F(iAn).output{1}(ienvFInds(iepienv));
                ppF(ic).time = iout.time;   % (all time vectors are the same..)
                ppF(ic).frtime = iout.frtime;
                ppF(ic).psth = [ppF(ic).psth ; iout.psth];
                ppF(ic).frhist = [ppF(ic).frhist; iout.frhist];
                ppF(ic).instantFR = [ppF(ic).instantFR ; iout.instantFR];
                ppF(ic).posteventmatrix = [ppF(ic).posteventmatrix ; iout.posteventmatrix];
                ppF(ic).eventduration = [ppF(ic).eventduration ; iout.eventduration];
                ppF(ic).eventtags = [ppF(ic).eventtags ; iout.eventtags];
                ppF(ic).nospikes = ppF(ic).nospikes + iout.nospikes;
                ppF(ic).noevents = ppF(ic).noevents + iout.noevents;
                ppF(ic).epochs = [ppF(ic).epochs epID];
                % epoch-by-epoch outputs
                ppF(ic).epoch_types{epID} = iout.epoch_type;
                ppF(ic).epoch_envs{epID} = iout.epoch_environment;
                ppF(ic).epoch_nospikes(epID) = sum(iout.nospikes);
                ppF(ic).epoch_noevents(epID) = sum(iout.noevents);
                ppF(ic).epoch_noeventspikes(epID) = sum(sum(iout.psth));
            end
            % if a no data for this cell, ignore
            if isempty(ppF(ic).psth)
                fprintf('psth is empty \n')
                continue
            end
            ppF(ic).psthsum = sum(ppF(ic).psth,1);
            ppF(ic).instantFRmean = mean(ppF(ic).instantFR,1);
        end
        fprintf(sprintf('combined epochs by env type for cell %d of %d \n',ic, ...
            size(daytetcells,1)))
    end
end
end

% ---------------- Save combined epochs ---------------------------------------------------
if saveCombinedEpochs == 1;
    save_filter_output(ppF, paths.filtOutputDirectory, paths.filenamesave, ...
        'filetail', '_combEps')
end
end

function ppF = initialize_riptrigspiking_output(ppF, ic, numeps, daytetcell)
%% ---------- for each epoch env type, concat the data --------------------------------
% initialize output .fields
ppF(ic).dtc = daytetcell;
ppF(ic).epoch_types = [];            % indices here correspond to epoch #s
ppF(ic).epoch_envs = [];             % indices here correspond to epoch #s
% all detected events, regardless of epoch or state
ppF(ic).psth = [];
ppF(ic).frhist = [];
ppF(ic).instantFR = [];
ppF(ic).psthsum = [];
ppF(ic).instantFRmean = [];
ppF(ic).posteventmatrix = [];
ppF(ic).eventduration = [];
ppF(ic).eventtags = [];
ppF(ic).eventtags_descript = '[ epochnum eventtime epochtype]';
ppF(ic).nospikes = 0;
ppF(ic).noevents = 0;                  % number of events reported for the epochs in which the unit was clustered ("events experienced")
ppF(ic).epochs = [];
ppF(ic).run_epochs = [];
ppF(ic).run_nospikes = nan(1,numeps);
ppF(ic).run_noevents = nan(1,numeps);
ppF(ic).sleep_epochs = [];
ppF(ic).sleep_nospikes = nan(1,numeps);
ppF(ic).sleep_noevents = nan(1,numeps);
ppF(ic).time = [];
end