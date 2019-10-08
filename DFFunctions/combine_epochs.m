
function ppF = combine_epochs(F,Fp, saveCombinedEpochs, paths, varargin)
% combine the data for a given day tet cell across epochs of the same type
% switch flag for different filter output types
% paths is an object from the constructor make_paths.m
ppF = struct;
if ~isempty(varargin)
    assign(varargin{:})
end

% switch Fp.filtfunction
%     case Fp.filtfunction
ippF = struct;
for ian = 1:length(Fp.animals)
    
    % load anim/tet info
    animalinfo = animaldef(lower(Fp.animals{ian}));
    animal = animalinfo{1,3}; %use anim prefix for name
    FFanimdir =  sprintf('%s',animalinfo{1,2});
    load([FFanimdir, animal, 'tetinfo']);
    tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
    
    % find the unique day/tet outputs
    idata = [F.output{1}{:}];
    matInds = cell2mat({idata.index}');
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
        
        load(sprintf('%s%s%s%02d.mat',FFanimdir, animal, 'task', day));
        % ---------- get unique epoch environment types  --------------------------------
        eptypes = [];
        for iepoch = 1:numeps;
            epoch = epochInds(iepoch,2);
            eptypes{iepoch,1} = task{day}{epoch}.environment;
        end
        [envtypes,~,IndC] = unique(eptypes, 'stable');
        if ic ==1
            ippF = initialize_riptrigspiking_output(daytetcells(ic,:));
        else
            ippF = [ippF initialize_riptrigspiking_output(daytetcells(ic,:))];
        end
        for ienv = 1:size(envtypes,1) %for each environment type
            ienvTypeInds = find(ienv == IndC);
            ienvFInds = icellFoutInds(ienvTypeInds);
            
            % the business end
            for iepienv = 1:length(ienvFInds)
                epID = epochInds(ienvTypeInds(iepienv),2);
                iout = idata(ienvFInds(iepienv));
                ippF(ic).time = iout.time;   % (all time vectors are the same..)
                ippF(ic).frtime = iout.frtime;
                ippF(ic).psth = [ippF(ic).psth ; iout.psth];
                ippF(ic).frhist = [ippF(ic).frhist; iout.frhist];
                ippF(ic).instantFR = [ippF(ic).instantFR ; iout.instantFR];
                ippF(ic).posteventmatrix = [ippF(ic).posteventmatrix ; iout.posteventmatrix];
                ippF(ic).eventduration = [ippF(ic).eventduration ; iout.eventduration];
                ippF(ic).eventtags = [ippF(ic).eventtags ; iout.eventtags];
                ippF(ic).nospikes = ippF(ic).nospikes + iout.nospikes;
                ippF(ic).noevents = ippF(ic).noevents + iout.noevents;
                ippF(ic).epochs = [ippF(ic).epochs epID];
                % epoch-by-epoch outputs
                ippF(ic).epoch_types{epID} = iout.epoch_type;
                ippF(ic).epoch_envs{epID} = iout.epoch_environment;
                ippF(ic).epoch_nospikes(epID) = sum(iout.nospikes);
                ippF(ic).epoch_noevents(epID) = sum(iout.noevents);
                ippF(ic).epoch_noeventspikes(epID) = sum(sum(iout.psth));
            end
            % if a no data for this cell, ignore
            if isempty(ippF(ic).psth)
                fprintf('psth is empty \n')
                continue
            end
            ippF(ic).psthsum = sum(ippF(ic).psth,1);
            ippF(ic).instantFRmean = mean(ippF(ic).instantFR,1);
        end
        fprintf(sprintf('combined epochs by env type for cell %d of %d \n',ic, ...
            size(daytetcells,1)))
    end
    ppF(ian).animal = animal;
    ppF(ian).data = ippF;
end
% end
% ---------------- Save combined epochs ---------------------------------------------------
if saveCombinedEpochs
    save_data(ppF, paths.filtOutputDirectory, paths.filenamesave, ...
        'filetail', '_combEps')
end
end
function p = initialize_riptrigspiking_output(daytetcell)
%% ---------- for each epoch env type, concat the data --------------------------------
% initialize output .fields
p.time = [];
p.frtime = [];
p.psth = [];
p.frhist = [];
p.instantFR = [];
p.posteventmatrix = [];
p.eventduration = [];
p.eventtags = [];
p.nospikes = 0;
p.noevents = 0;                  % number of events reported for the epochs in which the unit was clustered ("events experienced")
p.epochs = [];
p.epoch_types = [];            % indices here correspond to epoch #s
p.epoch_envs = [];             % indices here correspond to epoch #s
p.epoch_nospikes = [];
p.epoch_noevents = [];
p.epoch_noeventspikes = [];
p.psthsum = [];
p.instantFRmean = [];

p.dtc = daytetcell;
p.eventtags_descript = '[ epochnum eventtime epochtype]';

end
