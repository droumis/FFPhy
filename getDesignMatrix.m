
function designmat = getDesignMatrix(ripstruct, varargin)
    % Given a set of animal, day, epoch, timestamps.. return design matrix
    % of experiment variables: 
    % Demetris Roumis 2019
    
    if ~isempty(varargin)
        assign(varargin{:})
    end
    
    expvars = {'timeSinceDay', 'timeSinceEpoch','xpos', 'ypos', 'headdirection',  'speed',...
        'performance', 'learningrate', 'timeSinceLastReward', 'timeUntilNextReward'};
    
    numvars = length(expvars);
    designmat = struct;
    for ian = 1:length(ripstruct)
        animal = ripstruct(ian).animal;
        andef = animaldef(animal);
        designmat(ian).animal = animal;
        numrips = length(ripstruct(ian).ripStartTime);
        designmat(ian).ripStartTime = ripstruct(ian).ripStartTime;
        designmat(ian).ripEndTime = ripstruct(ian).ripEndTime;
        designmat(ian).dayeps = [ripstruct(ian).day ripstruct(ian).epoch];
        designmat(ian).expvars = expvars;
        designmat(ian).dm = nan(numrips,numvars);
        % timeSinceDay is just the ripstarttimes
        designmat(ian).dm(:,1) = ripstruct(ian).ripStartTime;
        % get position, speed vars
        dayeps = unique(designmat(ian).dayeps,'rows', 'stable');
        load(sprintf('%s/%sBehaveState.mat',andef{1,2}, animal));
        allbound = BehaveState.statespace.allbound;
        for ide = 1:length(dayeps(:,1))
            day = dayeps(ide,1);
            epoch = dayeps(ide,2);
            iderips = ismember(designmat(ian).dayeps, [day epoch], 'rows');
            ideripstarts = designmat(ian).ripStartTime(iderips);
            load(sprintf('%s%s%s%02d.mat',andef{1,2}, animal, 'pos', day));
            ripidx = knnsearch(pos{day}{epoch}.data(:,1), ideripstarts);
            % time since epoch
            epochstartime = pos{day}{epoch}.data(1,1);
            designmat(ian).dm(iderips,2) = ideripstarts - epochstartime;
            
            % pos vars
            designmat(ian).dm(iderips,[3 4 5 6]) = pos{day}{epoch}.data(ripidx,[6 7 8 9]);
       
            % get performance and reward vars       
            allbound_inep_inds = ismember(allbound(:,[5 6]), [day epoch], 'rows');
            bscorrect = logical(BehaveState.statespace.allepsMat(allbound_inep_inds,7));
            bstimes = BehaveState.statespace.allepsMat(allbound_inep_inds,3);
            if isempty(bstimes)
                fprintf('no behavestate %s day %d epoch %d. skipping\n',animal,day,epoch);
                continue
            end
            bsmode = allbound(allbound_inep_inds,1);
            bsmodediff = abs(diff(allbound(allbound_inep_inds,1)));
            designmat(ian).dm(iderips,[7]) = interp1(bstimes, bsmode, ideripstarts);
            designmat(ian).dm(iderips,[8]) = interp1(bstimes(2:end), bsmodediff,ideripstarts);
            
            % time since last reward
            iseidxrips = find(iderips);
            rewardtimes = bstimes(bscorrect);
            timesincelast = arrayfun(@(x) x-rewardtimes(max(find(rewardtimes<x))), ideripstarts, 'un', 0);
            usel = cellfun(@(x) ~isempty(x),timesincelast,'un',1);
            designmat(ian).dm(iseidxrips(usel),9) = cell2mat(timesincelast(usel));
            % time until next reward
            timeuntilnext = arrayfun(@(x) x-rewardtimes(find(rewardtimes>x,1)), ideripstarts, 'un', 0);
            usen = cellfun(@(x) ~isempty(x),timeuntilnext,'un',1);
            designmat(ian).dm(iseidxrips(usen),10) = cell2mat(timeuntilnext(usen));
        end
    end
end