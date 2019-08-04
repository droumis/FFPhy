
function out = makeExpvarContDesignMat(lfpstack, Fp, varargin)
    % Given a set of animal, day, epoch, timestamps.. return design matrix
    % of experiment variables: 
    % Demetris Roumis 2019
    
    pconf = paramconfig;
    saveout = 1;
    outdir = 'expvarCont'; 
    if ~isempty(varargin)
        assign(varargin{:})
    end
    outpath = [pconf.andef{2},outdir,'/'];    
    expvars = {'speed','performance', 'learningrate', 'day'}; 
    %,, 'epoch''ripnum',  'timeSinceDay', 'timeSinceEpoch','xpos', 'ypos', 'headdirection', 'timeSinceLastReward', 'timeUntilNextReward', 'ripnum'};
    
    
    out = struct;
    for ian = 1:length(lfpstack)
        animal = lfpstack(ian).animal;
        andef = animaldef(animal);
        out(ian).animal = animal;
        numrips = length(lfpstack(ian).ripStartTime);
        out(ian).ripStartTime = lfpstack(ian).ripStartTime;
        out(ian).ripEndTime = lfpstack(ian).ripEndTime;
        out(ian).dayeps = [lfpstack(ian).day lfpstack(ian).epoch];
        out(ian).expvars = expvars;
        out(ian).dims = {'ripple', 'expvar'};
        out(ian).dm = nan(numrips,length(expvars));
        % timeSinceDay is just the ripstarttimes
        out(ian).dm(:,1) = lfpstack(ian).ripStartTime;
        % get position, speed vars
        dayeps = unique(out(ian).dayeps,'rows', 'stable');
        load(sprintf('%s/%sBehaveState.mat',andef{1,2}, animal));
        allbound = BehaveState.statespace.allbound;
        out(ian).dm(:,9) = 1:numrips;
        out(ian).dm(:,10) = lfpstack(ian).day;
        out(ian).dm(:,11) = lfpstack(ian).epoch;
        for ide = 1:length(dayeps(:,1))
            day = dayeps(ide,1);
            epoch = dayeps(ide,2);
            iderips = ismember(out(ian).dayeps, [day epoch], 'rows');
            ideripstarts = out(ian).ripStartTime(iderips);
            load(sprintf('%s%s%s%02d.mat',andef{1,2}, animal, 'pos', day));
            ripidx = knnsearch(pos{day}{epoch}.data(:,1), ideripstarts);
            % time since epoch
            epochstartime = pos{day}{epoch}.data(1,1);
            out(ian).dm(iderips,2) = ideripstarts - epochstartime;
            
            % pos vars (x, y, hd, sp) ~ (6,7,8,9)
            out(ian).dm(iderips,[3 4 5 6]) = pos{day}{epoch}.data(ripidx,[6 7 8 9]);
       
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
            out(ian).dm(iderips,[7]) = interp1(bstimes, bsmode, ideripstarts);
            out(ian).dm(iderips,[8]) = interp1(bstimes(2:end),bsmodediff,ideripstarts);
            % these have too many nans and trim the data too much
%             % time since last reward
%             iseidxrips = find(iderips);
%             rewardtimes = bstimes(bscorrect);
%             timesincelast = arrayfun(@(x) x-rewardtimes(max(find(rewardtimes<x))), ideripstarts, 'un', 0);
%             usel = cellfun(@(x) ~isempty(x),timesincelast,'un',1);
%             designmat(ian).dm(iseidxrips(usel),9) = cell2mat(timesincelast(usel));
%             % time until next reward
%             timeuntilnext = arrayfun(@(x) x-rewardtimes(find(rewardtimes>x,1)), ideripstarts, 'un', 0);
%             usen = cellfun(@(x) ~isempty(x),timeuntilnext,'un',1);
%             designmat(ian).dm(iseidxrips(usen),10) = cell2mat(timeuntilnext(usen));
        end
    end
    if saveout
        save_data(out, outpath, [outdir,'_',Fp.epochEnvironment]); 
    end
end