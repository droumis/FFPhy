

% given an animal, day, epoch, timestamp.. evaluate the state of each type

function ripstate = getStateFilters(lfpstack, varargin)
    statefilters = {'all', 'onlywdays', 'rewarded', 'unrewarded', 'inbound' , 'outbound', 'rewarded_inbound', ...
        'unrewarded_inbound', 'rewarded_outbound', 'unrewarded_outbound'};
    
    if ~isempty(varargin)
        assign(varargin{:})
    end
    
    ripstate = struct;
    for ian = 1:length(lfpstack)
        animal = lfpstack(ian).animal;
        ripstate(ian).animal = animal;
        ripstate(ian).ripStartTime = lfpstack(ian).ripStartTime;
        ripstate(ian).ripEndTime = lfpstack(ian).ripEndTime;
        ripstate(ian).dayeps = [lfpstack(ian).day lfpstack(ian).epoch];
        ripstate(ian).statesetsfields = statefilters;
        ripstate(ian).statesets = zeros(length(ripstate(ian).ripStartTime), ...
            length(statefilters));
        
        % get a bunch of the behavioral states first.. probably should make
        % this into their own timefilter calls
        [ripstate(ian).state, ripstate(ian).statefields] = getState(animal, ...
            ripstate(ian).dayeps, ripstate(ian).ripStartTime);
        colfields = ripstate(ian).statefields;
        outbcol = find(cellfun(@(x) strcmp(x,'outbound'), colfields, 'un', 1));
        corrcol = find(cellfun(@(x) strcmp(x,'correct'), colfields, 'un', 1));
        
        for ss = 1:length(statefilters)
            switch statefilters{ss}
                case 'onlywdays'
                    % from the ripple times in the lfpstack.. determine which satisfy the
                    % conditions based on the F.epochs and F.excludetimes{t}{d}{e}
                    Fp.filtfunction = 'dfa_riptriglfp';
                    Fp.add_params = {'wtrackdays', 'excludeNoise','excludePriorFirstWell'};
                    Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
                    F = createfilter('animal', {ripstate(ian).animal}, 'epochs', ...
                        Fp.epochfilter, 'excludetime', Fp.timefilter);
                        ripstate(ian).ripStartTime = lfpstack(ian).ripStartTime;
                        includerips = zeros(length(ripstate(ian).ripStartTime),1);
                        for de = 1:length(F.epochs{1}(:,1))
                            day = F.epochs{1}(de,1);
                            epoch = F.epochs{1}(de,2);
                            ieprips = ismember(ripstate(ian).dayeps, [day epoch], 'rows');
                            ripstate(ian).statesets(ieprips,ss) = ...
                                ~isExcluded(ripstate(ian).ripStartTime(ieprips), ...
                                F.excludetime{1}{de});
                        end
                case 'all'
                    ripstate(ian).statesets(:,ss) = 1+ripstate(ian).statesets(:,ss);
                case 'rewarded'
                    ripstate(ian).statesets(:,ss) = ripstate(ian).state(:,corrcol);
                case 'unrewared'
                    ripstate(ian).statesets(:,ss) = ~ripstate(ian).state(:,corrcol);
                case 'outbound'
                    ripstate(ian).statesets(:,ss) = ripstate(ian).state(:,outbcol);
                case 'inbound'
                    ripstate(ian).statesets(:,ss) = ~ripstate(ian).state(:,outbcol);
                case 'rewarded_outbound'
                    ripstate(ian).statesets(:,ss) = all([ripstate(ian).state(:,outbcol) ...
                        ripstate(ian).state(:,corrcol)],2);
                case 'unrewarded_outbound'
                    ripstate(ian).statesets(:,ss) = all([ripstate(ian).state(:,outbcol) ...
                        ~ripstate(ian).state(:,corrcol)],2);
                case 'rewarded_inbound'
                    ripstate(ian).statesets(:,ss) = all([~ripstate(ian).state(:,outbcol) ...
                        ripstate(ian).state(:,corrcol)],2);
                case 'unrewarded_inbound'
                    ripstate(ian).statesets(:,ss) = all([~ripstate(ian).state(:,outbcol) ...
                        ~ripstate(ian).state(:,corrcol)],2);
            end
        end
    end
end


function [state, fields] = getState(animal, dayeps, times)
%% this hsould be a time filter
andef = animaldef(animal);
behavestate = load(sprintf('%s/%s%s.mat',andef{2}, animal, 'BehaveState'));
% linpos = loaddatastruct(andef{2}, animal, 'linpos');
% pos = loaddatastruct(andef{2}, animal, 'pos');

unqdayeps = unique(dayeps, 'rows');
state = zeros(length(times), 2);
fields = {'correct', 'outbound'}; %, 'linvelocity', '2dvelocity'};
for de = 1:length(unqdayeps(:,1))
    day =  unqdayeps(de,1);
    ep = unqdayeps(de,2);
    deidx = find(ismember(dayeps, [day ep], 'rows'));
    detimes = times(deidx);
    trialIO = behavestate.BehaveState.statechanges{day}{ep}.statechangeseq;
    trialIOfields = behavestate.BehaveState.statechanges{day}{ep}.fields;
    colfields = strsplit(trialIOfields, ' ');
    corrcol = find(cellfun(@(x) strcmp(x,'correct'), colfields, 'un', 1));
    portcol = find(cellfun(@(x) strcmp(x,'timeportout'), colfields, 'un', 1));
    outBcol = find(cellfun(@(x) strcmp(x,'outbound'), colfields, 'un', 1));
    inBcol = find(cellfun(@(x) strcmp(x,'inbound'), colfields, 'un', 1));
    lastTcol = find(cellfun(@(x) strcmp(x,'lasttime'), colfields, 'un', 1));
    currTcol = find(cellfun(@(x) strcmp(x,'currenttime'), colfields, 'un', 1));
    
    corrStartEnd = trialIO(trialIO(:,corrcol)==1,[lastTcol, currTcol]);
    corr = isExcluded(detimes, corrStartEnd);
    
    outBStartEnd = trialIO(trialIO(:,outBcol)==1,[lastTcol, currTcol]);
    outB = isExcluded(detimes, outBStartEnd);

    state(deidx,:) = [corr outB];
end
end