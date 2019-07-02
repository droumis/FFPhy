

% given an animal, day, epoch, timestamp.. evaluate the state of each type

function ripstate = getStateFilters(lfpstack)
    statefilters = {'all', 'rewarded', 'unrewarded', 'inbound' , 'outbound', 'rewarded_inbound', ...
        'unrewarded_inbound', 'rewarded_outbound', 'unrewarded_outbound'};
    
    ripstate = struct;
    for ian = 1:length(lfpstack)
        animal = lfpstack(ian).animal;
        ripstate(ian).animal = animal;
        ripstate(ian).ripStartTime = lfpstack(ian).ripStartTime;
        ripstate(ian).ripEndTime = lfpstack(ian).ripEndTime;
        ripstate(ian).dayeps = [lfpstack(ian).day lfpstack(ian).epoch];
        [ripstate(ian).state, ripstate(ian).statefields] = getState(animal, ...
            ripstate(ian).dayeps, ripstate(ian).ripStartTime);
        
        ripstate(ian).statesetsfields = statefilters;
        ripstate(ian).statesets = zeros(length(ripstate(ian).ripStartTime), ...
            length(statefilters));
        for ss = 1:length(statefilters)
            colfields = ripstate(ian).statefields;
            outbcol = find(cellfun(@(x) strcmp(x,'outbound'), colfields, 'un', 1));
            corrcol = find(cellfun(@(x) strcmp(x,'correct'), colfields, 'un', 1));
            switch statefilters{ss}
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