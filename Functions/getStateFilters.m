
function ripstate = getStateFilters(lfpstack, varargin)
    % given a set of animal, day, epoch, timestamp.. evaluate the state filters
    % Demetris Roumis 2019
    
    statefilters = {'onlywdays', 'rewarded', 'unrewarded', 'inbound' , 'outbound'};
    
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
        Fp = struct;
        defaults = {'wtrackdays', 'excludeNoise','excludePriorFirstWell'};
        for ss = 1:length(statefilters)
            switch statefilters{ss}
%                 case 'all'
%                     Fp.add_params = {'excludeNoise','excludePriorFirstWell'};
                case 'onlywdays'
                    Fp.add_params = defaults;
                case 'rewarded'
                    Fp.add_params = {defaults{:}, 'correcttrials'};
                case 'unrewared'
                    Fp.add_params = {defaults{:}, 'errortrials'};
                case 'outbound'
                    Fp.add_params = {defaults{:}, 'outbound'};
                case 'inbound'
                    Fp.add_params = {defaults{:}, 'inbound'};
                case 'rewarded_outbound'
                    Fp.add_params = {defaults{:}, 'correcttrials', 'outbound'};
                case 'unrewarded_outbound'
                    Fp.add_params = {defaults{:}, 'errortrials', 'outbound'};
                case 'rewarded_inbound'
                    Fp.add_params = {defaults{:}, 'correcttrials', 'inbound'};
                case 'unrewarded_inbound'
                    Fp.add_params = {defaults{:}, 'errortrials', 'inbound'};
            end
            Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
            F = createfilter('animal', {ripstate(ian).animal}, 'epochs', ...
                Fp.epochfilter, 'excludetime', Fp.timefilter);
            for de = 1:length(F.epochs{1}(:,1))
                day = F.epochs{1}(de,1);
                epoch = F.epochs{1}(de,2);
                ieprips = ismember(ripstate(ian).dayeps, [day epoch], 'rows');
                ripstate(ian).statesets(ieprips,ss) = ...
                    ~isExcluded(ripstate(ian).ripStartTime(ieprips), ...
                    F.excludetime{1}{de});
            end
        end
    end
end

% function [state, fields] = getState(animal, dayeps, times)
% %% this hsould be a time filter
% andef = animaldef(animal);
% behavestate = load(sprintf('%s/%s%s.mat',andef{2}, animal, 'BehaveState'));
% % linpos = loaddatastruct(andef{2}, animal, 'linpos');
% % pos = loaddatastruct(andef{2}, animal, 'pos');
% 
% unqdayeps = unique(dayeps, 'rows');
% state = zeros(length(times), 2);
% fields = {'correct', 'outbound'}; %, 'linvelocity', '2dvelocity'};
% for de = 1:length(unqdayeps(:,1))
%     day =  unqdayeps(de,1);
%     ep = unqdayeps(de,2);
%     deidx = find(ismember(dayeps, [day ep], 'rows'));
%     detimes = times(deidx);
%     trialIO = behavestate.BehaveState.statechanges{day}{ep}.statechangeseq;
%     trialIOfields = behavestate.BehaveState.statechanges{day}{ep}.fields;
%     colfields = strsplit(trialIOfields, ' ');
%     corrcol = find(cellfun(@(x) strcmp(x,'correct'), colfields, 'un', 1));
%     portcol = find(cellfun(@(x) strcmp(x,'timeportout'), colfields, 'un', 1));
%     outBcol = find(cellfun(@(x) strcmp(x,'outbound'), colfields, 'un', 1));
%     inBcol = find(cellfun(@(x) strcmp(x,'inbound'), colfields, 'un', 1));
%     lastTcol = find(cellfun(@(x) strcmp(x,'lasttime'), colfields, 'un', 1));
%     currTcol = find(cellfun(@(x) strcmp(x,'currenttime'), colfields, 'un', 1));
%     
%     corrStartEnd = trialIO(trialIO(:,corrcol)==1,[lastTcol, currTcol]);
%     corr = isExcluded(detimes, corrStartEnd);
%     
%     outBStartEnd = trialIO(trialIO(:,outBcol)==1,[lastTcol, currTcol]);
%     outB = isExcluded(detimes, outBStartEnd);
% 
%     state(deidx,:) = [corr outB];
% end
% end