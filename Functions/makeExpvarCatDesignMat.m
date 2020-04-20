
function [out] = makeExpvarCatDesignMat(data, expvars, varargin)
% [out] = makeExpvarCatDesignMat(data, expvars, varargin)
% - evaluate expvars timefilter groups for data.animal, day, epoch, evStart
%   where expvars label timefilt funcs and filtstrings from load_filter_params
% - the output design matrix used to filter event-indexed data array
%
% args:
% - data: struct: animal, evStart, day, epoch
% - expvars: cell array of strings
% varargs:
%
% output:
% - out.dm: design matrix of eventSet (event, set)

%{
Beer
FFPhy
@DR
%}

pconf = paramconfig;
% expvars = {'all'}; %, 'rewarded', 'unrewarded', 'inbound' , 'outbound', ...
%          'proximalWell', 'distalWell', 'rewarded_outbound', 'rewarded_inbound'};
saveout = 1;
outdir = 'expvarCat';
defaults = {'wtrackdays', 'excludePriorFirstWell', 'excludeAfterLastWell'};
eventType = 'swr';

if ~isempty(varargin)
    assign(varargin{:})
end
fprintf('defaults: %s\n', defaults{:})
outpath = [pconf.andef{2},outdir,'/'];

% if data comes from the spiking pipeline, need to convert to F with dayep, evtimes
if strcmp(data(1).iterator, 'singleDayCellAnal')
    for a = 1:length(data)
        dF(a).animal = data(a).animal;
        % get dayeps idx of one cell per day
        idx = cell2mat({data(a).output{1}.index}');
        [~, dayUnqIdx] = unique(idx(:,1)); 
        days = idx(dayUnqIdx,[1]);
        eps = idx(dayUnqIdx,[4:5]);
        % get num events per dayep. 
        numEvPerDay = cell2mat({data(a).output{1}(dayUnqIdx).numEventsPerEp})';
        % collect eventTimes across all days
        dF(a).evStart = cell2mat({data(a).output{1}(dayUnqIdx).eventTimes}');
        % label each the event time with day ep
        dF(a).day = [];
        dF(a).epoch = [];
        for d = 1:size(eps,1)
            for e = 1:size(eps,2)
                dF(a).day = [dF(a).day; repmat(days(d), numEvPerDay(d,e), 1)];
                dF(a).epoch = [dF(a).epoch; repmat(eps(d,e), numEvPerDay(d,e),1)];
            end
        end
    end
    data = dF;
end

for ian = 1:length(data)
    %         t = find(strcmp(lfpstack(ian).lfptypes, lfptype));
    try
        animal = data(ian).animal{3};
    catch
        animal = animaldef(data(ian).animal);
    end
    out(ian).animal = animal;
    out(ian).dims = {'event', 'expvar'};
    out(ian).evStart = data(ian).evStart;
    %         out(ian).evEnd = lfpstack(ian).evEnd{t};
    out(ian).dayeps = [data(ian).day data(ian).epoch];
    out(ian).expvars = expvars;
    out(ian).dm = zeros(length(out(ian).evStart), length(expvars));
    Fp = struct;
    for ss = 1:length(expvars)
        Fp.params = defaults;
        switch expvars{ss}
            case 'all'
                % all defaults
            case 'wetLickBursts'
                Fp.params{end+1} = 'wetLickBursts';
            case 'dryLickBursts'
                Fp.params{end+1} = 'dryLickBursts';                
            case 'lickbouts'
                Fp.params{end+1} = 'lickbouts';
            case 'nolickbouts'
                Fp.params{end+1} = 'nolickbouts';
            case 'rewarded'
                Fp.params{end+1} = 'correcttrials';
            case 'unrewarded'
                Fp.params{end+1} = 'errortrials';
                %                 case 'outbound'
                %                     Fp.params{end+1} = 'outbound';
                %                 case 'inbound'
                %                     Fp.params{end+1} = 'inbound';
            case 'rewarded_outbound'
                Fp.params{end+1} = 'correcttrials';
                Fp.params{end+1} = 'outbound';
            case 'unrewarded_outbound'
                Fp.params{end+1} = 'errortrials';
                Fp.params{end+1} = 'outbound';
            case 'rewarded_inbound'
                Fp.params{end+1} = 'correcttrials';
                Fp.params{end+1} = 'inbound';
            case 'unrewarded_inbound'
                Fp.params{end+1} = 'errortrials';
                Fp.params{end+1} = 'inbound';
                %                 case 'proximalWell'
                %                     Fp.params{end+1} = 'proximalWell';
                %                 case 'distalWell'
                %                     Fp.params{end+1} = 'distalWell';
            otherwise
                Fp.params{end+1} = expvars{ss};
        end
        fprintf('============ %s =============\n', expvars{ss});
        Fp = load_filter_params(Fp);
        F = createfilter('animal', {animal}, 'epochs', ...
            Fp.epochfilter, 'excludetime', Fp.timefilter);
        for de = 1:length(F.epochs{1}(:,1))
            day = F.epochs{1}(de,1);
            epoch = F.epochs{1}(de,2);
            % events in current epoch
            ievs = ismember(out(ian).dayeps, [day epoch], 'rows');
            % events in current timefilter
            out(ian).dm(ievs,ss) = ...
                ~isExcluded(out(ian).evStart(ievs), ...
                F.excludetime{1}{de});
        end
    end
end
if saveout
    save_data(out, outpath, [outdir,'_',Fp.env,'_',eventType]);
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