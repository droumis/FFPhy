
function out = makeExpvarCatDesignMat(lfpstack, varargin)
    % given a set of animal, day, epoch, timestamp.. evaluate the state filters
    % return design matrix of experiment variables [rip var]
    % Demetris Roumis 2019
    pconf = paramconfig;
    expvars = {'onlywdays', 'rewarded', 'unrewarded', 'inbound' , 'outbound', ...
         'proximalWell', 'distalWell'};
    saveout = 1;
    outdir = 'expvarCat';
    defaults = {'wtrackdays', 'excludeNoise','excludePriorFirstWell'};
    if ~isempty(varargin)
        assign(varargin{:})
    end
    
    outpath = [pconf.andef{2},outdir,'/'];
    for ian = 1:length(lfpstack)
        animal = lfpstack(ian).animal;
        out(ian).animal = animal;
        out(ian).dims = {'ripple', 'expvar'};
        out(ian).ripStartTime = lfpstack(ian).ripStartTime;
        out(ian).ripEndTime = lfpstack(ian).ripEndTime;
        out(ian).dayeps = [lfpstack(ian).day lfpstack(ian).epoch];
        out(ian).expvars = expvars;
        out(ian).dm = zeros(length(out(ian).ripStartTime), length(expvars));
        Fp = struct;
        for ss = 1:length(expvars)
            switch expvars{ss}
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
                case 'proximalWell'
                    Fp.add_params = {defaults{:}, 'proximalWell'};
                case 'distalWell'
                    Fp.add_params = {defaults{:}, 'distalWell'};
            end
            Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
            F = createfilter('animal', {out(ian).animal}, 'epochs', ...
                Fp.epochfilter, 'excludetime', Fp.timefilter);
            for de = 1:length(F.epochs{1}(:,1))
                day = F.epochs{1}(de,1);
                epoch = F.epochs{1}(de,2);
                ieprips = ismember(out(ian).dayeps, [day epoch], 'rows');
                out(ian).dm(ieprips,ss) = ...
                    ~isExcluded(out(ian).ripStartTime(ieprips), ...
                    F.excludetime{1}{de});
            end
        end
    end
    if saveout
        save_data(out, outpath, [outdir,'_',Fp.epochEnvironment]); 
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