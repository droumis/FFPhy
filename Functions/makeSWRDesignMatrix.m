
function out = makeSWRDesignMatrix(lfpstack, Fp, varargin)
% Given a set of animal, day, epoch, timestamps.. return design matrix
% of experiment variables related to swr
% Demetris Roumis 2019

pconf = paramconfig;
saveout = 1;
outdir = 'swrvarCont';
if ~isempty(varargin)
    assign(varargin{:})
end
outpath = [pconf.andef{2},outdir,'/'];
expvars = {'std', 'duration'}; %'total_energy', 
for ian = 1:length(lfpstack)
    animal = lfpstack(ian).animal;
    andef = animaldef(animal);
    out(ian).animal = animal;
    numrips = length(lfpstack(ian).ripStartTime);
    out(ian).ripStartTime = lfpstack(ian).ripStartTime;
    out(ian).ripEndTime = lfpstack(ian).ripEndTime;
    out(ian).dayeps = [lfpstack(ian).day lfpstack(ian).epoch];
    out(ian).expvars = expvars;
    out(ian).dims = {'ripple', 'swrvar'};
    out(ian).dm = nan(numrips,length(expvars));
    ripks = loaddatastruct(andef{1,2}, animal, 'ca1rippleskons');
    dayeps = unique(out(ian).dayeps,'rows', 'stable');
    for ide = 1:length(dayeps(:,1))
        day = dayeps(ide,1);
        epoch = dayeps(ide,2);
        iderips = ismember(out(ian).dayeps, [day epoch], 'rows');
        ideripstarts = out(ian).ripStartTime(iderips);
        ripidx = knnsearch(ripks{day}{epoch}{1}.starttime, ideripstarts);
        % var std
        out(ian).dm(iderips,1) = ripks{day}{epoch}{1}.maxthresh(ripidx);
%         out(ian).dm(iderips,2) = ripks{day}{epoch}{1}.total_energy(ripidx);
        % var duration
        out(ian).dm(iderips,2) = ripks{day}{epoch}{1}.endtime(ripidx) - ...
            ripks{day}{epoch}{1}.starttime(ripidx);
    end
end
if saveout
    save_data(out, outpath, [outdir,'_',Fp.epochEnvironment]);
end
end