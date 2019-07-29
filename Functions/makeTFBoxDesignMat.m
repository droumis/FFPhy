function out = makeTFBoxDesignMat(rawpwr,Fp, varargin)
% makes a design matrix of feature values per ripple
pconf = paramconfig;
tfboxlabels = {'swrRipple', 'swrFastGamma', 'swrSlowGamma','swrBeta','swrTheta',...
    'preSwrRipple', 'preSwrFastGamma', 'preSwrSlowGamma',  'preSwrBeta', 'preSwrTheta',...
    'postSwrRipple', 'postSwrFastGamma', 'postSwrSlowGamma', 'postSwrBeta', ...
    'postSwrTheta', 'prePreSwrRipple', 'prePreSwrFastGamma', 'prePreSwrSlowGamma', ...
    'prePreSwrBeta', 'prePreSwrTheta', 'postPostSwrRipple', 'postPostSwrFastGamma', ...
    'postPostSwrSlowGamma', 'postPostSwrBeta', 'postPostSwrTheta'};

saveout = 1;
outdir = 'tfbvarCont'; 
if ~isempty(varargin)
    assign(varargin{:});
end
outpath = [pconf.andef{2},outdir,'/'];
animals = {rawpwr.animal};
for ani = 1:length(animals)
    [out(ani).freq, out(ani).time] = tfboxlookup(tfboxlabels);
    out(ani).animal = animals{ani};
    out(ani).expvars = {};
    out(ani).dims = {'ripple', 'tfbvar', 'ntrode'};
    out(ani).expvars = tfboxlabels;
    for tfb = 1:length(tfboxlabels)
        out(ani).timeidx{tfb} = dsearchn(rawpwr(ani).time', out(ani).time{tfb}');
        out(ani).freqidx{tfb} = dsearchn(rawpwr(ani).frequency', out(ani).freq{tfb}');
        out(ani).realfreq{tfb} = rawpwr(ani).frequency(out(ani).freqidx{tfb}(1):...
            out(ani).freqidx{tfb}(2));
        out(ani).dm(:,tfb,:) = permute(squeeze(mean(mean(rawpwr(ani).pwr(:, ...
        out(ani).timeidx{tfb}(1):out(ani).timeidx{tfb}(2),:, ...
        out(ani).freqidx{tfb}(1):out(ani).freqidx{tfb}(2)), 2), 4)), [2 1 3]); % mean TF
    end
end
if saveout
    save_data(out, outpath, [outdir,'_',Fp.epochEnvironment]);
end

function [freq, time] = tfboxlookup(tfboxlabels)

frexlookup.Ripple = [125 250];
frexlookup.FastGamma = [125 250];
frexlookup.SlowGamma = [25 50];
frexlookup.Beta = [12 20];
frexlookup.Theta = [6 12];

timelookup.prePreSwr = [-.6 -.2];
timelookup.preSwr = [-.2 0];
timelookup.Swr = [0 .2];
timelookup.postSwr = [.2 .4];
timelookup.postPostSwr = [.4 .8];

freq = {};
time = {};
for i = 1:length(tfboxlabels)
    if strfind(tfboxlabels{i},'Ripple')
        freq{i} = frexlookup.Ripple;   
    elseif strfind(tfboxlabels{i},'FastGamma')
        freq{i} = frexlookup.FastGamma;
    elseif strfind(tfboxlabels{i},'SlowGamma')
        freq{i} = frexlookup.SlowGamma;
    elseif strfind(tfboxlabels{i},'Beta')
        freq{i} = frexlookup.Beta;
    elseif strfind(tfboxlabels{i},'Theta')
        freq{i} = frexlookup.Theta;
    end
    if strfind(tfboxlabels{i},'prePreSwr')
        time{i} = timelookup.prePreSwr;
    elseif strfind(tfboxlabels{i},'preSwr')
        time{i} = timelookup.preSwr;
    elseif strfind(tfboxlabels{i},'postPostSwr')
        time{i} = timelookup.postPostSwr;
    elseif strfind(tfboxlabels{i},'postSwr')
        time{i} = timelookup.postSwr ;
    elseif strfind(tfboxlabels{i},'swr') % lowercase swr
        time{i} = timelookup.Swr;
    end
end