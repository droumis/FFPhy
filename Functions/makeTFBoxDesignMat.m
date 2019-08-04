function out = makeTFBoxDesignMat(rawpwr,Fp, varargin)
% makes a design matrix of feature values per ripple
pconf = paramconfig;
tfboxlabels = {'swrRipple', 'postPostSwrTheta', 'prePreSwrTheta', 'postPostSwrFastGamma', ...
    'prePreSwrFastGamma'};

%,'swrFastGamma'}; %, ...
%     'swrSlowGamma','swrBeta','swrTheta'};%,...
%     'preSwrRipple', 'preSwrFastGamma', 'preSwrSlowGamma',  'preSwrBeta', 'preSwrTheta',...
%     'postSwrRipple', 'postSwrFastGamma', 'postSwrSlowGamma', 'postSwrBeta', ...
%     'postSwrTheta', 'prePreSwrRipple',  'prePreSwrSlowGamma', ...
%     'prePreSwrBeta', 'postPostSwrRipple', , ...
%     'postPostSwrSlowGamma', 'postPostSwrBeta', };

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
        timeidx = dsearchn(rawpwr(ani).time', out(ani).time{tfb}');
        freqidx = dsearchn(rawpwr(ani).frequency', out(ani).freq{tfb}');
        out(ani).timeidx{tfb} = timeidx;
        out(ani).freqidx{tfb} = freqidx;
        out(ani).realfreq{tfb} = rawpwr(ani).frequency(out(ani).freqidx{tfb}(1):...
            out(ani).freqidx{tfb}(2));
        % rawpwr dims [ntrode time ripple freq].. take mean in time, frex ranges
        pwrtfb = rawpwr(ani).pwr(:,timeidx(1):timeidx(2),:,freqidx(1):freqidx(2));
        mpwrtfb = nanmean(nanmean(pwrtfb, 2), 4);
        out(ani).dm(:,tfb,:) = squeeze(mpwrtfb)'; %reshape to [ripple ntrode]
        % dm dims [rip tfb ntrode]
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

timelookup.prePreSwr = [-1 -.5];
timelookup.preSwr = [-.3 0];
timelookup.Swr = [0 .1];
timelookup.postSwr = [.1 .4];
timelookup.postPostSwr = [.5 1];

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