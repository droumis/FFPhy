function out = makeTFBoxDesignMat(rawpwr, varargin)
% makes a design matrix of time/frequency boxes' mean power per ripple

pconf = paramconfig;
tfboxlabels = {'Ripple', 'FastGamma', 'Theta'};
epEnv = 'wtrack';
saveout = 1;
outdir = 'tfbvarCont';
frexlookup.Ripple = [125 250];
frexlookup.FastGamma = [65 115];
frexlookup.SlowGamma = [25 55];
% frexlookup.Beta = [12 20];
frexlookup.Theta = [6 12];
timelookup.baseline = [-1 -.5];
timelookup.event = [0 .15];
timelookup.postEvent= [.5 1];

if ~isempty(varargin)
    assign(varargin{:});
end
outpath = [pconf.andef{2},outdir,'/'];
animals = {rawpwr.animal};
for ian = 1:length(animals)
    out(ian).freq = {}; out(ian).time = {};  out(ian).baseline = {};
    for i = 1:length(tfboxlabels)
        out(ian).freq{i} = eval(['frexlookup.', tfboxlabels{i}]);
        if any(strfind(tfboxlabels{i},'FastGamma')) || any(strfind(tfboxlabels{i},'Theta'))
            out(ian).baseline{i} = timelookup.baseline;
            out(ian).time{i} = timelookup.postEvent;
        elseif any(strfind(tfboxlabels{i},'SlowGamma')) || any(strfind(tfboxlabels{i},'Ripple'))
            out(ian).baseline{i} = timelookup.baseline;
            out(ian).time{i} = timelookup.event;
        end
    end
    %     [out(ian).freq, out(ian).time out(ian).baseline] = tfboxlookup(tfboxlabels);
    out(ian).animal = animals{ian};
    out(ian).expvars = {};
    out(ian).dims = {'ripple', 'tfbvar', 'ntrode'};
    out(ian).expvars = tfboxlabels;
    for tfb = 1:length(tfboxlabels)
        timeidx = dsearchn(rawpwr(ian).time', out(ian).time{tfb}');
        baselidx = dsearchn(rawpwr(ian).time', out(ian).baseline{tfb}');
        freqidx = dsearchn(rawpwr(ian).frequency', out(ian).freq{tfb}');
        out(ian).timeidx{tfb} = timeidx;
        out(ian).baselidx{tfb} = baselidx;
        out(ian).freqidx{tfb} = freqidx;
        out(ian).realfreq{tfb} = rawpwr(ian).frequency(out(ian).freqidx{tfb}(1):...
            out(ian).freqidx{tfb}(2));
        % rawpwr dims [ntrode time ripple freq].. take mean in time, frex ranges
        pwrtfb = rawpwr(ian).pwr(:,timeidx(1):timeidx(2),:,freqidx(1):freqidx(2));
        baselinepwrtfb = rawpwr(ian).pwr(:,baselidx(1):baselidx(2),:,freqidx(1):freqidx(2));
        mpwrtfb = zscore(nanmean(nanmean(pwrtfb, 2), 4)-nanmean(nanmean(baselinepwrtfb,2),4),[],3);
        out(ian).dm(:,tfb,:) = squeeze(mpwrtfb)'; %reshape to [ripple ntrode]
        % dm dims [rip tfb ntrode]
    end
end
if saveout
    save_data(out, outpath, [outdir,'_',epEnv]); end
