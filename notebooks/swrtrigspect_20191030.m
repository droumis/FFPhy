

%{
just get a per-animal, all-swr-trig spect.. focusing on the higher
frequencies
%}

% get data
pconf = paramconfig;
create_filter = 1;
run_ff = 1;
load_ffdata = 0;
stack_swrLFP = 1;
load_swrLFPstack = 0;
% run cmwc
make_ASig = 1;
load_rawpwr = 0;
load_phase = 0;
% create condition design mat
make_expvarCat = 1;
load_expvarCat = 0;
makeNoiseEvents = 1;
loadNoiseEvents = 0;
% compute per condition
make_expvarCatMeanPwr = 1;
load_expvarCatMeanPwr = 0;
% combine per area
combineArea = 0;
% plot
plot_expvarCatMeanPwr = 0;
plot_ByArea = 0;
pausefigs = 0;
savefigs = 1;

%% 
Fp.animals = {'D10','D13'}; %, 'JZ4'}; %;{'D10'}; %
Fp.filtfunction = 'dfa_riptriglfp';
% Fp.add_params = {'sleepwtrackdays', 'excludeNoise', '<4cm/s', 'wavelets4-300Hz'};
Fp.params = {'referenced','wtrackdays', 'excludePriorFirstWell','<4cm/s', ...
    '4-300Hz_focusSWR','excludeAfterLastWell', Fp.filtfunction};
Fp = load_filter_params(Fp);
wp = getWaveParams(Fp.waveSet);
rs = ''; % ripstate
area = {{'ca1', 'd'}, {'mec', 'deep'}, {'mec', 'supf'}};
%%
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', ...
        Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F), 'un', 1);
    F = runfilter(F);
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.env '_' Fp.eventSourceArea Fp.eventtype]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', ['_' Fp.env '_' Fp.eventSourceArea Fp.eventtype]);
end
%% vectorize swr-trig lfp
if stack_swrLFP
    lfpstack = stack_riptriglfp(F, Fp);
end
if load_swrLFPstack
    lfpstack = load_data(Fp.paths.resultsDirectory, ...
        ['riptriglfpstack_',Fp.env], Fp.animals);
end
%% make rawpwr (all trials) [ntrode time rip freq]
if make_ASig
    [rawpwr, ~] = computeAnalyticSignal(lfpstack, 'waveSet', Fp.waveSet, 'saveOutput',1, ...
        'lfptype', Fp.uselfptype, 'env', Fp.env);
end
if load_rawpwr
    rawpwr = load_data(sprintf('%s/analyticSignal/', pconf.andef{2}), ...
        sprintf('AS_waveSet-%s_%s_%s_power', wp.waveSet, Fp.uselfptype, Fp.env), ...
        Fp.animals);
end
%% make design mat to slice the rawpwr trials
if make_expvarCat
    outdir = 'expvarCatBurst';
    expvarCat = makeExpvarCatDesignMat(lfpstack, 'outdir', outdir, 'expvars', {'all', ...
        'lickbouts', 'nolickbouts'}, 'lfptype', Fp.uselfptype);
end
if load_expvarCat
    outdir = 'expvarCatBurst';
    outpath = [pconf.andef{2},outdir,'/'];
    expvarCat = load_data(outpath, [outdir,'_',Fp.env], Fp.animals);
end
%% exclude noise
if makeNoiseEvents
    noiseEvents = make_noiseEvents(lfpstack);
end
if loadNoiseEvents
    noiseEvents = load_data('filterframework','noiseEvents', Fp.animals, 'animpos', 0);
end
%% get mean power per condition
if make_expvarCatMeanPwr % :expvarCat @mean /ntTF $time
    expvarCatMeanPwr = getPower(expvarCat, rawpwr, Fp, 'run_perm', 0, 'noiseEvents',...
        noiseEvents);
end
if load_expvarCatMeanPwr
    outdir = 'expvarCatMeanPwr';
    outpath = [pconf.andef{2},outdir,'/'];
    expvarCatMeanPwr = load_data(outpath, [outdir,'_', Fp.env] ,Fp.animals);
end
%% mean per area per condition
if combineArea
    for ani = 1:length(expvarCatMeanPwr) % for each animal
        animal = expvarCatMeanPwr(ani).animal;
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        for ia = 1:length(area)
            ntsInArea = evaluatefilter(ntinfo,...
                sprintf('isequal($area,''%s'') && isequal($subarea,''%s'')', area{ia}{1}, ...
                area{ia}{2}));
            ntsA = unique(ntsInArea(:,3));
%             ntsAIdx = %find(1&sort(ntsA));
            for iv = 1:length(expvarCatMeanPwr(ani).expvars)
                areaData = expvarCatMeanPwr(ani).meandbpower{iv}.pwr_mean_db(ntsA,:,:);
                meanPwrArea{ia}{iv} = squeeze(nanmean(areaData,1))';
            end
        end
    end
end

%% swr example lfp trace over ca1 area all heatmap


%% plot per area
if plot_ByArea
    figname = 'expvarCatMeanPwrByArea';
    for ani = 1:length(expvarCatMeanPwr) % for each animal
        animal = expvarCatMeanPwr(ani).animal;
        
        for iv = 1:length(expvarCatMeanPwr(ani).expvars)
            
            Pp=load_plotting_params({'defaults','areaTFspect'}, 'pausefigs', pausefigs, ...
                    'savefigs', savefigs);
            f = 0;
            for ia = 1:length(area)
                f = f+1;
                sf{ia} = subaxis(1,length(area),f);
                frequency = round(expvarCatMeanPwr(ani).wp.frex);
                time = wp.win(1):1/(wp.srate/wp.dsamp):wp.win(2);
                trim = knnsearch(time',Pp.win(1)):knnsearch(time', Pp.win(2));
                ttime = time(trim);
                tdata = meanPwrArea{ia}{iv}(:,trim);
                contourf(ttime,frequency,tdata,40,'linecolor','none')
                set(gca,'ydir','normal','yscale','log');
                colormap(Pp.usecolormap)
                caxis(sf{ia}, 'auto')
                cax = caxis;
                caxis([0 cax(2)])
                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', Pp.tickFsize)
                hold on
                title(sprintf('%s %s',area{ia}{1},area{ia}{2}))
                yl = ylim;
                xl = xlim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                colorbar
            end
            ylabel(sf{1},'frequency Hz')
            xlabel(sf{1}, 'time s')
            %%
            stit = sprintf('%s %s %s %s %s',figname, expvarCatMeanPwr(ani).expvars{iv}, animal, ...
                Fp.env, Fp.waveSet);
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                save_figure([pconf.andef{4} figname], stit);
            end
        end

    end
    end

%% plot MeanPower varCat TFzmap /nt
if plot_expvarCatMeanPwr
    figname = 'expvarCatMeanPwr';
    for ani = 1:length(expvarCatMeanPwr) % for each animal
        Pp=load_plotting_params({'defaults','powerTFmap'}, 'pausefigs', pausefigs, ...
            'savefigs', savefigs);
        animal = expvarCatMeanPwr(ani).animal;
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo, 'area');
        matidx = unique(den.index(:,3));
        matidx = circshift(matidx,-1);
        anidx = find(strcmp({expvarCatMeanPwr.animal}, animal));
        evanidx = find(strcmp({expvarCatMeanPwr.animal}, animal));
        % exclude invalid tets
        invalidtets = evaluatefilter(ntinfo, 'isequal($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        if isempty(rs)
            rs = expvarCatMeanPwr(evanidx).expvars;
        end
        for i = 1:length(rs);
            iv = find(strcmp(rs{i}, expvarCatMeanPwr(evanidx).expvars));
            numcols = 8;
            numrows = 4; %ceil(length(ntrodes) / numcols);
            for nti = 1:length(ntrodes)
                nt = ntrodes(nti);
                if ismember(nt, invalidtets)
                    continue
                end
                area = ntinfo{1}{1}{nt}.area;
                subarea = ntinfo{1}{1}{nt}.subarea;
                cann = ntinfo{1}{1}{nt}.cannula;
                ntxy = ntinfo{1}{1}{nt}.ntxy;
                if cann == 'ca1'
                    ntxy(1) = ntxy(1)+4; % offset ca1 to the right
                end
                ntp = (ntxy(2)-1)*8+ntxy(1);
                sf = subaxis(numrows,numcols,ntp);
                if isnumeric(subarea)
                    subarea = num2str(subarea);
                end
                ntidx = find(matidx == nt);
                idata2plot = squeeze(...
                    expvarCatMeanPwr(anidx).meandbpower{iv}.pwr_mean_db(ntidx,:,:))';
                idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                    'dsamp', expvarCatMeanPwr(anidx).wp.dsamp);
                time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                contourf(sf, time, wp.frex, idata2plot, Pp.contourRes, ...
                    'linecolor','none');
                set(gca,'ydir','normal','yscale','log');
                colormap(hot); %Pp.usecolormap)
                caxis(sf, 'auto')
                hold on;
                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %%
            %         allAxesInFigure = findall(gcf,'type','axes');
            %         linkaxes(allAxesInFigure, 'xy');
            stit = sprintf('%s %s %s %s %s',figname, expvarCatMeanPwr(anidx).expvars{iv}, animal, ...
                Fp.env, Fp.waveSet);
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                save_figure([pconf.andef{4} figname], stit);
            end
        end
    end
end
