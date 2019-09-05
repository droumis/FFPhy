Fp.animals = {'JZ3'}; %{'D10','D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'JZ4'}; %;{'D10'}; %
Fp.filtfunction = 'dfa_riptriglfp';
% Fp.add_params = {'sleepwtrackdays', 'excludeNoise', '<4cm/s', 'wavelets4-300Hz'};
Fp.add_params = {'referenced','wtrackdays','excludeNoise','excludePriorFirstWell','<4cm/s', ...
    'wavelets4-300Hz','excludeAfterLastWell'};
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
pconf = paramconfig;

% LFP
make_swrLFP = 0;
save_swrLFP = make_swrLFP;
load_swrLFP = 0;
stack_swrLFP = 0; %make_swrLFP;
load_swrLFPstack = 0;
makenoiseEvents = 1;
% make.load Power, Phase
make_powerPhase = 0;
load_rawpwr = 0;
load_phase = 0;
% plot Power, Phase strips
plot_LFPstrips = 0;
plot_pwrStrips_timeXrip = 0;
plot_phaseStrips_timeXrip = 0;

% make.load Design Matrices: (expvarCat.expvarCatDiff.expvarCont.swrvarCont.tfbvarCont)
makeDesMats = 0;
make_expvarCat = makeDesMats;
make_expvarCont = makeDesMats;
make_swrvarCont = makeDesMats;
make_tfbvarCont = makeDesMats;
% load design matrices
loadDesMats = 0;
load_expvarCat = loadDesMats;
load_expvarCont = loadDesMats;
load_swrvarCont = loadDesMats;
load_tfbvarCont = loadDesMats;

% make Results: Mean,Corr,ITPC / covariate
make_pwr_results = 0;
make_phase_results = 0;
make_expvarCatMeanPwr = make_pwr_results;
make_expvarCatMeanPwrDiff = make_pwr_results;
make_varContPwrCorr = make_pwr_results; % expvarCont, swrvarCont tfbvarCont
make_expvarCatITPC = make_phase_results;
make_expvarCatITPCDiff = make_phase_results;
make_varContPhaseCorr = make_phase_results;
% load results
load_pwr_results = 0;
load_phase_results = 0;
load_expvarCatMeanPwr = load_pwr_results;
load_expvarCatMeanPwrDiff = load_pwr_results;
load_varContPwrCorr = load_pwr_results; % expvarCont, swrvarCont tfbvarCont
load_expvarCatITPC = load_phase_results;
load_expvarCatITPCDiff = load_phase_results;
% plot results
plot_pwr_results = 0;
plot_phase_results = 0;
plot_expvarCatMeanPwr = 0; plot_pwr_results;
plot_expvarCatMeanPwrDiff = plot_pwr_results;
plot_varCont = plot_pwr_results; % expvarCont, swrvarCont tfbvarCont
plot_expvarCatITPC = plot_phase_results;
plot_expvarCatITPCDiff = plot_phase_results;

% plot stats
plot_ContFit = 0; % fitLM
plot_CatDiffBars = 0; % KS
plot_combined_animalsareas = 0;
plot_CombCorr = 0;
plot_CombDiff = 0;

% plot options
usez = 1;
pauseb4supertit = 0;
pausefigs = 0;
savefigs = 1;

%% ripLFP [ntrode time rip]
if make_swrLFP
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', ...
        Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    F = runfilter(F); for d = 1:length(F); F(d).datafilter_params = Fp; end; end
if save_swrLFP
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s', Fp.epochEnvironment)); end
if load_swrLFP
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s', Fp.epochEnvironment)); end
% stack.loadstack riptrigLFP [nt times rip]
if stack_swrLFP; lfpstack = stack_riptriglfp(F, Fp); clear F; end
if load_swrLFPstack; lfpstack = load_data(Fp.paths.resultsDirectory, ...
        ['riptriglfpstack_',Fp.epochEnvironment], Fp.animals); end
%% Noise Events
if makenoiseEvents; noiseEvents = make_noiseEvents(lfpstack); else
    noiseEvents=load_data('filterframework','noiseEvents', Fp.animals, 'animpos', 0); end
%% Power, Phase [ntrode time rip freq]
if make_powerPhase; computeAnalyticSignal(lfpstack,'waveSet',Fp.waveSet,'saveOutput',1, ...
        'uselfptype', Fp.uselfptype, 'epochEnvironment', Fp.epochEnvironment); end
if load_rawpwr; rawpwr = load_data(sprintf('%s/analyticSignal/',pconf.andef{2}), ...
        sprintf('AS_waveSet-%s_%s_%s_power',Fp.wp.waveSet,Fp.uselfptype,Fp.epochEnvironment), ...
        Fp.animals); end
if load_phase; phase = load_data(sprintf('%s/analyticSignal/', pconf.andef{2}), ...
        sprintf('AS_waveSet-%s_%s_%s_phase',Fp.wp.waveSet,Fp.uselfptype,Fp.epochEnvironment), ...
        Fp.animals);end
%% Design Matrices [rip var (nt) (ntB)]
% :expvarCat @mean @ITPC /rip $time [rip var] ...
% {all, outbound, inbound, rewarded, unrewarded, proximalWell, distalWell}
if make_expvarCat; outdir = 'expvarCat';
    expvarCat = makeExpvarCatDesignMat(lfpstack, 'outdir', outdir); end
if load_expvarCat; outdir = 'expvarCat'; outpath = [pconf.andef{2},outdir,'/'];
    expvarCat = load_data(outpath, [outdir,'_',Fp.epochEnvironment], Fp.animals);end
% :expvarCont @corr /rip $swr [rip var] {speed, performance, learningrate}
if make_expvarCont; outdir = 'expvarCont';
    expvarCont = makeExpvarContDesignMat(lfpstack, Fp, 'outdir', outdir); end
if load_expvarCont; outdir = 'expvarCont'; outpath = [pconf.andef{2},outdir,'/'];
    expvarCont = load_data(outpath, [outdir,'_',Fp.epochEnvironment], Fp.animals); end
% :swrvarCont @corr /rip $swr [rip var] {std, duration}
if make_swrvarCont; outdir = 'swrvarCont';
    swrvarCont = makeSWRDesignMatrix(lfpstack, Fp, 'outdir', outdir); end
if load_swrvarCont; outdir = 'swrvarCont'; outpath = [pconf.andef{2},outdir,'/'];
    swrvarCont = load_data(outpath, [outdir,'_',Fp.epochEnvironment], Fp.animals); end
% :tfbvarCont @corr /rip $swr [rip var nt] {Ripple FastGamma Theta}
if make_tfbvarCont; outdir = 'tfbvarCont';
    tfbvarCont = makeTFBoxDesignMat(rawpwr,'epEnv',Fp.epochEnvironment,'outdir',outdir);end
if load_tfbvarCont; outdir = 'tfbvarCont'; outpath = [pconf.andef{2},outdir,'/'];
    tfbvarCont = load_data(outpath, [outdir,'_',Fp.epochEnvironment], Fp.animals); end

%% Mean, ITPC, Corr
% don't run perm test for means bc time shuffling doesn't work for tonic signal
if make_expvarCatMeanPwr % :expvarCat @mean /ntTF $time
    expvarCatMeanPwr = getPower(expvarCat, rawpwr,Fp, 'run_perm', 0, 'noiseEvents',...
        noiseEvents); end
if make_expvarCatMeanPwrDiff % :expvarCatDiff @dmean /ntTF $var
    expvarCatMeanPwrDiff = getPowerDiff(expvarCat,rawpwr,Fp, 'run_perm', 0, 'noiseEvents', ...
        noiseEvents, 'tfboxes', tfbvarCont); end

if make_varContPwrCorr % @corr (expvarCont tfbvarCont swrvarCont) /ntTF $swr
    expvarContPwrCorr = runDesignDataRegression(expvarCont,rawpwr,Fp,'outdir','expvarCont', ...
        'noiseEvents', noiseEvents, 'tfboxes', tfbvarCont);
    swrvarContPwrCorr = runDesignDataRegression(swrvarCont,rawpwr,Fp,'outdir','swrvarCont', ...
        'noiseEvents', noiseEvents, 'tfboxes', tfbvarCont);
    tfbvarContPwrCorr = runDesignDataRegression(tfbvarCont,rawpwr,Fp,'outdir','tfbvarCont', ...
        'noiseEvents', noiseEvents, 'tfboxes', tfbvarCont); end

if make_expvarCatITPC  % :expvarCat @itpc /ntTF $time
    expvarCatITPC = getITPC(expvarCat, phase, Fp,'run_perm', 0,'noiseEvents', noiseEvents); end
if make_expvarCatITPCDiff % :expvarCatDiff @ditpc /ntTF $var
    expvarCatITPCDiff = makeExpvarCatITPCDiff(expvarCat, phase, Fp, 'run_perm', 1, ...
        'noiseEvents', noiseEvents); end
if make_varContPhaseCorr % @corr (expvarCont tfbvarCont swrvarCont) /ntTF $swr
    expvarContPwrCorr = runDesignDataRegression(expvarCont,rawpwr,Fp,'outdir','expvarCont', ...
        'noiseEvents', noiseEvents, 'tfboxes', tfbvarCont);
    swrvarContPwrCorr = runDesignDataRegression(swrvarCont,rawpwr,Fp,'outdir','swrvarCont', ...
        'noiseEvents', noiseEvents, 'tfboxes', tfbvarCont);
    tfbvarContPwrCorr = runDesignDataRegression(tfbvarCont,rawpwr,Fp,'outdir','tfbvarCont', ...
        'noiseEvents', noiseEvents, 'tfboxes', tfbvarCont); end

%% loading Mean, ITPC, Corr
if load_expvarCatMeanPwr; outdir = 'expvarCatMeanPwr'; outpath = [pconf.andef{2},outdir,'/'];
    expvarCatMeanPwr = load_data(outpath, [outdir,'_', Fp.epochEnvironment] ,Fp.animals); end
if load_expvarCatMeanPwrDiff; outdir = 'expvarCatMeanPwrDiff'; outpath = [pconf.andef{2},outdir,'/'];
    expvarCatMeanPwrDiff = load_data(outpath, [outdir,'_',Fp.epochEnvironment], Fp.animals);end
if load_varContPwrCorr
    outdir = 'expvarCont'; outfile='expvarContCorr'; outpath=[pconf.andef{2},outdir,'/'];
    expvarContPwrCorr=load_data(outpath,[outfile,'_', Fp.epochEnvironment],Fp.animals);
    outdir = 'swrvarCont'; outfile='swrvarContCorr'; outpath=[pconf.andef{2},outdir,'/'];
    swrvarContPwrCorr=load_data(outpath,[outfile,'_', Fp.epochEnvironment],Fp.animals);
    outdir = 'tfbvarCont'; outfile='tfbvarContCorr'; outpath=[pconf.andef{2},outdir,'/'];
    tfbvarContPwrCorr=load_data(outpath,[outfile,'_', Fp.epochEnvironment],Fp.animals); end
if load_expvarCatITPC; outdir = 'expvarCatITPC'; outpath = [pconf.andef{2},outdir,'/'];
    expvarCatITPC = load_data(outpath, [outdir,'_', Fp.epochEnvironment] ,Fp.animals); end
if load_expvarCatITPCDiff; outdir = 'expvarCatITPCDiff'; outpath = [pconf.andef{2},outdir,'/'];
    expvarCatITPCDiff = load_data(outpath, [outdir,'_', Fp.epochEnvironment] ,Fp.animals); end

%% plot LFP (1-400Hz) heatRast strips /nt
if plot_LFPstrips; Pp=load_plotting_params({'defaults','powerheatRast'});
    for ian = 1:numel(Fp.animals)
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        anidx = find(strcmp({lfpstack.animal}, animal));
        ntrodes = lfpstack(anidx).ntrodes;
        dayep = [lfpstack(anidx).day lfpstack(anidx).epoch];
        evanidx = find(strcmp({expvarCat.animal}, animal));
        % exclude invalid tets
        invalidtets = evaluatefilter(ntinfo, 'isequal($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        % exclude invalid rips
        userips = ones(length(lfpstack(anidx).day),1);
        if ~isempty(noiseEvents)
            noiseanidx = find(strcmp({noiseEvents.animal}, animal));
            if ~isempty(noiseEvents(noiseanidx).events)
                invalidrips = find(ismember([lfpstack(anidx).day lfpstack(anidx).epoch ...
                    lfpstack(anidx).ripStartTime], ...
                    noiseEvents(noiseanidx).events, 'rows'));
                userips(invalidrips) = 0; end; end
        rs = 'onlywdays'; % ripstate
        istateidx = find(strcmp(rs, expvarCat(evanidx).expvars));
        include_rips = expvarCat(ian).dm(:,istateidx);
        include_rips = all([include_rips userips], 2);
        if savefigs && ~pausefigs; close all
            ifig = figure('Visible','off','units','normalized','position', Pp.position);
        else; ifig = figure('units','normalized','position',Pp.position); end
        set(gcf,'color','white')
        for nti = 1:length(ntrodes)
            ntrode = ntrodes(nti);
            area = ntinfo{1}{1}{ntrode}.area;
            if ismember(ntrode, invalidtets); continue; end
            sf = subaxis(2,ceil(max(ntrodes)/2),nti, 'SpacingVert', Pp.SpVt, ...
                'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
            h = zoom; h.Motion = 'vertical'; h.Enable = 'on';
            il = find(strcmp(lfpstack(anidx).lfptypes, Fp.uselfptype));
            excld_stack = squeeze(double(lfpstack(anidx).data{il}...
                (nti,:,find(include_rips))))';
            idata2plot = trim2win(excld_stack, Fp.srate, Pp.pwin);
            ptime = linspace(-Pp.pwin(1),Pp.pwin(2),size(idata2plot,2));
            z = idata2plot;
            if usez
                m = nanmean(idata2plot,2);
                s = nanstd(idata2plot, [], 2);
                z = (idata2plot-m)./s;
            end
            imagesc(ptime, 1:size(z,1), z)
            colormap(Pp.cmap); caxis(sf, 'auto')
            % epoch and day bounds
            exdayep = dayep(find(include_rips),:); de = unique(exdayep, 'rows');
            daybounds = find(diff(exdayep(:,1))); daybounds = [1;daybounds];
            epbounds = find(abs(diff(exdayep(:,2))));
            line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color',[0 0 0])
            line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [.9 .9 .9])
            line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
            title(sprintf('nt%d %s', ntrode, area), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                'FontName', Pp.FontNm)
            xlabel('time', 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
            ylabel('ripnum','FontSize',Pp.FontS, 'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
            if nti ~= 1
                xlabel('')
                ylabel('')
                set(gca, 'ytick', []);
                set(gca, 'xtick', []);
            else
                days = num2str(unique(de(:,1)));
                for id = 1:length(days)
                    text(-1.2, daybounds(id), days(id), 'FontSize', 8, 'Color','r');
                end
            end
            allAxesInFigure = findall(ifig,'type','axes');
            linkaxes(allAxesInFigure, 'y');
        end
        %% super
        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
        if usez
            sprtit = sprintf('zLFPHeatRast %s %s %s %s', animal, Fp.uselfptype, rs, ...
                Fp.epochEnvironment);
        else
            sprtit = sprintf('LFPHeatRaster %s %s %s %s', animal, Fp.uselfptype, rs, ...
                Fp.epochEnvironment); end
        iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
        set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center','FontSize', 12);
        spylabel = text(.02, .5, sprintf('ripnum'), 'Parent', sprtitleax, 'Units', ...
            'normalized', 'Rotation', 90);
        set(spylabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center', 'FontSize', 12);
        spxlabel = text(.5, .02, sprintf('time'), 'Parent', sprtitleax, 'Units', ...
            'normalized');
        set(spxlabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center', 'FontSize', 12);
        %% ---- pause, save figs ----
        if pausefigs; pause; end; if savefigs; save_figure([pconf.andef{4},'/LFPHeatRaster/'],sprtit); end
    end; end
%% plot power strips heatRast /nt /freq
if plot_pwrStrips_timeXrip; Pp=load_plotting_params({'defaults','powerheatRast'});
    for ian = 1:numel(Fp.animals)
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        anidx = find(strcmp({rawpwr.animal}, animal));
        ntrodes = rawpwr(anidx).ntrode;
        dayep = [rawpwr(anidx).day rawpwr(anidx).epoch];
        evanidx = find(strcmp({expvarCat.animal}, animal));
        % exclude invalid tets
        invalidtets = evaluatefilter(ntinfo, 'isequal($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        % exclude invalid rips
        userips = ones(length(rawpwr(anidx).day),1);
        if ~isempty(noiseEvents)
            noiseanidx = find(strcmp({noiseEvents.animal}, animal));
            if ~isempty(noiseEvents(noiseanidx).events)
                invalidrips = find(ismember([rawpwr(anidx).day rawpwr(anidx).epoch ...
                    rawpwr(anidx).ripStartTime], ...
                    noiseEvents(noiseanidx).events, 'rows'));
                userips(invalidrips) = 0;
            end
        end
        for fx = 1:length(Pp.plot_frex)
            frexidx = knnsearch(rawpwr.wp.frex', Pp.plot_frex(fx));
            knnfrex = round(rawpwr.wp.frex(frexidx));
            % plot per ripstate
            for rs = 1 %:length(Fp.use)
                %         use_filts = find(any(cell2mat(cellfun(@(x) strcmp(x, lfpstack(anidx).filterfields), ...
                %             use_filters, 'un', 0)), 2));
                %                 istate = Fp.useripstates{rs};
                istateidx = find(strcmp('onlywdays', expvarCat(evanidx).expvars));
                include_rips = expvarCat(ian).dm(:,istateidx);
                %                 include_rips = all([include_rips userips], 2);
                %% ---- init fig----
                if savefigs && ~pausefigs
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                
                for nti = 1:length(ntrodes)
                    ntrode = ntrodes(nti);
                    if ismember(ntrode, invalidtets)
                        continue
                    end
                    %                     2,ceil(max(ntrodes)/2),
                    sf = subaxis(2,ceil(max(ntrodes)/2),nti, 'SpacingVert', Pp.SpVt, ...
                        'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                        Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    h = zoom;
                    h.Motion = 'vertical';
                    h.Enable = 'on';
                    exdayep = dayep(find(include_rips),:);
                    de = unique(exdayep, 'rows');
                    daybounds = find(diff(exdayep(:,1)));
                    daybounds = [1;daybounds];
                    epbounds = find(abs(diff(exdayep(:,2))));
                    
                    excld_stack=squeeze(double(rawpwr(anidx).pwr(nti,:,find(include_rips),...
                        frexidx)))';
                    idata2plot = trim2win(excld_stack, Fp.srate, Pp.pwin, ...
                        'dsamp', rawpwr(anidx).wp.dsamp);
                    idata2plot(~userips,:) = 0;
                    % mididx = ceil(size(excld_stack,2)/2); % right now assumes center is rip start
                    ptime = linspace(-Pp.pwin(1),Pp.pwin(2),size(idata2plot,2));
                    z = idata2plot;
                    if usez
                        m = nanmean(idata2plot,2);
                        s = nanstd(idata2plot, [], 2);
                        z = (idata2plot-m)./s;
                    end
                    imagesc(ptime, 1:size(z,1), z)
                    colormap(Pp.cmap)
                    caxis(sf, 'auto')
                    
                    line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color',[0 0 0])
                    line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [.9 .9 .9])
                    line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
                    title(sprintf('%d', ntrode), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                        'FontName', Pp.FontNm)
                    % caxis([-1 1])
                    xlabel('time s', 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', ...
                        Pp.FontNm)
                    ylabel('ripnum (day-w epoch-b)','FontSize',Pp.FontS, ...
                        'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                    if nti ~= 1
                        xlabel('')
                        ylabel('')
                        set(gca, 'ytick', []);
                        set(gca, 'xtick', []);
                    else
                        days = num2str(unique(de(:,1)));
                        for id = 1:length(days)
                            text(-1.2, daybounds(id), days(id), 'FontSize', 8, 'Color','r');
                        end
                    end
                    allAxesInFigure = findall(ifig,'type','axes');
                    linkaxes(allAxesInFigure, 'y');
                end
                
                %% super
                if pauseb4supertit
                    pause
                end
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                if usez
                    sprtit = sprintf('Zmeandb %s all%d %s %dHz %s', animal, rs, Fp.uselfptype, ...
                        knnfrex, Fp.epochEnvironment);
                else
                    sprtit = sprintf('meandb %s all%d %s %dHz %s', animal, rs, Fp.uselfptype, ...
                        knnfrex, Fp.epochEnvironment);
                end
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                spylabel = text(.02, .5, sprintf('ripnum'), 'Parent', sprtitleax, 'Units', ...
                    'normalized', 'Rotation', 90);
                set(spylabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center', 'FontSize', 12);
                
                spxlabel = text(.5, .02, sprintf('time'), 'Parent', sprtitleax, 'Units', ...
                    'normalized');
                set(spxlabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center', 'FontSize', 12);
                %% ---- pause, save figs ----
                if pausefigs; pause; end; if savefigs;
                    save_figure([pconf.andef{4},'/powerHeatRast/'],sprtit); end
            end; end; end; end
%% plot phase strips heatRast /nt /freq
if plot_phaseStrips_timeXrip; Pp=load_plotting_params({'defaults','powerheatRast'});
    for ian = 1:numel(Fp.animals)
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        anidx = find(strcmp({rawpwr.animal}, animal));
        ntrodes = lfpstack(anidx).ntrodes;
        dayep = [lfpstack(anidx).day lfpstack(anidx).epoch];
        
        % exclude invalid tets
        invalidtets = evaluatefilter(ntinfo, 'isequal($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        
        for fx = 1:length(Pp.plot_frex)
            frexidx = knnsearch(rawpwr.wp.frex', Pp.plot_frex(fx));
            knnfrex = round(rawpwr.wp.frex(frexidx));
            % plot per ripstate
            for rs = 1 %:length(Fp.use)
                %         use_filts = find(any(cell2mat(cellfun(@(x) strcmp(x, lfpstack(anidx).filterfields), ...
                %             use_filters, 'un', 0)), 2));
                istate = Fp.useripstates{rs};
                istateidx = find(strcmp(istate, expvarCat.expvars));
                include_rips = expvarCat.statesets(:,istateidx);
                %% ---- init fig----
                if savefigs && ~pausefigs
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                
                for nti = 1:length(ntrodes)
                    ntrode = ntrodes(nti);
                    if ismember(ntrode, invalidtets)
                        continue
                    end
                    sf = subaxis(2,ceil(max(ntrodes)/2), nti, 'SpacingVert', Pp.SpVt, ...
                        'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                        Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    
                    exdayep = dayep(find(include_rips),:);
                    de = unique(exdayep, 'rows');
                    daybounds = find(diff(exdayep(:,1)));
                    epbounds = find(abs(diff(exdayep(:,2))));
                    %                     [d,i] = max(squeeze(double(pwr(anidx).pwr(11, :, find(include_rips),frexidx))));
                    %                     [j,k] = max(d);
                    %                     include_rips(k) = 0;
                    %
                    excld_stack = squeeze(double(rawpwr(anidx).pwr(nti, :, find(include_rips),frexidx)))';
                    idata2plot = trim2win(excld_stack, Fp.srate, Pp.pwin, ...
                        'dsamp', rawpwr(anidx).wp.dsamp);
                    %mididx = ceil(size(excld_stack,2)/2); % right now assumes center is rip start
                    ptime = linspace(-Pp.pwin(1),Pp.pwin(2),size(idata2plot,2));
                    z = idata2plot;
                    %                     m = nanmean(idata2plot,2);
                    %                     s = nanstd(idata2plot, [], 2);
                    %                     z = (idata2plot-m)./s;
                    imagesc(ptime, 1:size(z,1), z)
                    colormap(parula)
                    caxis(sf, 'auto')
                    
                    line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color',[.9 .9 .9])
                    line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
                    line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
                    title(sprintf('%d', ntrode), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                        'FontName', Pp.FontNm)
                    %                         caxis([-1 1])
                    xlabel('time s', 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', ...
                        Pp.FontNm)
                    ylabel('ripnum (day-b epoch-w)','FontSize',Pp.FontS, ...
                        'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                    if nti ~= 1
                        xlabel('')
                        ylabel('')
                        set(gca, 'ytick', []);
                        set(gca, 'xtick', []);
                    end
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('meandb %s %s %s %dHz %s', animal, istate, Fp.uselfptype, ...
                    knnfrex, Fp.epochEnvironment);
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                spylabel = text(.02, .5, sprintf('ripnum'), 'Parent', sprtitleax, 'Units', ...
                    'normalized', 'Rotation', 90);
                set(spylabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center', 'FontSize', 12);
                
                spxlabel = text(.5, .02, sprintf('time'), 'Parent', sprtitleax, 'Units', ...
                    'normalized');
                set(spxlabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center', 'FontSize', 12);
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    save_figure([pconf.andef{4},'/powerHeatRast/'],sprtit)
                    close all
                end
            end; end; end; end
%% plot MeanPower varCat TFzmap /nt
if plot_expvarCatMeanPwr; Pp = load_plotting_params({'defaults', 'powerTFmap'});
    for ian = 1:length(expvarCatMeanPwr) % for each animal
        animal = expvarCatMeanPwr(ian).animal;
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
        rs = 'onlywdays'; % ripstate
        iv = find(strcmp(rs, expvarCatMeanPwr(evanidx).expvars));
        if savefigs && ~pausefigs
            close all
            ifig =figure('Visible','off','units','normalized','position', ...
                Pp.position);
        else
            ifig = figure('units','normalized','position',Pp.position);
        end
        set(gcf,'color','white')
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
            sf = subaxis(numrows,numcols,ntp, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                Pp.MgTp, 'MarginBottom', Pp.MgBm);
            if isnumeric(subarea)
                subarea = num2str(subarea);
            end
            ntidx = find(matidx == nt);
            idata2plot = squeeze(...
                expvarCatMeanPwr(anidx).meandbpower{iv}.pwr_mean_db(ntidx,:,:))';
            idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                'dsamp', expvarCatMeanPwr(anidx).wp.dsamp);
            time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
            contourf(sf, time, Fp.wp.frex, idata2plot, Pp.contourRes, ...
                'linecolor','none');
            set(gca,'ydir','normal','yscale','log');
            colormap(jet); %Pp.usecolormap)
            caxis(sf, 'auto')
            hold on
%             % thresholded single pix zmask
%             if ~isempty(fieldnames(expvarCatMeanPwr(anidx).meandbpower{iv}.permt))
%                 zmask2plot = squeeze(expvarCatMeanPwr(anidx).meandbpower{iv}.permt.threshmean(ntidx,:,:))';
%                 zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', Fp.wp.dsamp);
%                 try
%                     [~,h] = contour(sf, time, Fp.wp.frex, logical(zmask2plot), 1);
%                     h.LineColor = 'black';
%                 catch
%                     fprintf('invalid zmask\n')
%                 end
%             end
            hold on;
            ytickskip = 2:4:Fp.wp.numfrex;
            set(gca,'ytick', round(Fp.wp.frex(ytickskip)), 'FontSize', 8)
            title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
            yl = ylim;
            line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
        end
        %% super
        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
        sprtit = sprintf('meanpwr %s %s %s', expvarCatMeanPwr(anidx).expvars{iv}, animal, ...
            Fp.epochEnvironment);
        iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
        set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center','FontSize', 16);        
        %% ---- pause, save figs ----
        if pausefigs;  pause; end
        if savefigs; save_figure(sprintf('%s/expvarCatMeanPwr/',pconf.andef{4}), sprtit); end
    end; end
%% plot meanpowerDiff varCat TFzmap /nt
if plot_expvarCatMeanPwrDiff; Pp = load_plotting_params({'defaults', 'powerTFmap'});
    for ian = 1:length(expvarCatMeanPwrDiff) % for each animal
        animal = expvarCatMeanPwrDiff(ian).animal;
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo, 'area');
        matidx = unique(den.index(:,3));
        matidx = circshift(matidx,-1);
        anidx = find(strcmp({expvarCatMeanPwrDiff.animal}, animal));
        for iv = 1:length(expvarCatMeanPwrDiff(anidx).expvars)
            if savefigs && ~pausefigs; close all
                ifig =figure('Visible','off','units','normalized','position', Pp.position);
            else; ifig = figure('units','normalized','position',Pp.position); end
            set(gcf,'color','white')
            numcols = 8; numrows = 4; %ceil(length(ntrodes) / numcols);
            for nti = 1:length(ntrodes)
                nt = ntrodes(nti); area = ntinfo{1}{1}{nt}.area; 
                subarea = ntinfo{1}{1}{nt}.subarea; cann = ntinfo{1}{1}{nt}.cannula;
                ntxy = ntinfo{1}{1}{nt}.ntxy;
                if cann == 'ca1'; ntxy(1) = ntxy(1)+4; end
                ntp = (ntxy(2)-1)*8+ntxy(1);
                sf = subaxis(numrows,numcols,ntp, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                    'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                    Pp.MgTp, 'MarginBottom', Pp.MgBm);
                if isnumeric(subarea); subarea = num2str(subarea); end
%                 ntidx = find(matidx == nt);
                idata2plot = squeeze(expvarCatMeanPwrDiff(anidx).meandbpowerDiff{iv}.pwr_meandb_diff(nti,:,:))';
                idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                    'dsamp', expvarCatMeanPwrDiff(anidx).wp.dsamp);
                time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                contourf(sf, time, Fp.wp.frex, idata2plot, Pp.contourRes, ...
                    'linecolor','none');
                set(gca,'ydir','normal','yscale','log'); colormap(Pp.usecolormap)
                caxis(sf, 'auto'); hold on
                ytickskip = 2:4:Fp.wp.numfrex;
                set(gca,'ytick', round(Fp.wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                yl = ylim; xl = xlim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);                
                tfbani = find(strcmp({tfbvarCont.animal}, animal));
                for itfb = 1:length(expvarCatMeanPwrDiff(anidx).tfb)
                    tfb = tfbvarCont(tfbani).expvars{itfb};
                    realfreq = tfbvarCont(tfbani).realfreq{itfb};
                    mAv = mean(expvarCatMeanPwrDiff(anidx).Avals{nti, itfb, iv});
                    mBv = mean(expvarCatMeanPwrDiff(anidx).Bvals{nti, itfb, iv});
                    P = expvarCatMeanPwrDiff(anidx).P(nti, itfb, iv);
                    ln = line([xl(1) xl(1)],[min(realfreq) max(realfreq)],'LineWidth',20);
                    if any(strcmp(tfb, 'Theta')) ||  any(strcmp(tfb, 'FastGamma'))
                        if P < .01 && mAv > mBv; ln.Color = [0 1 0];
                        elseif P < .01 && mAv < mBv; ln.Color = [1 0 1];
                        else; ln.Color = [.5 .5 .5]; end
                    else; if P < .01 && coef > 0; ln.Color = [1 0 0]; 
                        elseif P < .01 && coef < 0; ln.Color = [0 0 1];
                        else; ln.Color = [.5 .5 .5]; end; end; end
                            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('meanpwrdiff %s-%s %s %s', expvarCatMeanPwrDiff(anidx).expvars{iv}{1},...
                expvarCatMeanPwrDiff(anidx).expvars{iv}{2}, animal, Fp.epochEnvironment);
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 16);
            %% ---- pause, save figs ----
            if pausefigs; pause; end
            if savefigs; save_figure(sprintf('%s/expvarCatMeanPwrDiff/',pconf.andef{4}), sprtit)
            end
            %                 % thresholded zmask
%                 if ~isempty(fieldnames(expvarCatMeanPwrDiff(anidx).meandbpowerDiff{iv}.permt))
%                     zmask2plot = squeeze(expvarCatMeanPwrDiff(anidx).meandbpowerDiff{iv}.permt.threshmean(ntidx,:,:))';
%                     zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', Fp.wp.dsamp);
%                     try
%                         [~,h] = contour(sf, time, Fp.wp.frex, logical(zmask2plot), 1);
%                         h.LineColor = 'black';
%                     catch
%                         fprintf('invalid zmask\n')
%                     end
%                 end                
%                 v = [-.4 yl(1); .4 yl(1); 0 25]; f = [1 2 3];
%                 patch('Faces',f,'Vertices',v, 'EdgeColor','black','FaceColor','black','FaceAlpha', .3, 'LineWidth',1);
%                 usetfboxes = {'prePreSwrTheta', 'postPostSwrTheta', 'swrRipple', 'preSwrFastGamma', 'postSwrFastGamma'};
%                 for tfb = 1:length(usetfboxes)
%                     tfbani = find(strcmp({tfbvarCont.animal}, animal));
%                     ppthidx = find(strcmp(tfbvarCont(tfbani).expvars, usetfboxes{tfb}));
%                     if isempty(ppthidx)
%                         fprintf('could not find %s\n', usetfboxes{tfb});
%                         continue
%                     end
%                     ppthyidx = [min(tfbvarCont(tfbani).realfreq{ppthidx}); ...
%                         max(tfbvarCont(tfbani).realfreq{ppthidx})];
%                     ppthxidx = tfbvarCont(tfbani).time{ppthidx};
%                     xw = ppthxidx(2)-ppthxidx(1);
%                     yh = ppthyidx(2)-ppthyidx(1);
%                     cx = ppthxidx(1);
%                     cy = ppthyidx(1);
%                     rectangle('Position',[cx cy xw yh])
%                     hold on
%                 end
        end; end; end
%% corrpower varCont TFzmap /nt
if plot_varCont; Pp = load_plotting_params({'defaults', 'powerTFmap'});
    for ian = 1:numel(tfbvarContPwrCorr) % for each animal
        animal = tfbvarContPwrCorr(ian).animal;
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        conds = {'expvarContPwrCorr', 'swrvarContPwrCorr', 'tfbvarContPwrCorr'}; %,
        for icorr = 1:length(conds)
            PV = eval(conds{icorr});
            pvanidx = find(strcmp({PV.animal}, animal));
            for iv = 1:length(PV(pvanidx).expvars)
                if savefigs && ~pausefigs; close all
                    ifig =figure('Visible','off','units','normalized','position',Pp.position);
                else; ifig = figure('units','normalized','position',Pp.position); end
                set(gcf,'color','white')
                numcols = 8; numrows = 4; %ceil(length(ntrodes) / numcols);
                for nti = 1:length(ntrodes)
                    nt = ntrodes(nti);
                    area = ntinfo{1}{1}{nt}.area;
                    subarea = ntinfo{1}{1}{nt}.subarea;
                    cann = ntinfo{1}{1}{nt}.cannula;
                    ntxy = ntinfo{1}{1}{nt}.ntxy;
                    if cann == 'ca1'
                        ntxy(1) = ntxy(1)+4; % offset ca1 to the right
                    end
                    ntp = (ntxy(2)-1)*8+ntxy(1);
                    sf = subaxis(numrows,numcols,ntp, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                        Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    zmap = squeeze(PV(pvanidx).zmap(nti,:,:,iv));
                    zmapthresh = squeeze(PV(pvanidx).clusterZmapThresh(nti, :,:,iv));
                    contourf(PV(pvanidx).time,PV(pvanidx).frequency,zmap,40,'linecolor','none')
                    hold on
                    %                 [~,h] = contour(PV(pvanidx).time,PV(pvanidx).frequency,logical(zmapthresh),1);
                    %                 h.LineColor = 'black';
                    set(gca,'ydir','normal','yscale','log');
                    caxis(sf, 'auto')
                    colormap(Pp.usecolormap)
                    ytickskip = 2:4:Fp.wp.numfrex;
                    set(gca,'ytick', round(Fp.wp.frex(ytickskip)), 'FontSize', Pp.tickFsize)
                    title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',Pp.sfTitFsize,...
                        'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                    yl = ylim;
                    xl = xlim;
                    line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                    tfbani = find(strcmp({tfbvarCont.animal}, animal));
                    for itfb = 1:length(tfbvarCont(tfbani).expvars)
                        tfb = tfbvarCont(tfbani).expvars{itfb};
                        realfreq = tfbvarCont(tfbani).realfreq{itfb};
                        P = PV(pvanidx).P(nti, itfb, iv);
                        R = PV(pvanidx).R(nti, itfb, iv);
                        coef = PV(pvanidx).coef(nti, itfb, iv);
                        tfboxYidx = [min(realfreq); max(realfreq)];
                        tfboxXidx = tfbvarCont(tfbani).time{itfb};
                        xw = tfboxXidx(2)-tfboxXidx(1);
                        yh = tfboxYidx(2)-tfboxYidx(1);
                        cx = tfboxXidx(1);
                        cy = tfboxYidx(1);
                        rl = rectangle('Position',[cx cy xw yh]);
                        ln = line([xl(1) xl(1)],tfboxYidx,'LineWidth',10);
                        if any(strcmp(tfb, 'Theta')) ||  any(strcmp(tfb, 'FastGamma'))
                            if P < .01 && coef > 0
                                ln.Color = [0 1 0];
                            elseif P < .01 && coef < 0
                                ln.Color = [1 0 1];
                            else
                                ln.Color = [.5 .5 .5];
                            end
                        else
                            if P < .01 && coef > 0
                                ln.Color = [1 0 0];
                            elseif P < .01 && coef < 0
                                ln.Color = [0 0 1];
                            else
                                ln.Color = [.5 .5 .5];
                            end
                        end
                    end
                end
                %% super ax
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('Corr %s %s', animal, PV(ian).expvars{iv});
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', Pp.stitFsize);
                %% pause, save
                if pausefigs; pause; end
                if savefigs; save_figure(sprintf('%s/%s/',pconf.andef{4}, conds{icorr}), sprtit); end
            end; end; end; end
% ITPC varCat TFzmap /nt
if plot_expvarCatITPC; Pp = load_plotting_params({'defaults', 'powerTFmap'});
    animals = {expvarCatITPC.animal};
    for ian = 1:length(animals) % for each animal
        animal = animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo, 'area');
        matidx = unique(den.index(:,3));
        anidx = find(strcmp({expvarCatITPC.animal}, animal));
        for iv = 1:length(expvarCatITPC(anidx).expvars)
            if savefigs && ~pausefigs
                close all
                ifig =figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            numcols = 8;
            numrows = 4; %ceil(length(ntrodes) / numcols);
            for nti = 1:length(ntrodes)
                nt = ntrodes(nti);
                area = ntinfo{1}{1}{nt}.area;
                subarea = ntinfo{1}{1}{nt}.subarea;
                %                 try
                cann = ntinfo{1}{1}{nt}.cannula;
                %                 catch
                %                     fprintf('%d\n',nt)
                %                 end
                ntxy = ntinfo{1}{1}{nt}.ntxy;
                if cann == 'ca1'
                    ntxy(1) = ntxy(1)+4; % offset ca1 to the right
                end
                ntp = (ntxy(2)-1)*8+ntxy(1);
                sf = subaxis(numrows,numcols,ntp, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                    'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                    Pp.MgTp, 'MarginBottom', Pp.MgBm);
                
                if isnumeric(subarea)
                    subarea = num2str(subarea);
                end
                ntidx = find(matidx == nt);
                idata2plot = squeeze(expvarCatITPC(anidx).ITPC{iv}.ITPC_db(ntidx,:,:))';
                idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                    'dsamp', expvarCatITPC(anidx).wp.dsamp);
                time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                contourf(sf, time, Fp.wp.frex, idata2plot, 40, ...
                    'linecolor','none');
                set(gca,'ydir','normal','yscale','log');
                
                colormap(Pp.usecolormap)
                caxis(sf, 'auto')
                %                 colorbar
                
                hold on
                % thresholded single pix zmask
                if ~isempty(fieldnames(expvarCatITPC(anidx).ITPC{iv}.permt))
                    zmask2plot = squeeze(expvarCatITPC(anidx).ITPC{iv}.permt.threshmean(ntidx,:,:))';
                    zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', Fp.wp.dsamp);
                    try
                        [~,h] = contour(sf, time, Fp.wp.frex, logical(zmask2plot), 1);
                        h.LineColor = 'black';
                    catch
                        fprintf('invalid zmask\n')
                    end
                end
                hold on;
                ytickskip = 2:4:Fp.wp.numfrex;
                set(gca,'ytick', round(Fp.wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('itpc %s %s %s', expvarCatITPC(anidx).expvars{iv}, animal, ...
                Fp.epochEnvironment);
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 16);
            
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(sprintf('%s/expvarCatITPC/',pconf.andef{4}), sprtit)
                close all
            end
            close all;
            %                 end
        end; end; end
% ITPCDiff varCat TFzmap /nt
if plot_expvarCatITPCDiff; Pp = load_plotting_params({'defaults', 'powerTFmap'});
    animals = {expvarCatITPCDiff.animal};
    for ian = 1:length(animals) % for each animal
        animal = animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo, 'area');
        matidx = unique(den.index(:,3));
        anidx = find(strcmp({expvarCatITPCDiff.animal}, animal));
        for iv = 1:length(expvarCatITPCDiff(anidx).expvars)
            if savefigs && ~pausefigs
                close all
                ifig =figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            numcols = 8;
            numrows = 4; %ceil(length(ntrodes) / numcols);
            for nti = 1:length(ntrodes)
                nt = ntrodes(nti);
                area = ntinfo{1}{1}{nt}.area;
                subarea = ntinfo{1}{1}{nt}.subarea;
                %                 try
                cann = ntinfo{1}{1}{nt}.cannula;
                %                 catch
                %                     fprintf('%d\n',nt)
                %                 end
                ntxy = ntinfo{1}{1}{nt}.ntxy;
                if cann == 'ca1'
                    ntxy(1) = ntxy(1)+4; % offset ca1 to the right
                end
                ntp = (ntxy(2)-1)*8+ntxy(1);
                sf = subaxis(numrows,numcols,ntp, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                    'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                    Pp.MgTp, 'MarginBottom', Pp.MgBm);
                if isnumeric(subarea)
                    subarea = num2str(subarea);
                end
                ntidx = find(matidx == nt);
                idata2plot = squeeze(expvarCatITPCDiff(anidx).ITPCDiff{iv}.ITPC_diff(ntidx,:,:))';
                idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                    'dsamp', expvarCatITPCDiff(anidx).wp.dsamp);
                time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                contourf(sf, time, Fp.wp.frex, idata2plot, Pp.contourRes, ...
                    'linecolor','none');
                set(gca,'ydir','normal','yscale','log');
                colormap(Pp.usecolormap)
                caxis(sf, 'auto')
                hold on
                % thresholded single pix zmask
                if ~isempty(fieldnames(expvarCatITPCDiff(anidx).ITPCDiff{iv}.permt))
                    zmask2plot = squeeze(expvarCatITPCDiff(anidx).ITPCDiff{iv}.permt.threshmean(ntidx,:,:))';
                    zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', Fp.wp.dsamp);
                    try
                        [~,h] = contour(sf, time, Fp.wp.frex, logical(zmask2plot), 1);
                        h.LineColor = 'black';
                    catch
                        fprintf('invalid zmask\n')
                    end
                end
                hold on;
                ytickskip = 2:4:Fp.wp.numfrex;
                set(gca,'ytick', round(Fp.wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('itpcdiff %s-%s %s %s', expvarCatITPCDiff(anidx).expvars{iv}{1},...
                expvarCatITPCDiff(anidx).expvars{iv}{2}, animal, Fp.epochEnvironment);
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 16);
            if pausefigs; pause; end
            if savefigs; save_figure(sprintf('%s/expvarCatITPCDiff/',pconf.andef{4}), sprtit);
            end; end; end; end
%% Plot Stats
% varCont corr plots... testing linear fit
if plot_ContFit; Pp = load_plotting_params({'defaults', 'fitLM'});
    for ian = 1:length(tfbvarCont) % for each animal
        animal = tfbvarCont(ian).animal;
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        conds = {'swrvarContPwrCorr','expvarContPwrCorr'};
        for icorr = 1:length(conds)
            PV = eval(conds{icorr});
            pvanidx = find(strcmp({PV.animal}, animal));
            for iv = 1:length(PV(pvanidx).expvars)
                %             Xvar = PV(pvanidx).dm(:,iv);
                if savefigs && ~pausefigs
                    close all
                    ifig =figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                numcols = 8;
                numrows = 4; %ceil(length(ntrodes) / numcols);
                for itfb = 1:size(PV(pvanidx).tfb)
                    tfb = PV(pvanidx).tfb{itfb};
                    for nti = 1:length(ntrodes)
                        nt = ntrodes(nti);
                        area = ntinfo{1}{1}{nt}.area;
                        subarea = ntinfo{1}{1}{nt}.subarea;
                        cann = ntinfo{1}{1}{nt}.cannula;
                        ntxy = ntinfo{1}{1}{nt}.ntxy;
                        if cann == 'ca1'
                            ntxy(1) = ntxy(1)+4; % offset ca1 to the right
                        end
                        ntp = (ntxy(2)-1)*8+ntxy(1);
                        sf = subaxis(numrows,numcols,ntp, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                            'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                            Pp.MgTp, 'MarginBottom', Pp.MgBm);
                        %                 swrpwr = squeeze(tfbvarCont(ian).dm(:,itfb,nti));
                        %                 [Rc,~] = corr(Xvar, swrpwr);
                        %                 a = fitlm(Xvar, swrpwr);
                        a = PV(pvanidx).fitlm{nti, itfb, iv};
                        b = a.plot; b(1).Marker = 'o'; b(1).MarkerSize = 1; b(1).Color = [.8 .8 .8];
                        legend off
                        axis tight
                        P = a.coefTest; R = a.Rsquared.Ordinary;
                        % conf bounds fill
                        %                 x_axis = b(3).XData';
                        %                 x_plot =[x_axis, fliplr(x_axis)];
                        %                 y_plot = [b(3).ydata', flipud(b(4).y_data)'];
                        %                 y_plot=[CI(:,1)', flipud(CI(:,2))'];
                        %                 hold on
                        %                 plot(x_axis, mu_diff, 'black', 'linewidth', 1)
                        %                 fill(x_plot, y_plot, 1,'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.4);
                        if P < .01 && a.Coefficients.Estimate(end) < 0
                            b(2).Color = [0 0 1];
                            b(3).Color = [0 0 1];
                            b(4).Color = [0 0 1];
                        elseif P < .01 && a.Coefficients.Estimate(end) > 0
                            b(2).Color = [1 0 0 ];
                            b(3).Color = [1 0 0];
                            b(4).Color = [1 0 0];
                        else
                            b(2).Color = [0 0 0];
                            b(3).Color = [0 0 0];
                            b(4).Color = [0 0 0];
                        end
                        title(sprintf('nt%d %s%s r:%.02f p:%.03f',nt,area,subarea,R,P), ...
                            'FontSize', 8);
                        xlabel('')
                        ylabel('')
                        if all(ntxy == [1 1])
                            xlabel(PV(pvanidx).expvars{iv},'FontSize',Pp.FontS, ...
                                'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                            ylabel(tfb,'FontSize',Pp.FontS, ...
                                'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                        end
                    end
                    %% super
                    sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                    sprtit = sprintf('fitLM %s %s VS %s dB', animal, PV(pvanidx).expvars{iv}, tfb);
                    iStitle = text(.5, .98, {sprtit},'Parent',sprtitleax,'Units','normalized');
                    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                        'horizontalAlignment', 'center','FontSize', Pp.stitFsize);
                    %% ---- pause, save figs ----
                    if pausefigs; pause; end;
                    if savefigs; save_figure([pconf.andef{4},'/fitLM/'],sprtit); end
                end
            end
        end
    end
end
% expvarCatDiff Box whisker plots.. testing two distributions
if plot_CatDiffBars; Pp = load_plotting_params({'defaults', 'fitLM'});
    for ian = 1:length(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        invalidtets = evaluatefilter(ntinfo, 'isequal($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        conds = {'expvarCatMeanPwrDiff'};
        for icorr = 1:length(conds)
            PV = eval(conds{icorr});
            pvanidx = find(strcmp({PV.animal}, animal));
            for iv = 1:length(PV(pvanidx).expvars)
                Xvar = logical(PV(pvanidx).dm.dm(:,iv));
                for itfb = 1:length(PV(pvanidx).tfb)
                    if savefigs && ~pausefigs
                        close all
                        ifig =figure('Visible','off','units','normalized','position', ...
                            Pp.position);
                    else
                        ifig = figure('units','normalized','position',Pp.position);
                    end
                    set(gcf,'color','white')
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
                        sf = subaxis(numrows,numcols,ntp, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                            'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                            Pp.MgTp, 'MarginBottom', Pp.MgBm);
                        % i could either have both theta and gamma plots, or i
                        % could compute a xfc measure
                        %                 prethetapwr = squeeze(tfbvarCont(ian).dm(:,3,nti));
                        %                 postthetapwr = squeeze(tfbvarCont(ian).dm(:,2,nti));
                        %                 diffPostPreTheta = postthetapwr - prethetapwr;
                        %                 [R,P] = corr(Xvar, diffPostPreTheta);
                        %                 a = fitlm(Xvar, diffPostPreTheta);
                        %                 b = a.plot;
                        % Xvar vector specifying the rip groupings for this diff cond
                        Avals = PV(pvanidx).Avals{nti, itfb, iv}';
                        Bvals = PV(pvanidx).Bvals{nti, itfb, iv}';
                        mAv = mean(Avals);
                        mBv = mean(Bvals);
                        Astd = std(Avals)/sqrt(length(Avals)); %/numel(Avals);
                        Bstd = std(Bvals)/sqrt(length(Avals)); %/numel(Bvals);
                        br = bar([1 2], [mAv, mBv], 'FaceColor', [0 0 0]);
                        hold on
                        er = errorbar([1 2], [mAv, mBv], [Astd, Bstd]);
                        er.Color = [.3 .3 .3];
                        er.LineStyle = 'none';
                        %                 bp = boxplot([Avals; Bvals], [ones(length(Avals),1); zeros(length(Bvals),1)]);
                        P = PV(pvanidx).P(nti, itfb, iv);
                        %                 ks = PV(pvanidx).kstat(nti, itfb, iv);
                        if P < .01 && mAv > mBv
                            br.FaceColor = [1 0 0];
                            br.EdgeColor = [1 0 0];
                        elseif P < .01 && mAv < mBv
                            
                            br.FaceColor = [0 0 1];
                            br.EdgeColor = [0 0 1];
                        else
                            br.FaceColor = [.5 .5 .5];
                            br.EdgeColor = [.5 .5 .5];
                        end
                        title(sprintf('nt%d %s%s p:%.03f',nt,area,subarea,P), ...
                            'FontSize', 10);
                    end
                    %% super
                    sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                    sprtit = sprintf('tfbVarDiffKS %s %sVS%s %s', animal, PV(pvanidx).expvars{iv}{1}, ...
                        PV(pvanidx).expvars{iv}{2}, PV(pvanidx).tfb{itfb});
                    iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                        'horizontalAlignment', 'center','FontSize', Pp.stitFsize);
                    %% ---- pause, save figs ----
                    if pausefigs
                        pause
                    end
                    if savefigs
                        save_figure(sprintf('%s/%s/',pconf.andef{4}, conds{icorr}), sprtit)
                        close all
                    end
                    close all;
                end
            end
        end
    end
end
% choose Conditions, tfboxes, ntrode areas.
% combine ntrodes by area within animal
% plot TFzmap and either fitLM or diffbars
if plot_CombCorr; Pp = load_plotting_params({'defaults', 'combinedAreasTFstats'});
    for ian = 1:numel(tfbvarContPwrCorr) % for each animal
        animal = tfbvarContPwrCorr(ian).animal;
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        conds = {'expvarContPwrCorr', 'swrvarContPwrCorr'}; %,
        for icorr = 1:length(conds)
            PV = eval(conds{icorr});
            pvanidx = find(strcmp({PV.animal}, animal));
            
            for iv = 1:length(PV(pvanidx).expvars)
                if savefigs && ~pausefigs; close all
                    ifig =figure('Visible','off','units','normalized','position',Pp.position);
                else; ifig = figure('units','normalized','position',Pp.position); end
                set(gcf,'color','white')
                numcols = 3; numrows = 2;  %ceil(length(ntrodes) / numcols);
                areas = {{'mec', 'sup'}, {'mec', 'deep'}, {'ca1', 'd'}};
                c = 1;
                for ia = 1:length(areas)
                    ntsInArea = evaluatefilter(ntinfo, sprintf('isequal($area, ''%s'')', areas{ia}{1}));
                    ntsA = unique(ntsInArea(:,3));
                    meanzmap = squeeze(mean(squeeze(PV(pvanidx).zmap(ntsA,:,:,iv)),1));
                    
                    sftop = subaxis(numrows,numcols,ia, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                        Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    %                 zmap = squeeze(PV(pvanidx).zmap(nti,:,:,iv));
                    %                 zmapthresh = squeeze(PV(pvanidx).clusterZmapThresh(nti, :,:,iv));
                    contourf(PV(pvanidx).time,PV(pvanidx).frequency,meanzmap,40,'linecolor','none')
                    if c
                        xlabel('time')
                        ylabel('freq')
                    else
                        xlabel('');
                        ylabel('');
                    end
                    hold on
                    %                 [~,h] = contour(PV(pvanidx).time,PV(pvanidx).frequency,logical(zmapthresh),1);
                    %                 h.LineColor = 'black';
                    set(gca,'ydir','normal','yscale','log');
                    caxis(sftop, 'auto')
                    colormap(Pp.usecolormap)
                    ytickskip = 2:4:Fp.wp.numfrex;
                    set(gca,'ytick', round(Fp.wp.frex(ytickskip)), 'FontSize', Pp.tickFsize)
                    title(sprintf('%s %s',areas{ia}{1},areas{ia}{2}), 'FontSize',Pp.sfTitFsize,...
                        'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                    yl = ylim;
                    xl = xlim;
                    line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                    userips = ones(length(tfbvarCont(tfbanim).dm(:,1,1)),1);
                    if ~isempty(noiseEvents)
                        % exclude invalid rips
                        noiseanidx = find(strcmp({noiseEvents.animal}, animal));
                        if ~isempty(noiseEvents(noiseanidx).events)
                            invalidrips = ismember([swrvarContPwrCorr(ian).dm.dayeps swrvarContPwrCorr(ian).dm.ripStartTime], ...
                                noiseEvents(noiseanidx).events, 'rows');
                            userips(invalidrips) = 0;
                        end
                    end
                    userips = find(userips);
                    if length(size(PV(pvanidx).dm.dm)) == 3
                        Xvar = mean(PV(pvanidx).dm.dm(userips,iv,ntsInArea),3)';
                    else
                        Xvar = PV(pvanidx).dm.dm(userips,iv)';
                    end
                    tfbani = find(strcmp({tfbvarCont.animal}, animal));
                    for itfb = 1; %:length(tfbvarCont(tfbani).expvars)
                        sftop = subaxis(numrows,numcols,ia, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                            'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                            Pp.MgTp, 'MarginBottom', Pp.MgBm);
                        realfreq = tfbvarCont(tfbani).realfreq{itfb};
                        tfboxYidx = [min(realfreq); max(realfreq)];
                        tfboxXidx = tfbvarCont(tfbani).time{itfb};
                        xw = tfboxXidx(2)-tfboxXidx(1);
                        yh = tfboxYidx(2)-tfboxYidx(1);
                        cx = tfboxXidx(1);
                        cy = tfboxYidx(1);
                        rl = rectangle('Position',[cx cy xw yh]);
                        %                 plot this areas combined fitLM in the second row
                        sfbot = subaxis(numrows,numcols,ia+numcols, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                            'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                            Pp.MgTp, 'MarginBottom', Pp.MgBm);
                        %                     how to get a mean of the fitlm object? or should i
                        %                     plot all nts? no i need to run the LM using combined
                        %                     tfmaps.. so rerun the fit now
                        
                        tfbanim = strcmp({tfbvarCont.animal}, animal);
                        itfbpwr = squeeze(mean(tfbvarCont(tfbanim).dm(userips, itfb, ntsInArea),3))';
                        a = fitlm(Xvar, itfbpwr);
                        %                     a = PV(pvanidx).fitlm{ntsA, itfb, iv};
                        b = a.plot; b(1).Marker = 'o'; b(1).MarkerSize = 1; b(1).Color = [.8 .8 .8];
                        legend off
                        %                     plot(Xvar, )
                        %                     scatter(Xvar', itfbpwr', 'filled');
                        axis tight
                        P = a.coefTest; R = a.Rsquared.Ordinary;
                        if P < .01 && a.Coefficients.Estimate(end) < 0
                            b(2).Color = [0 0 1];
                            b(3).Color = [0 0 1];
                            b(4).Color = [0 0 1];
                        elseif P < .01 && a.Coefficients.Estimate(end) > 0
                            b(2).Color = [1 0 0 ];
                            b(3).Color = [1 0 0];
                            b(4).Color = [1 0 0];
                        else
                            b(2).Color = [0 0 0];
                            b(3).Color = [0 0 0];
                            b(4).Color = [0 0 0];
                        end
                        title(sprintf('r:%.02f p:%.03f',R,P), 'FontSize', 8);
                        xlabel('')
                        ylabel('')
                        
                        if c
                            xlabel(PV(pvanidx).expvars{iv},'FontSize',Pp.FontS, ...
                                'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                            ylabel(tfbvarCont(tfbani).expvars{itfb},'FontSize',Pp.FontS, ...
                                'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                        end
                        c = 0;
                    end
                end
                %% super ax
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('combCorr %s %s', animal, PV(ian).expvars{iv});
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', Pp.stitFsize);
                %% pause, save
                if pausefigs; pause; end
                if savefigs; save_figure(sprintf('%s/%s/',pconf.andef{4}, conds{icorr}), sprtit); end
            end
        end
    end
end

if plot_CombDiff; Pp = load_plotting_params({'defaults', 'combinedAreasTFstats'});
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        conds = {'expvarCatMeanPwrDiff'}; %,
        for icorr = 1:length(conds)
            PV = eval(conds{icorr});
            pvanidx = find(strcmp({PV.animal}, animal));
            for iv = 1:length(PV(pvanidx).expvars)
                if ~strcmp(PV(pvanidx).expvars{iv}{1}, 'distalWell')
                    continue
                end
                if savefigs && ~pausefigs; close all
                    ifig =figure('Visible','off','units','normalized','position',Pp.position);
                else; ifig = figure('units','normalized','position',Pp.position); end
                set(gcf,'color','white')
                numcols = 3; numrows = 3;  %ceil(length(ntrodes) / numcols);
                areas = {{'mec', 'supf'}, {'mec', 'deep'}, {'ca1', 'd'}};
                %                     c = 1;
                for ia = 1:length(areas)
                    % get nt's in this area
                    ntsInArea = evaluatefilter(ntinfo, sprintf(...
                        '(isequal($area, ''%s'')) && (isequal($subarea, ''%s''))', ...
                        areas{ia}{1}, areas{ia}{2}));
                    ntsA = unique(ntsInArea(:,3));
                    [~, ntsA] = ismember(ntsA, ntrodes); %bc D13 doesn't have an nt14
                    ntsA = ntsA(ntsA>0);
                    sftop = subaxis(numrows,numcols,ia, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                        Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    
                    expvarAnix = find(strcmp({expvarCatMeanPwr.animal}, animal));
                    allmean = squeeze(nanmean(expvarCatMeanPwr(...
                        expvarAnix).meandbpower{1}.pwr_mean_db(ntsA,:,:),1))';
                    allmean = trim2win(allmean, Fp.srate, Pp.pwin, ...
                            'dsamp', expvarCatMeanPwr(expvarAnix).wp.dsamp);
                    time = linspace(-Pp.pwin(1),Pp.pwin(2),size(allmean,2));
                    contourf(time,PV(pvanidx).wp.frex,allmean,Pp.contourRes,'linecolor','none')
                    set(gca,'ydir','normal','yscale','log');
                    caxis(sftop,'auto')
                    xlabel('time')
                    ylabel('freq')
                    hold on
                    ytickskip = 2:4:Fp.wp.numfrex;
                    set(gca,'ytick', round(Fp.wp.frex(ytickskip)), 'FontSize', Pp.tickFsize)
                    title(sprintf('%s %s',areas{ia}{1},areas{ia}{2}), 'FontSize',Pp.sfTitFsize,...
                        'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                    yl = ylim; xl = xlim;
                    line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                    tfbani = find(strcmp({tfbvarCont.animal}, animal));
                    cm = hsv(3);
                    for itfb = 1:length(PV(pvanidx).tfb)
                        %                         tfb = PV(pvanidx).tfb{itfb};
                        realfreq = tfbvarCont(tfbani).realfreq{itfb};
                        tfboxYidx = [min(realfreq); max(realfreq)];
                        tfboxXidx = tfbvarCont(tfbani).time{itfb};
                        xw = tfboxXidx(2)-tfboxXidx(1);
                        yh = tfboxYidx(2)-tfboxYidx(1);
                        cx = tfboxXidx(1);
                        cy = tfboxYidx(1);
                        rl = rectangle('Position',[cx cy xw yh], 'EdgeColor', cm(itfb,:));
                    end
                    sfmid = subaxis(numrows,numcols,ia+numcols, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                        Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    meandiffmap=squeeze(nanmean(...
                        PV(pvanidx).meandbpowerDiff{iv}.pwr_meandb_diff(ntsA,:,:),1))';
                    meandiffmap = trim2win(meandiffmap, Fp.srate, Pp.pwin, ...
                            'dsamp', expvarCatMeanPwr(expvarAnix).wp.dsamp);
                    time = linspace(-Pp.pwin(1),Pp.pwin(2),size(meandiffmap,2));
                    contourf(time,PV(pvanidx).wp.frex,meandiffmap,40,'linecolor','none')
                    %                         if c
                    xlabel('time')
                    ylabel('freq')
                    %                         else
                    %                             xlabel('');
                    %                             ylabel('');
                    %                         end
                    hold on
                    %                 [~,h] = contour(PV(pvanidx).time,PV(pvanidx).frequency,logical(zmapthresh),1);
                    %                 h.LineColor = 'black';
                    set(gca,'ydir','normal','yscale','log');
                    caxis(sfmid, 'auto')
                    colormap(Pp.usecolormap)
                    ytickskip = 2:4:Fp.wp.numfrex;
                    set(gca,'ytick', round(Fp.wp.frex(ytickskip)), 'FontSize', Pp.tickFsize)
                    title('diff', 'FontSize',Pp.sfTitFsize,'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                    yl = ylim;
                    xl = xlim;
                    line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                    %                         tfb = PV(pvanidx).tfb{itfb};
                    %                         realfreq = tfbvarCont(tfbani).realfreq{itfb};
                    %                         tfboxYidx = [min(realfreq); max(realfreq)];
                    %                         tfboxXidx = tfbvarCont(tfbani).time{itfb};
                    %                         xw = tfboxXidx(2)-tfboxXidx(1);
                    %                         yh = tfboxYidx(2)-tfboxYidx(1);
                    %                         cx = tfboxXidx(1);
                    %                         cy = tfboxYidx(1);
                    cm = hsv(3);
                    for itfb = 1:length(PV(pvanidx).tfb)
                        %                         tfb = PV(pvanidx).tfb{itfb};
                        realfreq = tfbvarCont(tfbani).realfreq{itfb};
                        tfboxYidx = [min(realfreq); max(realfreq)];
                        tfboxXidx = tfbvarCont(tfbani).time{itfb};
                        xw = tfboxXidx(2)-tfboxXidx(1);
                        yh = tfboxYidx(2)-tfboxYidx(1);
                        cx = tfboxXidx(1);
                        cy = tfboxYidx(1);
                        rl = rectangle('Position',[cx cy xw yh], 'EdgeColor', cm(itfb,:));
                    end
                    %                     userips = ones(length(tfbvarCont(tfbanim).dm(:,1,1)),1);
                    %                     if ~isempty(noiseEvents)
                    %                         % exclude invalid rips
                    %                         noiseanidx = find(strcmp({noiseEvents.animal}, animal));
                    %                         if ~isempty(noiseEvents(noiseanidx).events)
                    %                             invalidrips = ismember([swrvarContPwrCorr(ian).dm.dayeps swrvarContPwrCorr(ian).dm.ripStartTime], ...
                    %                                 noiseEvents(noiseanidx).events, 'rows');
                    %                             userips(invalidrips) = 0;
                    %                         end
                    %                     end
                    %                     userips = find(userips);
                    %                     if length(size(PV(pvanidx).dm.dm)) == 3
                    %                         Xvar = mean(PV(pvanidx).dm.dm(userips,iv,ntsInArea),3)';
                    %                     else
                    %                         Xvar = PV(pvanidx).dm.dm(userips,iv)';
                    %                     end
                    %                     tfbani = find(strcmp({tfbvarCont.animal}, animal));
                    %                     for itfb = 1; %:length(tfbvarCont(tfbani).expvars)
                    %                         sftop = subaxis(numrows,numcols,ia, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                    %                             'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                    %                             Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    %                         realfreq = tfbvarCont(tfbani).realfreq{itfb};
                    %                         tfboxYidx = [min(realfreq); max(realfreq)];
                    %                         tfboxXidx = tfbvarCont(tfbani).time{itfb};
                    %                         xw = tfboxXidx(2)-tfboxXidx(1);
                    %                         yh = tfboxYidx(2)-tfboxYidx(1);
                    %                         cx = tfboxXidx(1);
                    %                         cy = tfboxYidx(1);
                    %                         rl = rectangle('Position',[cx cy xw yh]);
                    sfbot = subaxis(numrows,numcols,ia+(numcols*2), 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                        Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    
                    %                     how to get a mean of the fitlm object? or should i
                    %                     plot all nts? no i need to run the LM using combined
                    %                     tfmaps.. so rerun the fit now
                    
                    %                         tfbanim = strcmp({tfbvarCont.animal}, animal);
                    %                         itfbpwr = squeeze(mean(tfbvarCont(tfbanim).dm(userips, itfb, ntsInArea),3))';
                    %                         tfbani = find(strcmp({tfbvarCont.animal}, animal));
                    mAv = [];
                    mBv = [];
                    Asem = [];
                    Bsem = [];
                    for itfb = 1:length(PV(pvanidx).tfb)
                        Avals = cell2mat(PV(pvanidx).Avals(ntsA, itfb, iv)');
                        Bvals = cell2mat(PV(pvanidx).Bvals(ntsA, itfb, iv)');
                        mAv(itfb) = nanmean(Avals);
                        mBv(itfb) = nanmean(Bvals);
                        Asem(itfb) = std(Avals)/sqrt(length(Avals)); %/numel(Avals);
                        Bsem(itfb) = std(Bvals)/sqrt(length(Avals)); %/numel(Bvals);
                    end
                    br = bar([mAv; mBv], 'Barwidth', 1); hold on
                    colormap(sfbot, hsv(3))
                    xticks([1 2]);
                    xticklabels({'distal', 'proximal'});
                    hold on;
                    hAx=gca;
                    X=cell2mat(get(br,'XData')).'+[br.XOffset];  
%                     hEB = errorbar(X,y,std_dev,'.') 
                    er = errorbar(X, [mAv; mBv], [Asem; Bsem], 'k.');
                    %                         er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; er(3).LineStyle = 'none';
                    %                 bp = boxplot([Avals; Bvals], [ones(length(Avals),1); zeros(length(Bvals),1)]);
                    %                         P = nanmean(PV(pvanidx).P(ntsA, itfb, iv));
                    %                 ks = PV(pvanidx).kstat(nti, itfb, iv);
                    %                         if P < .01 && mAv > mBv
                    %                             br.FaceColor = [1 0 0];
                    %                             br.EdgeColor = [1 0 0];
                    %                         elseif P < .01 && mAv < mBv
                    %
                    %                             br.FaceColor = [0 0 1];
                    %                             br.EdgeColor = [0 0 1];
                    %                         else
                    %                             br.FaceColor = [.5 .5 .5];
                    %                             br.EdgeColor = [.5 .5 .5];
                    %                         end
%                     title(sprintf('%s%s',areas{ia}{1}, areas{ia}{2}), ...
%                         'FontSize', 10);
                end
                %% super ax
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('tfbVarDiffKS %s %sVS%s', animal, PV(pvanidx).expvars{iv}{1}, ...
                    PV(pvanidx).expvars{iv}{2});
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', Pp.stitFsize);
                %% pause, save
                if pausefigs; pause; end
                if savefigs; save_figure(sprintf('%s/byArea/%s/',pconf.andef{4}, conds{icorr}), sprtit); end
            end; end; end; end