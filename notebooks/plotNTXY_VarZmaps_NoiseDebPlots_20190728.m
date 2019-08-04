

%{
Part 1
Hypothesis: TFcontext, expVar --> swrpropagation:(FastGamma, Ripple, MUspikesSWR)

Monday:
    - fix whatever is happening to mislabel the ntrode numbers, as apparent with D10

    - then noise debugging use the all rips mean power plots and the strips
plots.. fix references and exclude noise, then assign ntrodes and go from
there

    - plot diff of expvarcats
    - plot itpc for expvarcat,  diffs

    - combine frequiency ranges for the strips plots.. i.e. use the freq lookups
from tdbvarCont maker to get mean of frequencies within each tfbox (gammas, theta, etc)  +
(250-300Hz noise)
        - ?? add the single ripple loader and plotter to this.. need to print out all
        the individual rips again to comb through for noise debugging

%}
% LFP
make_swrLFP = 0;
save_swrLFP = make_swrLFP;
load_swrLFP = 0;
stack_swrLFP = make_swrLFP;
load_swrLFPstack = 0;

% Power Phase
make_powerPhase = 0;
load_rawpwr = 0;
load_phase = 0;

% Des Mats
loadDesMats = 1;
makeDesMats = 0;
make_expvarCat = makeDesMats;
load_expvarCat = loadDesMats;
make_expvarCont = makeDesMats;
load_expvarCont = loadDesMats;
make_swrvarCont = makeDesMats;
load_swrvarCont = loadDesMats;
make_tfbvarCont = makeDesMats;
load_tfbvarCont = loadDesMats;

% Agg ops and permtest ********
make_expvarCatMeanPwr = 0;
make_expvarCatMeanPwrDiff = 0; % need to finish writing this, with perm testing
make_varContPwrCorr = 0; % three calls, expvarCont, swrvarCont tfbvarCont

make_itpc = 0;
make_expvarCatITPCDiff = 0;

runStats = 0;

% load results
load_expvarCatMeanPwr = 0;
load_expvarCatDiff = 0;
load_varContPwrCorr = 0;
load_itpc = 0;
load_itpcDiff = 0;

% plot
plot_meanpwr = 0; % agg
plot_meanpwrDiff = 0; %agg 
plot_corrpwr = 0; %agg
plot_pwr_timeXrip = 0; % strips

plot_itpc = 0; %agg
plot_ITPCDiff = 0; %agg
plot_phase_timeXrip = 0; %strips

plot_ContFit = 0;
plot_CatDiffBars = 1;

pausefigs = 1;
savefigs = 0;

Fp.animals = {'D10', 'D12', 'JZ1', 'JZ4'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp.add_params = {'wtrackdays', 'excludeNoise','excludePriorFirstWell', '<4cm/s', ...
    'wavelets4-300Hz',  'excludeAfterLastWell'};
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
Fp.uselfptype = 'eeg';
% Fp.useExpvars = {'onlywdays','rewarded', 'unrewarded', 'inbound' , 'outbound', ...
%     'proximalWell', 'distalWell'};
% Fp.useDiffTrialTypes = {'diffRewardedUnrewarded', 'diffOutboundInbound'};

pconf = paramconfig;
me = animaldef('Demetris');
wp = getWaveParams(Fp.waveSet);
%% riptrigLFP m.s.l.
    % make.save.load riptrigLFP
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
if stack_swrLFP; lfpstack = stack_riptriglfp(F); clear F
    save_data(lfpstack, Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack'); end
if load_swrLFPstack
    lfpstack = load_data(Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack', Fp.animals);
end

%% Power, Phase m.s.l
    % make.save.load rawpower, phase [nt times rip freq]
if make_powerPhase; computeAnalyticSignal(lfpstack,'waveSet',Fp.waveSet,'saveOutput',1, ...
        'uselfptype', Fp.uselfptype); end
if load_rawpwr; rawpwr = load_data(sprintf('%s/analyticSignal/', me{2}), ...
        sprintf('AS_waveSet-%s_%s_power', wp.waveSet, Fp.uselfptype), Fp.animals);end
if load_phase; phase = load_data(sprintf('%s/analyticSignal/', me{2}), ...
        sprintf('AS_waveSet-%s_%s_phase', wp.waveSet, Fp.uselfptype), Fp.animals);end

%% Design Matrices [rip var (nt) (ntB)]
    % :expvarCat @mean @ITPC /rip $time [rip var] ...
        % all, outbound, inbound, rewarded,  unrewarded, proximalWell, distalWell
if make_expvarCat; outdir = 'expvarCat'; 
    expvarCat = makeExpvarCatDesignMat(lfpstack, 'outdir', outdir); end
if load_expvarCat; outdir = 'expvarCat'; outpath = [pconf.andef{2},outdir,'/'];
    expvarCat = load_data(outpath, [outdir,'_',Fp.epochEnvironment], Fp.animals);end

    % :expvarCont @corr /rip $swr [rip var] ...
        % hd, speed, xpos, ypos, learn, perform, timeDay, timeEp, ripnum, day, epoch
if make_expvarCont; outdir = 'expvarCont'; 
    expvarCont = makeExpvarContDesignMat(lfpstack, Fp, 'outdir', outdir); end
if load_expvarCont; outdir = 'expvarCont'; outpath = [pconf.andef{2},outdir,'/'];
    expvarCont = load_data(outpath, [outdir,'_',Fp.epochEnvironment], Fp.animals); end
    
    % :swrvarCont @corr /rip $swr [rip var] std, total_energy, duration
if make_swrvarCont; outdir = 'swrvarCont';
    swrvarCont = makeSWRDesignMatrix(lfpstack, Fp, 'outdir', outdir); end
if load_swrvarCont; outdir = 'swrvarCont'; outpath = [pconf.andef{2},outdir,'/'];
    swrvarCont = load_data(outpath, [outdir,'_',Fp.epochEnvironment], Fp.animals); end

    % :tfbvarCont @corr /rip $swr [rip var nt] ...
        % {pre prepre swr post postpost, rip gammas beta theta}
if make_tfbvarCont; outdir = 'tfbvarCont';
    tfbvarCont = makeTFBoxDesignMat(rawpwr, Fp, 'outdir', outdir); end
if load_tfbvarCont; outdir = 'tfbvarCont'; outpath = [pconf.andef{2},outdir,'/'];
    tfbvarCont = load_data(outpath, [outdir,'_',Fp.epochEnvironment], Fp.animals); end

%% Aggregate Operations on rawpower, phase data by Design Matrices
    % :expvarCat @mean /ntTF $time m.s.l
    % don't run perm test for means bc time shuffling doesn't work for tonic signal
if make_expvarCatMeanPwr; getPower(expvarCat,rawpwr,Fp, 'run_perm', 0); end

% :expvarCatDiff @dmean /ntTF $var m.s.l outbound-inbound, rewarded-unrewarded, ...
        % proximalWell-distalWell, wtrackEp2-wtrackEp4
if make_expvarCatMeanPwrDiff
    getPowerDiff(expvarCat,rawpwr,Fp, 'run_perm', 0); 
end

    % @corr (expvarCont tfbvarCont swrvarCont) /ntTF $swr m.s.l
if make_varContPwrCorr
    expvarContPwrCorr = runDesignDataRegression(expvarCont,rawpwr,Fp,'outdir','expvarCont');
    swrvarContPwrCorr = runDesignDataRegression(swrvarCont,rawpwr,Fp,'outdir','swrvarCont');
    tfbvarContPwrCorr = runDesignDataRegression(tfbvarCont,rawpwr,Fp,'outdir','tfbvarCont');
end

    % :expvarCat @itpc /ntTF $time m.s.l
if make_itpc; getITPC(expvarCat, phase, Fp,'run_perm', 0); end
%     :expvarCatDiff @ditpc /ntTF $var m.s.l
if make_expvarCatITPCDiff; makeExpvarCatITPCDiff(expvarCat, phase, Fp, 'run_perm', 1); end


%% loading and plotting Aggregate Results Exploratory
if load_expvarCatMeanPwr; outdir = 'expvarCatMeanPwr'; outpath = [pconf.andef{2},outdir,'/'];
    meanpwr = load_data(outpath, [outdir,'_', Fp.epochEnvironment] ,Fp.animals); end
if load_expvarCatDiff; outdir = 'expvarCatMeanPwrDiff'; outpath = [pconf.andef{2},outdir,'/'];
    expvarCatDiff = load_data(outpath, [outdir,'_',Fp.epochEnvironment], Fp.animals);end
if load_varContPwrCorr
    outdir = 'expvarCont'; outfile='expvarContCorr'; outpath=[pconf.andef{2},outdir,'/'];
    expvarContPwrCorr=load_data(outpath,[outfile,'_', Fp.epochEnvironment],Fp.animals);
    outdir = 'swrvarCont'; outfile='swrvarContCorr'; outpath=[pconf.andef{2},outdir,'/'];
    swrvarContPwrCorr=load_data(outpath,[outfile,'_', Fp.epochEnvironment],Fp.animals);
    outdir = 'tfbvarCont'; outfile='tfbvarContCorr'; outpath=[pconf.andef{2},outdir,'/'];
    tfbvarContPwrCorr=load_data(outpath,[outfile,'_', Fp.epochEnvironment],Fp.animals);
end

if load_itpc; outdir = 'expvarCatITPC'; outpath = [pconf.andef{2},outdir,'/'];
   itpc = load_data(outpath, [outdir,'_', Fp.epochEnvironment] ,Fp.animals); end
if load_itpcDiff; outdir = 'expvarCatITPCDiff'; outpath = [pconf.andef{2},outdir,'/'];
   itpcdiff = load_data(outpath, [outdir,'_', Fp.epochEnvironment] ,Fp.animals); end

    % meanpowerDiff varCat TFzmap /nt
if plot_meanpwrDiff; Pp = load_plotting_params({'defaults', 'powerTFmap'});
    animals = {expvarCatDiff.animal};
    for ian = 1:length(animals) % for each animal
        animal = animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo, 'area');
        matidx = unique(den.index(:,3));
        anidx = find(strcmp({expvarCatDiff.animal}, animal));
        for co = 1:length(expvarCatDiff(anidx).expvars)
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
                idata2plot = squeeze(expvarCatDiff(anidx).meandbpowerDiff{co}.pwr_meandb_diff(ntidx,:,:))';
                idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                    'dsamp', expvarCatDiff(anidx).wp.dsamp);
                time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                contourf(sf, time, wp.frex, idata2plot, Pp.contourRes, ...
                    'linecolor','none');
                set(gca,'ydir','normal','yscale','log');
                
                colormap(Pp.usecolormap)
                caxis(sf, 'auto')
%                 colorbar
                
                hold on
                % thresholded single pix zmask
                if ~isempty(fieldnames(expvarCatDiff(anidx).meandbpowerDiff{co}.permt))
                    zmask2plot = squeeze(expvarCatDiff(anidx).meandbpowerDiff{co}.permt.threshmean(ntidx,:,:))'; 
                    zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', wp.dsamp);
                    try
                        [~,h] = contour(sf, time, wp.frex, logical(zmask2plot), 1);
                        h.LineColor = 'black';
                    catch
                        fprintf('invalid zmask\n')
                    end
                end
                hold on;
                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)       
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('meanpwrdiff %s-%s %s %s', expvarCatDiff(anidx).expvars{co}{1},...
                expvarCatDiff(anidx).expvars{co}{2}, animal, Fp.epochEnvironment);
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 16);
            
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(sprintf('%s/expvarCatMeanPwrDiff/',pconf.andef{4}), 'expvarCatMeanPwrDiff', sprtit)
                close all
            end
            close all;
            %                 end
        end
    end
end
    % meanpower varCat TFzmap /nt
if plot_meanpwr; Pp = load_plotting_params({'defaults', 'powerTFmap'});
    animals = {meanpwr.animal};
    for ian = 1:length(animals) % for each animal
        animal = animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo, 'area');
        matidx = unique(den.index(:,3));
        anidx = find(strcmp({meanpwr.animal}, animal));
        for co = 1:length(meanpwr(anidx).expvars)
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
                idata2plot = squeeze(meanpwr(anidx).meandbpower{co}.pwr_mean_db(ntidx,:,:))';
                idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                    'dsamp', meanpwr(anidx).wp.dsamp);
                time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                contourf(sf, time, wp.frex, idata2plot, Pp.contourRes, ...
                    'linecolor','none');
                set(gca,'ydir','normal','yscale','log');
                
                colormap(Pp.usecolormap)
                caxis(sf, 'auto')
%                 colorbar
                
                hold on
                % thresholded single pix zmask
                if ~isempty(fieldnames(meanpwr(anidx).meandbpower{co}.permt))
                    zmask2plot = squeeze(meanpwr(anidx).meandbpower{co}.permt.threshmean(ntidx,:,:))'; 
                    zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', wp.dsamp);
                    try
                        [~,h] = contour(sf, time, wp.frex, logical(zmask2plot), 1);
                        h.LineColor = 'black';
                    catch
                        fprintf('invalid zmask\n')
                    end
                end
                hold on;
                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)       
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('meanpwr %s %s %s', meanpwr(anidx).expvars{co}, animal, ...
                Fp.epochEnvironment);
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 16);
            
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(sprintf('%s/expvarCatMeanPwr/',pconf.andef{4}), 'expvarCatMeanPwr', sprtit)
                close all
            end
            close all;
            %                 end
        end
    end
end
    % corrpower varCont TFzmap /nt
if plot_corrpwr; Pp = load_plotting_params({'defaults', 'powerTFmap'});
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        conds = {'tfbvarContPwrCorr','expvarContPwrCorr', 'swrvarContPwrCorr'};
        for icorr = 1:length(conds)
            PV = eval(conds{icorr});
        for iv = 1:length(PV(ian).expvars)
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
                
                zmap = squeeze(PV(ian).zmap(nti,:,:,iv));
                zmapthresh = squeeze(PV(ian).clusterZmapThresh(nti, :,:,iv));
                contourf(PV(ian).time,PV(ian).frequency,zmap,40,'linecolor','none')
                hold on
                [~,h] = contour(PV(ian).time,PV(ian).frequency,logical(zmapthresh),1);
                h.LineColor = 'black';
                set(gca,'ydir','normal','yscale','log');
                caxis(sf, 'auto')
                                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', Pp.tickFsize)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',Pp.sfTitFsize,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)       
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end            
                            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('powerVarCorr %s %s', animal, PV(ian).expvars{iv});
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', Pp.stitFsize);
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(sprintf('%s/%s/',pconf.andef{4}, conds{icorr}), conds{icorr}, sprtit)
                close all
            end
            close all;
        end
        end
    end   
end
    % ITPC varCat TFzmap /nt
if plot_itpc; Pp = load_plotting_params({'defaults', 'powerTFmap'});
        animals = {itpc.animal};
    for ian = 1:length(animals) % for each animal
        animal = animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo, 'area');
        matidx = unique(den.index(:,3));
        anidx = find(strcmp({itpc.animal}, animal));
        for co = 1:length(itpc(anidx).expvars)
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
                idata2plot = squeeze(itpc(anidx).ITPC{co}.ITPC_db(ntidx,:,:))';
                idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                    'dsamp', itpc(anidx).wp.dsamp);
                time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                contourf(sf, time, wp.frex, idata2plot, 40, ...
                    'linecolor','none');
                set(gca,'ydir','normal','yscale','log');
                
                colormap(Pp.usecolormap)
                caxis(sf, 'auto')
%                 colorbar
                
                hold on
                % thresholded single pix zmask
                if ~isempty(fieldnames(itpc(anidx).ITPC{co}.permt))
                    zmask2plot = squeeze(itpc(anidx).ITPC{co}.permt.threshmean(ntidx,:,:))'; 
                    zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', wp.dsamp);
                    try
                        [~,h] = contour(sf, time, wp.frex, logical(zmask2plot), 1);
                        h.LineColor = 'black';
                    catch
                        fprintf('invalid zmask\n')
                    end
                end
                hold on;
                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)       
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('itpc %s %s %s', itpc(anidx).expvars{co}, animal, ...
                Fp.epochEnvironment);
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 16);
            
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(sprintf('%s/expvarCatITPC/',pconf.andef{4}), 'expvarCatITPC', sprtit)
                close all
            end
            close all;
            %                 end
        end
    end
end
    % ITPCDiff varCat TFzmap /nt
if plot_ITPCDiff; Pp = load_plotting_params({'defaults', 'powerTFmap'});
    animals = {itpcdiff.animal};
    for ian = 1:length(animals) % for each animal
        animal = animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo, 'area');
        matidx = unique(den.index(:,3));
        anidx = find(strcmp({itpcdiff.animal}, animal));
        for co = 1:length(itpcdiff(anidx).expvars)
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
                idata2plot = squeeze(itpcdiff(anidx).ITPCDiff{co}.ITPC_diff(ntidx,:,:))';
                idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                    'dsamp', itpcdiff(anidx).wp.dsamp);
                time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                contourf(sf, time, wp.frex, idata2plot, Pp.contourRes, ...
                    'linecolor','none');
                set(gca,'ydir','normal','yscale','log');
                
                colormap(Pp.usecolormap)
                caxis(sf, 'auto')
%                 colorbar
                
                hold on
                % thresholded single pix zmask
                if ~isempty(fieldnames(itpcdiff(anidx).ITPCDiff{co}.permt))
                    zmask2plot = squeeze(itpcdiff(anidx).ITPCDiff{co}.permt.threshmean(ntidx,:,:))'; 
                    zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', wp.dsamp);
                    try
                        [~,h] = contour(sf, time, wp.frex, logical(zmask2plot), 1);
                        h.LineColor = 'black';
                    catch
                        fprintf('invalid zmask\n')
                    end
                end
                hold on;
                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)       
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('itpcdiff %s-%s %s %s', itpcdiff(anidx).expvars{co}{1},...
                itpcdiff(anidx).expvars{co}{2}, animal, Fp.epochEnvironment);
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 16);
            
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(sprintf('%s/expvarCatITPCDiff/',pconf.andef{4}), ...
                    'expvarCatITPCDiff', sprtit)
                close all
            end
            close all;
            %                 end
        end
    end
end
    % strips plot power heatRast /nt /freq
if plot_pwr_timeXrip; Pp=load_plotting_params({'defaults','powerheatRast'});
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
                    save_figure([pconf.andef{4},'/powerHeatRast/'],'powerHeatRast',sprtit)
                    close all
                end
            end
        end
    end
end

%% SwrDuration, SwrSTD corr mecRipPwr
% lin fit corr coef plots... testing linear relationship
if plot_ContFit; Pp = load_plotting_params({'defaults', 'powerTFmap'});
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        conds = {'swrvarCont'}; %'expvarCont', 
        for icorr = 1:length(conds)
            PV = eval(conds{icorr});
        for iv = 1:length(PV(ian).expvars)
            Xvar = PV(ian).dm(:,iv);
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
                cann = ntinfo{1}{1}{nt}.cannula;
                ntxy = ntinfo{1}{1}{nt}.ntxy;
                if cann == 'ca1'
                    ntxy(1) = ntxy(1)+4; % offset ca1 to the right
                end
                ntp = (ntxy(2)-1)*8+ntxy(1);
                sf = subaxis(numrows,numcols,ntp, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                    'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                    Pp.MgTp, 'MarginBottom', Pp.MgBm);

                swrpwr = squeeze(tfbvarCont(ian).dm(:,1,nti));
                [R,P] = corr(Xvar, swrpwr);
                a = fitlm(Xvar, swrpwr);
                b = a.plot;
                b(1).Marker = 'o';
                b(1).MarkerSize = 1;
                b(1).Color = [.8 .8 .8];
                if P < .05 && R < 0
                    b(2).Color = [0 0 1];
                    b(3).Color = [0 0 1];
                    b(4).Color = [0 0 1];
                elseif P < .05 && R > 0
                    b(2).Color = [1 0 0 ];
                    b(3).Color = [1 0 0];
                    b(4).Color = [1 0 0];
                else
                    b(2).Color = [0 0 0];
                    b(3).Color = [0 0 0];
                    b(4).Color = [0 0 0];
                end
                legend off
                axis tight
%                 p = a.coefTest;
                title(sprintf('nt%d %s%s r:%.02f p:%.03f',nt,area,subarea,R,P), ...
                    'FontSize', 10);
            end            
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('fitlm %s %s', animal, PV(ian).expvars{iv});
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', Pp.stitFsize);
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(sprintf('%s/%s/',pconf.andef{4}, conds{icorr}), conds{icorr}, sprtit)
                close all
            end
            close all;
        end
        end
    end
end

%% outbThetaPostPre VS inbThetaPostPre also distal vs proximal
% Box whisker plots.. testing two distributions

if plot_CatDiffBars; Pp = load_plotting_params({'defaults', 'powerTFmap'});
    expvarCatDiff = struct;
    for ian = 1:numel(Fp.animals) % for each animal
        expvarCatDiff(ian).expvars = {'rewarded', 'inbound', 'proximal'};
        expvarCatDiff(ian).dm = expvarCat(ian).dm(:,[2 4 6]);
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        conds = {'expvarCatDiff'};
        for icorr = 1:length(conds)
            PV = eval(conds{icorr});
        for iv = 1:length(PV(ian).expvars)
            Xvar = logical(PV(ian).dm(:,iv));
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
                prethetapwr = squeeze(tfbvarCont(ian).dm(:,3,nti));
                postthetapwr = squeeze(tfbvarCont(ian).dm(:,2,nti));
                diffPostPreTheta = postthetapwr - prethetapwr;
%                 [R,P] = corr(Xvar, diffPostPreTheta);
%                 a = fitlm(Xvar, diffPostPreTheta);
%                 b = a.plot;
                % Xvar vector specifying the rip groupings for this diff cond
                bp = bar(diffPostPreTheta, Xvar);
                [h,P,ks2stat] = kstest2(diffPostPreTheta(Xvar), diffPostPreTheta(~Xvar));
                if P < .05
                    pause
                end
                legend off
%                 p = a.coefTest;
                title(sprintf('nt%d %s%s p:%.03f',nt,area,subarea,R,P), ...
                    'FontSize', 10);
            end            
                            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('fitlm %s %s', animal, PV(ian).expvars{iv});
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', Pp.stitFsize);
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(sprintf('%s/%s/',pconf.andef{4}, conds{icorr}), conds{icorr}, sprtit)
                close all
            end
            close all;
        end
        end
    end
end

%{
%%%%%%%%%%%%%% TODO:
-- fix the expvarCont values.. they are wrong.. 
-- the duration values look discrete.. is this bc rip chains?
-- plot the pre and post TFboxes for the diff comparisons

-- need to choose ntrodes
    -- print out all the animals zmapXY
        -- outbrew-inbound thetapost-pre boxplots
        -- ? distal-prox thetapost-pre boxplots
        -- swr std corr mecrip
        -- swr duration corr mecrip
        -- perform corr mecrip

-- measure for theta/gamma xfc??
-- speed is wrong.. i think all of the expvarCont's are wrong ???
-- jz3 still bunk
-- D10 (and maybe others, offset nts by 1)

-- need Cat, rip chain or not to complement rip duration
-- get examples for each animal. get ntrode histo for each animal

-- need to show as much as i can with spikes as well
-- then finale with ca1 replay corr
-- replay correspondance corr?

-- mutual information/ entropy measure?


-- plot the wavelength info figs from that other paper using my wavelet
    params
-- change to time unit instead of ncycles.. ?


%%%%%%%%% RESULT OBSERVATIONS
Cont CORR
    
SWRcorr    
    - swrTFBOX : D10 high 0-10Hz for all animals.. but is this fast gamma?
        if it was then you'd expect high 125-250 for the rest of the window
        too.. but there are definately ntrodes where the 125-250 power is
        punctate (0-200ms) but the theta is high throughout
        *** theta power correlates with ripple band power in MEC
            during ca1 SWRs

    - STD: pretty consistent 20-50 Hz power increase across the entire period
        - significant 6-10 Hz decrease at the time of the swr.. but also
        extending during more of the window
        *** slow gamma power related to ca1 SWR std, but theta power
            inversely related to strength of ripple..

    - duration D10, JZ2, D13 increase in ripple band power. D13, JZ3
        decrease 0-40Hz
        *** longer ripples have more MEC ripple band power

EXPcorr
    - the speed corr is definately weird.. i thought before it was consitent
        with theta and fast gamma being high,, but that's not really that
        strong anymore.. 

    - performance high at 10Hz for D10, D13, JZ1, JZ4. maybe D12.
        *** theta power at swr time related to performance level

Cat Mean
    - outbound looks like it has more 40-85 Hz with 10 Hz stuff following the ripple
        - def true of D10, nt11 which has grid cells
        - D12, nt1, 2, 5, 8
        - not D13, JZ2, 3, 4
        - JZ1 nt1
        *** gamma/theta increase post ripple on outbound.. indicative of
            swr entering working memory?

    - distal vs proximal well mean.. more gamma, theta after rip
        - D10, lots of tets
        - D12, D13, JZ1, JZ2, a few 
        - JZ4 a couple
        *** gamma/theta increase post rip away from well, indicative of swr
            entering working memory?


%%%%%%%%%%% story: SWRS INTERACT WITH ONGOING MEMORY PROCESSING
- MEC response during ca1 swrs changes depending on the state of MEC, and
the state of the animal.

== SWR propagation
: MEC response to ca1 swr is heterogeneous
: SWR Propagation measured with ripple band power in MEC
: SWR STD is not correlated to SWR Propagation 
: SWR duration is correlated to SWR Propagation 

== SWR functional context 
: SWR Propagation correlated with theta outside of sharp wave (+-.5 s)

== SWR Integration
: SWR Integration measured with Post-Pre theta/gamma xfc
: Post-Pre theta/gamma xfc correlates with performance
: Post-Pre theta/gamma xfc greater during outboundRewarded vs inboundRewarded

controls
: Speed is NOT related to Post-Pre theta/gamma xfc

Other measures
?: theta phase reset? theta phase consistency Post vs Pre?



%}

