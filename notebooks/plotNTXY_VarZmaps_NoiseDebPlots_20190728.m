


%{
Part 1
Hypothesis: TFcontext, expVar --> swrpropagation:(FastGamma, Ripple, MUspikesSWR)

Monday:
#1.. get the meanpower and pwrcorr plots working (how does duration/std
vars look?)
then:
=== cXYplots: plot ntrodes at relative cannula position =====
    - fix whatever is happening to mislabel the ntrode numbers, as apparent with D10
    - add ntrode relative XY plot position to tetinfo struct based on each animals cannula
    - plot the ntrodes in the xy subfigure position specified in the tetinfostruct

then noise debugging use the all rips mean power plots and the strips
plots.. fix references and exclude noise, then assign ntrodes and go from
there

combine frequiency ranges for the strips plots.. i.e. use the freq lookups
from tdbvarCont maker to get mean of frequencies within each tfbox (gammas, theta, etc)  +
(250-300Hz noise)

- add the single ripple loader and plotter to this.. need to print out all
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
loadDesMats = 0;
makeDesMats = 0;
make_expvarCat = makeDesMats;
load_expvarCat = loadDesMats;
make_expvarCont = makeDesMats;
load_expvarCont = loadDesMats;
make_swrvarCont = makeDesMats;
load_swrvarCont = loadDesMats;
make_tfbvarCont = makeDesMats;
load_tfbvarCont = loadDesMats;

% Agg ops and permtest
make_expvarCatMeanPwr = 0;
% make_expvarCatMeanPwrDiff = 0; % need to finish writing this, with perm testing
make_varContPwrCorr = 0; % three calls, expvarCont, swrvarCont tfbvarCont

% make_itpc = 0;
% make_expvarCatMeanItpcDiff = 0;

% load results
load_expvarCatMeanPwr = 0;
load_varContPwrCorr = 0;

load_expvarCatDiff = 0;
load_itpc = 0;

% plot
plot_meanpwr_timeXfreq = 0; % agg
plot_corrpwr_timeXfreq = 1; %agg
plot_pwr_timeXrip = 0; % strips

plot_diffmeanpwer_timeXfreq = 0; %agg
plot_phase_timeXrip = 0; %strips
plot_itpc_timeXfreq = 0; %agg
plot_diffitpc_timeXfreq = 0; %agg

pausefigs = 0;
savefigs = 1;

Fp.animals = {'D10'};
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
    % don't run perm test for means bc time shuffling doesn't work for
    % tonic signal
if make_expvarCatMeanPwr; getPower(expvarCat,rawpwr,Fp, 'run_permutation_test', 0); end
    % :expvarCatDiff @dmean /ntTF $var m.s.l outbound-inbound, rewarded-unrewarded, ...
        % proximalWell-distalWell, wtrackEp2-wtrackEp4
% if make_expvarCatMeanPwrDiff; getExpVarDiff(expvarCatDiff); end

    % @corr (expvarCont tfbvarCont swrvarCont) /ntTF $swr m.s.l
if make_varContPwrCorr
    expvarContPwrCorr = runDesignDataRegression(expvarCont,rawpwr,Fp,'outdir','expvarCont');
    swrvarContPwrCorr = runDesignDataRegression(swrvarCont,rawpwr,Fp,'outdir','swrvarCont');
    tfbvarContPwrCorr = runDesignDataRegression(tfbvarCont,rawpwr,Fp,'outdir','tfbvarCont');
end

    % :expvarCat @itpc /ntTF $time m.s.l
% if make_itpc; getITPC(expvarCat, Fp, 'expvars', Fp.useExpvars); end
    % :expvarCatDiff @ditpc /ntTF $var m.s.l
% if make_expvarCatITPCDiff; makeExpvarCatITPCDiff(expvarCatDiff); end

% currently i'm not planning on doing phase corr for cont vars.. but should i?

%% loading and plotting Aggregate Results Exploratory
if load_expvarCatMeanPwr; outdir = 'expvarCatMeanPwr'; outpath = [pconf.andef{2},outdir,'/'];
    meanpwr = load_data(outpath, [outdir,'_', Fp.epochEnvironment] ,Fp.animals); end
% if load_itpc; savedir = sprintf('%s/itpc/', pconf.andef{2});
%     savestr=sprintf('/itpc_waveSet-%s_%s_%s',Fp.waveSet,Fp.uselfptype,Fp.epochEnvironment);
%     itpc = load_data(savedir, savestr, Fp.animals); end
% if load_expvarCatDiff; outdir = 'expvardiff'; outpath = [pconf.andef{2},outdir,'/'];
%     expvarCatDiff = load_data(outpath, [outdir,'_',Fp.epochEnvironment], Fp.animals);end
if load_varContPwrCorr
    outdir = 'expvarCont'; outfile='expvarContCorr'; outpath=[pconf.andef{2},outdir,'/'];
    expvarContPwrCorr=load_data(outpath,[outfile,'_', Fp.epochEnvironment],Fp.animals);
    outdir = 'swrvarCont'; outfile='swrvarContCorr'; outpath=[pconf.andef{2},outdir,'/'];
    swrvarContPwrCorr=load_data(outpath,[outfile,'_', Fp.epochEnvironment],Fp.animals);
%     outdir = 'tfbvarCont'; outfile='tfbvarContCorr'; outpath=[pconf.andef{2},outdir,'/'];
%     tfbvarContPwrCorr=load_data(outpath,[outfile,'_', Fp.epochEnvironment],Fp.animals);
end
    % agg plot meanpower TFzmap /nt
if plot_meanpwr_timeXfreq; Pp = load_plotting_params({'defaults', 'power'});
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
                pconf = animaldef('Demetris');
                save_figure(sprintf('%s/expvarCatMeanPwr/',pconf{4}), 'expvarCatMeanPwr', sprtit)
                close all
            end
            close all;
            %                 end
        end
    end
end
    % plot regression TFzmap /nt
if plot_corrpwr_timeXfreq; Pp = load_plotting_params({'defaults', 'powerVarCorr'});
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        corrs = {'expvarContPwrCorr', 'swrvarContPwrCorr'};
        for icorr = 1:2
            PV = eval(corrs{icorr});
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
                pconf = animaldef('Demetris');
                save_figure(sprintf('%s/%s/',pconf{4}, corrs{icorr}), corrs{icorr}, sprtit)
                close all
            end
            close all;
        end
        end
    end   
end

    % agg plot meanpowerDiff TFzmap /nt


% plot ITPC TFzmap /nt
if plot_itpc_timeXfreq; Pp = load_plotting_params({'defaults', 'power'});
    for ian = 1:numel(Fp.animals) % for each animal
%         animal = Fp.animals{ian};
%         
%         ntinfanidx = find(strcmp({ntinfo.animal}, animal));
%         ntrodes = evaluatefilter(ntinfo(ntinfanidx).ntinfo, 'strcmp($valid, ''yes'')');
%         ntrodes = unique(ntrodes(:,3));
%         den = cellfetch(ntinfo(ntinfanidx).ntinfo, 'area');
%         matidx = unique(den.index(:,3));
%         anidx = find(strcmp({pwr.animal}, animal));
%         for co = 1:length(Fp.useripstates)
%             if savefigs && ~pausefigs
%                 close all
%                 ifig =figure('Visible','off','units','normalized','position', ...
%                     Pp.position);
%             else
%                 ifig = figure('units','normalized','position',Pp.position);
%             end
%             set(gcf,'color','white')
%             numcols = 5;
%             numrows = ceil(length(ntrodes) / numcols);
%             for nti = 1:length(ntrodes)
%                 sf = subaxis(numrows,numcols,nti, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
%                     'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
%                     Pp.MgTp, 'MarginBottom', Pp.MgBm);
%                 nt = ntrodes(nti);
%                 area = ntinfo(ntinfanidx).ntinfo{1}{1}{nt}.area;
%                 subarea = ntinfo(ntinfanidx).ntinfo{1}{1}{nt}.subarea;
%                 if isnumeric(subarea)
%                     subarea = num2str(subarea);
%                 end
%                 ntidx = find(matidx == nt);
%                 idata2plot = squeeze(itpc(anidx).ITPC{co}.ITPC(ntidx,:,:))';
%                 idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
%                                         'dsamp', itpc(anidx).wp.dsamp);
%                 time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
%                 contourf(sf, time, wp.frex, idata2plot, Pp.contourRes, ...
%                     'linecolor','none');
%                 set(gca,'ydir','normal','yscale','log');
%                 colormap(Pp.usecolormap)
%                 caxis(sf, 'auto')
% %                 colorbar
%                 hold on                
%                 ytickskip = 2:4:wp.numfrex;
%                 set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 7)
%                 title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
%                     'FontWeight',Pp.FontW, 'FontName', ...
%                     Pp.FontNm)
% %                 xlabel('time s', 'FontSize',8,'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
% %                 ylabel('freq Hz','FontSize',8, 'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
% %                 else 
% %                     set(gca, 'xlabel', [])
% %                     set(gca, 'ylabel', [])
% %                 end
% %                                 
%                 yl = ylim;
%                 line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
%             end
%             %% super
%             sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
%             sprtit = sprintf('%s %s itpc %s %s %s %s', animal, 'allnts', ...
%                 Fp.uselfptype, Fp.useripstates{co}, Fp.add_params{1}, Fp.add_params{2});
%             iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
%             set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
%                 'horizontalAlignment', 'center','FontSize', 16);
%             
%             %% ---- pause, save figs ----
%             if pausefigs
%                 pause
%             end
%             if savefigs
%                 pconf = animaldef('Demetris');
%                 save_figure(sprintf('%s/itpc/',pconf{4}), 'itpc', sprtit)
%                 close all
%             end
%             close all;
%             %                 end
%         end
    end
end
%% noise debugging plots

    % agg plot meanpower TFzmap /nt

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
% plot each ripple


%% load.make.plot final mean.dmean.itpc.ditpc.corr **stats** **1kperms**









