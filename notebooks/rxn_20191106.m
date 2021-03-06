%{
continuation of exampleBursts_20191101.m

JZ3 has no lick events?
how do i create the lick data structure?
- get_licks a la lickXcorrSU_20190916

%}

pconf = paramconfig;
create_filter = 0;
run_ff = 0;
load_ff = 0;
createSummaryData = 0;

plotfigs = 1;
plotFullperAn = 1;
plotFullAllAn = 0;

plotPPCperAn = 0;
plotPPCAllAn = 0;

% plotTrace = 0;
% plotPCdemo = 0;

showfigs = 0;
pausefigs = 0;
savefigs = 1;
savefigas = {'pdf', 'png'};

%% FF
Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'}; %, 'JZ1', 'JZ4'};
Fp.filtfunction = 'dfa_reactivationPSTH';
Fp.params = {'wtrackdays', 'excludePriorFirstWell', 'excludeAfterLastWell', ...
    Fp.filtfunction};
Fp = load_filter_params(Fp);
%%
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'cells', Fp.cellFilter, 'iterator',Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F), 'un', 1);
    F = runfilter(F);
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.env]);
end
if load_ff
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', ['_' Fp.env]);
end
%% Concat and Create results across epochs per animal
if createSummaryData
    figname = 'RxnPPCAllAn'; % just to load params
    Pp = load_plotting_params({'defaults',figname});
    time = Fp.win(1):Fp.bin:Fp.win(2);
    time = time(1:end-1)+Fp.bin/2;
    winidx = knnsearch(time', Pp.winSE');
    time = time(winidx(1):winidx(2))';
    
    % PER ANIMAL =======================================
    D = struct;
    for a = 1:length(F)
        D(a).animal = F(a).animal;
        % PER ANIMAL FULL MODEL =============
        RswrIntraMean = cell2mat({F(a).output{1}.swrIntraBurstRXNZ_mean}');
        D(a).RswrIntraMeanDay = RswrIntraMean(:,winidx(1):winidx(2));  % clip win
        D(a).RswrIntraMeanAn = nanmean(D(a).RswrIntraMeanDay); % mean DayMEAN across days
        RswrIntraSem = cell2mat({F(a).output{1}.swrIntraBurstRXNZ_sem}');
        D(a).RswrIntraSemDay = RswrIntraSem(:,winidx(1):winidx(2)); % clip win
        D(a).RswrIntraSemAn = nanmean(D(a).RswrIntraSemDay); % mean DaySEM across days
        
        RswrExtraMean = cell2mat({F(a).output{1}.swrExtraBurstRXNZ_mean}');
        D(a).RswrExtraMeanDay = RswrExtraMean(:,winidx(1):winidx(2));
        D(a).RswrExtraMeanAn = nanmean(D(a).RswrExtraMeanDay);
        RswrExtraSem = cell2mat({F(a).output{1}.swrExtraBurstRXNZ_sem}');
        D(a).RswrExtraSemDay = RswrExtraSem(:,winidx(1):winidx(2));
        D(a).RswrExtraSemAn = nanmean(D(a).RswrExtraSemDay);
        
        RxpMean = cell2mat({F(a).output{1}.XPnoswrRXNZ_mean}');
        D(a).RxpMeanDay = RxpMean(:,winidx(1):winidx(2));
        D(a).RxpMeanAn = nanmean(D(a).RxpMeanDay);
        RxpSem = cell2mat({F(a).output{1}.XPnoswrRXNZ_sem}');
        D(a).RxpSemDay = RxpSem(:,winidx(1):winidx(2));
        D(a).RxpSemAn = nanmean(D(a).RxpSemDay);
        
        % PER ANIMAL FULL MODEL ZSCORE ======
        
        
        % PER ANIMAL PERPC MODEL ==========================
        try
            for d = 1:length(F(a).output{1})
                D(a).swrIntraRXNperPC_meanDay(:,:,d) = ...
                    F(a).output{1}(d).swrIntraBurstRXNperPC_mean(1:Pp.numPCs,:);
                D(a).swrExtraRXNperPC_meanDay(:,:,d) = ...
                    F(a).output{1}(d).swrExtraBurstRXNperPC_mean(1:Pp.numPCs,:);
                D(a).xpRXNperPC_meanDay(:,:,d) = ...
                    F(a).output{1}(d).xpRXNperPC_mean(1:Pp.numPCs,:);
            end
            % PER ANIMAL PERPC MODEL ZSCORE ======
            % i don't think i'm doing the zscore corretly here
            D(a).swrIntraRXNperPC_meanAnZ = zscore(squeeze(nanmean(D(a).swrIntraRXNperPC_meanDay, 3))');
            D(a).swrExtraRXNperPC_meanAnZ = zscore(squeeze(nanmean(D(a).swrExtraRXNperPC_meanDay, 3))');
            D(a).xpRXNperPC_meanAnZ = zscore(squeeze(nanmean(D(a).xpRXNperPC_meanDay, 3))');
        catch
            fprintf(' ***** no PPC MODEL?\n')
        end
    end
    
    % COMBINE ACROSS ANIMALS ==========================
    allAn = struct;
    % FULL MODEL =============
    allAn.RswrIntraMeanAllAn = nanmean(cell2mat({D.RswrIntraMeanAn}')); 
    allAn.RswrIntraSemAllAn = sem(cell2mat({D.RswrIntraMeanAn}'),1); % sem of animal means
    allAn.RswrExtraMeanAllAn = nanmean(cell2mat({D.RswrExtraMeanAn}')); 
    allAn.RswrExtraSemAllAn = sem(cell2mat({D.RswrExtraMeanAn}'),1); % sem of animal means
    allAn.RxpMeanAllAn = nanmean(cell2mat({D.RxpMeanAn}')); 
    allAn.RxpMeanAllAn = sem(cell2mat({D.RxpMeanAn}'),1); 
    
    % PERPC MODEL =============
    try
        allAn.swrIntraRXNperPC_meanAllAnZ = nanmean(cell2mat(permute({D.swrIntraRXNperPC_meanAnZ},[1 3 2])),3);
        allAn.swrExtraRXNperPC_meanAllAnZ  = nanmean(cell2mat(permute({D.swrExtraRXNperPC_meanAnZ},[1 3 2])),3);
        allAn.xpIntraRXNperPC_meanAllAnZ = nanmean(cell2mat(permute({D.xpRXNperPC_meanAnZ},[1 3 2])),3);
    catch
        
    end
end

%% ------------------------------plot------------------------------
if plotfigs
    %% plot perAnimal allEventGroups rxn ETA full
    if plotFullperAn
        figname = 'wRxnFullperAn';
        Pp=load_plotting_params({'defaults',figname});
        for a = 1:length(D)
            %% Fig start
            ifig = init_plot(showfigs, Pp.position); % init fig
            time = Fp.win(1):Fp.bin:Fp.win(2);
            time = time(1:end-1)+Fp.bin/2;
            winidx = knnsearch(time', Pp.winSE');
            time = time(winidx(1):winidx(2))';
            
            %% plot all swr rxn ETA full

            sf = subaxis(1,1,1,Pp.posparams{:});           
            Rm = D(a).RswrIntraMeanAn;
            Rsem = D(a).RswrIntraSemAn;
            f1 = plot(time, Rm, 'b', 'linewidth', 1, 'DisplayName','SWR inBurst');
            hold on;
            fill([time; flipud(time)],[Rm'+Rsem'; flipud(Rm'-Rsem')],'b', ...
                'linestyle','none', 'facealpha', .1)
            
            Rm = D(a).RswrExtraMeanAn;
            Rsem = D(a).RswrExtraSemAn;
            f2 = plot(time, Rm, 'k', 'linewidth', 1, 'DisplayName','SWR exBurst');
            fill([time; flipud(time)],[Rm'+Rsem'; flipud(Rm'-Rsem')],'k', ...
                'linestyle','none', 'facealpha', .1)
            
            Rm = D(a).RxpMeanAn;
            Rsem = D(a).RxpSemAn;
            f3 = plot(time, Rm, 'r', 'linewidth', 1, 'DisplayName','XP inBurst');
            fill([time; flipud(time)],[Rm'+Rsem'; flipud(Rm'-Rsem')],'r', ...
                'linestyle','none', 'facealpha', .1)
            hold off
            %Pp
            
            axis tight
            ylabel('rxn strength')
            xlabel('Time from Event Start (s)')
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
            legend([f1 f2 f3]);
%             title('Full Model RXN')

            %% End Fig
            try
                animal = F(a).animal{3};
            catch
                animal = F(a).animal;
            end
            stit = sprintf('%s %s %s', animal, figname);
            setSuperAxTitle(stit);
            
            if pausefigs; pause; end
            if savefigs
                save_figure(strjoin({pconf.andef{4}, figname}, '/'), stit,...
                    'savefigas', savefigas);
            end
        end
    end
end
    if plotFullAllAn
    figname = 'wRxnFullAllAn';
    Pp=load_plotting_params({'defaults',figname});
    ifig = init_plot(showfigs, Pp.position); % init fig
    time = Fp.win(1):Fp.bin:Fp.win(2);
    time = time(1:end-1)+Fp.bin/2;
    winidx = knnsearch(time', Pp.winSE');
    time = time(winidx(1):winidx(2))';
        
        %% plot all swr rxn ETA full
        sf = subaxis(1,1,1,Pp.posparams{:});
        Rm = allAn.RswrIntraMeanAllAn;
        Rsem = allAn.RswrIntraSemAllAn;
        f1 = plot(time,Rm,'b','linewidth', 1, 'DisplayName','SWR inBurst');
        hold on;
        fill([time; flipud(time)],[Rm'+Rsem'; flipud(Rm'-Rsem')],'b', ...
            'linestyle','none', 'facealpha', .1)
        
        Rm = allAn.RswrExtraMeanAllAn;
        Rsem = allAn.RswrExtraSemAllAn;
        f2 = plot(time, Rm, 'k', 'linewidth', 1, 'DisplayName','SWR exBurst');
        fill([time; flipud(time)],[Rm'+Rsem'; flipud(Rm'-Rsem')],'k', ...
            'linestyle','none', 'facealpha', .1)
        
        Rm = allAn.RxpMeanAllAn;
        Rsem = allAn.RxpMeanAllAn;
        f3 = plot(time, Rm, 'r', 'linewidth', 1, 'DisplayName','XP inBurst');
        fill([time; flipud(time)],[Rm'+Rsem'; flipud(Rm'-Rsem')],'r', ...
            'linestyle','none', 'facealpha', .1)
        hold off
        
        %Pp
        axis tight
        ylabel('rxn strength')
        xlabel('Time from Event Start (s)')
        line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
        legend([f1 f2 f3]);
        
        %% End Fig
        stit = sprintf('%s', figname);
        setSuperAxTitle(stit);
        
        if pausefigs; pause; end
        if savefigs
            save_figure(strjoin({pconf.andef{4}, figname}, '/'), stit, ...
            'savefigas', savefigas);
        end
    end
    
    %% perPC perAn RXN
    if plotPPCperAn
        figname = 'RxnPPCperAn';
        Pp=load_plotting_params({'defaults',figname});
        time = Fp.win(1):Fp.bin:Fp.win(2);
        time = time(1:end-1)+Fp.bin/2;
        winidx = knnsearch(time', Pp.winSE');
        time = time(winidx(1):winidx(2))';
        for a = 1:length(F) % animal
%             F(a).output{1}.eigValSortSig
            c = lines(Pp.numPCs);
%             c = magma(3);
            ifig = init_plot(showfigs, Pp.position); % init fig
            %% intraBurst SWR perPC perAn
            sf = subaxis(3,1,1,Pp.posparams{:});
            pcRs = D(a).swrIntraRXNperPC_meanAnZ(winidx(1):winidx(2),:);
            hold on;
            arrayfun(@(x) plot(time, pcRs(:,x), 'color', c(x,:), 'linewidth', 1), 1:Pp.numPCs, 'un', 0);
            %Pp
            axis tight
            hold off;
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', 1)
            ylabel('SWR inBurst')

            %% extraBurst SWR perPC perAn
            sf = subaxis(3,1,2,Pp.posparams{:});
            pcRs = D(a).swrExtraRXNperPC_meanAnZ(winidx(1):winidx(2),:);
            hold on;
            arrayfun(@(x) plot(time, pcRs(:,x), 'color', c(x,:), 'linewidth', 1), 1:Pp.numPCs, 'un', 0);
            %Pp
            axis tight
            hold off;
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', 1)
            ylabel('SWR exBurst')
            
            %% XP noswr perPC perAn
            sf = subaxis(3,1,3,Pp.posparams{:});
            pcRs = D(a).xpRXNperPC_meanAnZ(winidx(1):winidx(2),:);
            hold on;
            arrayfun(@(x) plot(time, pcRs(:,x), 'color', c(x,:), 'linewidth', 1), 1:Pp.numPCs, 'un', 0);
            %Pp
            axis tight
            hold off;
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', 1)
            ylabel('XP noswr')
            
            %% End Fig
            allAxesInFigure = findall(gcf,'type','axes');
            linkaxes(allAxesInFigure, 'y');
            try
                animal = D(a).animal{3};
            catch
                animal = D(a).animal;
            end
            stit = sprintf('%s %s %s', figname, animal);
            setSuperAxTitle(stit);
            
            if pausefigs; pause; end
            if savefigs
                save_figure(strjoin({pconf.andef{4}, figname}, '/'), stit, 'savefigas', savefigas);
            end
        end
    end
    
%% perPC All An RXN
    if plotPPCAllAn
        figname = 'RxnPPCAllAn';
        Pp=load_plotting_params({'defaults',figname});
        time = Fp.win(1):Fp.bin:Fp.win(2);
        time = time(1:end-1)+Fp.bin/2;
        winidx = knnsearch(time', Pp.winSE');
        time = time(winidx(1):winidx(2))';
        for a = 1:length(F) % animal
%             F(a).output{1}.eigValSortSig
            c = lines(Pp.numPCs);
%             c = magma(3);
            ifig = init_plot(showfigs, Pp.position); % init fig
            %% intraBurst SWR perPC perAn
            sf = subaxis(3,1,1,Pp.posparams{:});
            pcRs = allAn.swrIntraRXNperPC_meanAllAnZ(winidx(1):winidx(2),1:Pp.numPCs);
            hold on;
            arrayfun(@(x) plot(time, pcRs(:,x), 'color', c(x,:), 'linewidth', 1), 1:Pp.numPCs, 'un', 0);
            %Pp
            axis tight
            hold off;
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', 1)
            ylabel('SWR inBurst')

            %% extraBurst SWR perPC perAn
            sf = subaxis(3,1,2,Pp.posparams{:});
            pcRs = D(a).swrExtraRXNperPC_meanAnZ(winidx(1):winidx(2),1:Pp.numPCs);
            hold on;
            arrayfun(@(x) plot(time, pcRs(:,x), 'color', c(x,:), 'linewidth', 1), 1:Pp.numPCs, 'un', 0);
            %Pp
            axis tight
            hold off;
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', 1)
            ylabel('SWR exBurst')
            
            %% XP noswr perPC perAn
            sf = subaxis(3,1,3,Pp.posparams{:});
            pcRs = D(a).xpRXNperPC_meanAnZ(winidx(1):winidx(2),1:Pp.numPCs);
            hold on;
            arrayfun(@(x) plot(time, pcRs(:,x), 'color', c(x,:), 'linewidth', 1), 1:Pp.numPCs, 'un', 0);
            %Pp
            axis tight
            hold off;
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', 1)
            ylabel('XP noswr')
            
            %% End Fig
            allAxesInFigure = findall(gcf,'type','axes');
            linkaxes(allAxesInFigure, 'y');
            try
                animal = D(a).animal{3};
            catch
                animal = D(a).animal;
            end
            stit = sprintf('%s %s %s', figname, animal);
            setSuperAxTitle(stit);
            
            if pausefigs; pause; end
            if savefigs
                save_figure(strjoin({pconf.andef{4}, figname}, '/'), stit, 'savefigas', savefigas);
            end
        end
    end
    
    

%% plot spatial correspondance to the PC's perEpoch

%% plot perAnimal top 3 PC mean for swr in burst vs licks


