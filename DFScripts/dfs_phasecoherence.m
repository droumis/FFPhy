


%% Dashboard
%%% First create the processed data structure similar to dfs_riptriglfp.m
close all
runFilterFramework = 0;
; saveFilterOutput = runFilterFramework;
; loadFilterOutput = 0;

%%% Then run phasecoherence on the LFP traces

calcEventState = 0;
; saveEventState = calcEventState;
; loadEventState = 0;
calculateAnalyticSignal = 0;
; saveAnalyticSignal = calculateAnalyticSignal;
; loadAnalyticSignal = 0;
calculateIXPC = 0;
; saveIXPC = calculateIXPC;
% ; combineByArea = 0;%calculateIXPC;
; loadIXPCresults = 0;
runPermutationTesting = 0;
; savePermutationOutput = runPermutationTesting;

plotIXPC = 1;
plotISPC = 0;
savefigs = 1;
pausefigs = 0;
calcfunction = 'ISPC'; %ITPC power ISPC
% calcfunction = 'ITPC'; % 'power' or 'ITPC'
plotzmask = 0;
plotlogbased = 1; % if zero, plots percent baseline normalized
plotByArea = 1;
plotByNTrode = 0;

%% ----------------  preset params --------------------------
behavestruct = 'BehaveState';
colorSet = 'DR1';
figspecs = 'itpc2';
areafigspecs = 'AreasByDays'; %'itpcAreas'
ripSet = 'DR4'; % DR4 = .5 excl dur; 3std
waveSet = 'DR1';
indsSet = 'performanceByDay'; %'DR1'; byDay; performance; performanceByDay
% clims = [-8 8]; %zscore power
% coloraxis = 'auto';%[-10 10];

%% ---------------- Data Filters --------------------------
animals = {'JZ1'}; %{'JZ1'}; 
days = [];
filtfunction = 'riptriglfp';
LFPtypes = {'eeg'};%, 'ripple', 'theta', 'lowgamma', 'fastgamma'}; %
eventtype = 'rippleskons';
epochType = {'run'}; %sleep or run
epochEnvironment = {'wtrack'};% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
eventSourceArea = 'ca1';
ntAreas = {'ca1','sub','mec', 'por', 'v2l', 'ref'}; %which areas to include in datafilrer
rippms = getRipParams(ripSet); % set parameters for ripple detection like speed, numtetrodes, etc.

%% ----------------- wavelet parameters ------------------------------------
win = [-1.5 1.5]; %in seconds. The morlet wavelet needs to be centered at time zero.. so always make this window symmetric about zero

%% ---------------- Paths and Title strings ---------------------------------------------------
investInfo = animaldef(lower('Demetris'));
filtOutputDirectory = sprintf('%s%s/', investInfo{2}, filtfunction);
figdirectory = sprintf('%s%s/%s/', investInfo{4}, calcfunction);
% filenamesave = sprintf('%s%s_%s_%s_%s', eventSourceArea, eventtype, epochEnvironment, cell2mat(animals), cell2mat(LFPtypes));
filenamesave = sprintf('%s%sSTD%d_%s_%s_%s', eventSourceArea, eventtype, rippms.minstdthresh, strjoin(epochEnvironment,'-'),...
    strjoin(animals,'-'), strjoin(LFPtypes, '-'));%strjoin(arrayfun(@(x) num2str(x,'%-2d'),days,'UniformOutput',false),'-'));
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
resultfilename = filenamesave; %sprintf('%s_%s%s_%s_%s_D%s', calcfunction, eventSourceArea, eventtype, strjoin(epochEnvironment,'-'), cell2mat(animals), strrep(num2str(days, '%-2d'),' ', '-'));
filenameTitle = strrep(resultfilename,'_', ' ');
DataTypeFields = {'all', 'corrOut', 'mistOut', 'outB', 'inB','outB-inB', 'corrOut-mistOut'};


%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1;
    eptypeEnv = [epochType; epochEnvironment];
    epochfilter = sprintf('((isequal($type, ''%s'')) && (isequal($environment, ''%s''))) || ',eptypeEnv{:});
    yuck = strfind(epochfilter, '||');
    epochfilter = epochfilter(1:yuck(end)-1);  %cut off the trailing '||'
    iterator = 'multitetrodeanal'; %multitetrodeanal
    %tetfilter: the ntrodeID's that pass this filter get stored into f.eegdata{day}{epoch}[ntrodeIDs]
    tetfilter = sprintf('(isequal($area,''%s'')) || ', ntAreas{:});
    yuck = strfind(tetfilter, '||');
    tetfilter = tetfilter(1:yuck(end)-1); %cut off the trailing '||'
    timefilter{1} = {'getconstimes', '($cons == 1)', [eventSourceArea eventtype],1,'consensus_numtets',rippms.consensus_numtets,...
        'minstdthresh',rippms.minstdthresh,'exclusion_dur',rippms.exclusion_dur,'minvelocity',rippms.minvelocity,'maxvelocity',rippms.maxvelocity};
    %     timefilter{2} = {'excludenoiseevents', '($noise == 0)', [eventSourceArea,'noisekons'], 1, };
    
    %----------F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);--------
    F = createfilter('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);
    %----------f = setfilteriterator(f, funcname, loadvariables, options)--------
%     eval(['F = setfilterfunction(F, [''dfa_'' filtfunction], {[eventSourceArea eventtype],' strjoin(arrayfun(@(x) sprintf('''%s'',', cell2mat(x)), LFPtypes,'UniformOutput',false)) '},''eventtype'', [eventSourceArea eventtype], ''LFPtypes'', LFPtypes, ''win'', win);']);
    F = setfilterfunction(F, sprintf('dfa_%s', filtfunction), {[eventSourceArea eventtype], strjoin(arrayfun(@(x) sprintf('%s', cell2mat(x)), LFPtypes,'UniformOutput',false))},'eventtype', [eventSourceArea eventtype], 'LFPtypes', LFPtypes, 'win', win);
    %     eval(['F = setfilterfunction(F, [''dfa_'' filtfunction], {[eventSourceArea eventtype],' strjoin(reshape(repmat(arrayfun(@(x) sprintf('''%s'',', cell2mat(x)), LFPtypes,'UniformOutput',false), 2, 1), 1,length(LFPtypes)*2)) '},''eventtype'', [eventSourceArea eventtype], ''LFPtypes'', LFPtypes);']);
    %     tic
    F = runfilter(F);
    %     F(1).filterTimer = toc; F(1).filterTimer
    %     F(1).worldDateTime = clock;
    %     F(1).dataFilter = struct('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator, 'filename', filename);
end
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    if ~isdir(filtOutputDirectory);
        mkdir(filtOutputDirectory);
    end
    %save the entire workspace for filter provenance
    save(sprintf('%s/%s',filtOutputDirectory, filename), 'F','-v7.3');
    disp(sprintf('filteroutput saved to %s/%s',filtOutputDirectory, filename))
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s/%s',filtOutputDirectory, filename))
    disp(sprintf('filteroutput loaded: %s/%s',filtOutputDirectory, filename))
end
%% ---------------- Calc Event State ---------------------------------------------------
if calcEventState
    tic
    [eState] = calcEventPerformanceState(F, animals, behavestruct);
    %% ---------------- Save eventState Output ---------------------------------------------------
    if saveEventState == 1;
        dtypedir = 'eventState';
        if ~isdir(sprintf('%s%s/', investInfo{3}, dtypedir));
            mkdir(sprintf('%s%s/', investInfo{3}, dtypedir));
        end
%         savestr = sprintf('%s%s/%s_%dHz-%dHz_rCyc%.1f-%.1f.mat', investInfo{3}, calcfunction, resultfilename, min_freq, max_freq,range_cycles(1), range_cycles(2));
        savestr = sprintf('%s%s/%s_waveSet-%s_eventState.mat', investInfo{3}, dtypedir, resultfilename, waveSet);
        save(savestr, 'eState','-v7.3');
        disp(sprintf('SAVED EVENTSTATE RESULTS ++++++++++ %s',savestr))
    end
    eventtoc = toc;
end
if loadEventState == 1;
    dtypedir = 'eventState';
%     load(sprintf('%s%s/%s_%dHz-%dHz_rCyc%.1f-%.1f.mat', investInfo{3}, calcfunction, resultfilename, min_freq, max_freq,range_cycles(1), range_cycles(2)));
    load(sprintf('%s%s/%s_waveSet-%s_eventState.mat', investInfo{3}, dtypedir, resultfilename, waveSet));
end

%% Calculate Analytic Signal
if calculateAnalyticSignal
    tic
    numanimals = length(F); 
    for ian = 1:length(numanimals)
        ixpc.dataByDay = []; ixpc.dataByDay{ian} = [];
        ixpc.animals{ian} = F(ian).animal{1};
        ixpc.andays{ian} = find(~cellfun(@isempty,F(ian).output)); %get nonempty eps
        [ixpc.index{ian}, allNTDataCat, ixpc.dataByDay{ian}] = gatherRipSnips(...
            ixpc.andays{ian}, F(ian).output);
%         [ixpc.IndTypes{ian}, ixpc.datatypes{ian}, ixpc.datatypesMat{ian}] = getIndsByType(indsSet,eState, ian, 'dataByDay', ixpc.dataByDay{ian});
        ixpc.wp = getWaveParams(waveSet, allNTDataCat);
        computeAnalyticSignal(allNTDataCat, ixpc.wp,ixpc.animals{ian}, resultfilename, 'saveAnalyticSignal', saveAnalyticSignal);
        clear allNTDataCat F
    end
    AStoc = toc;
end
%% ---------------- LOAD Analytic Signal Output ---------------------------------------------------
if loadAnalyticSignal == 1;
    dtypedir = 'analyticSignal';
    loadstr = sprintf('%s%s/%s_waveSet-%s_analyticSignal.mat', investInfo{3}, dtypedir, resultfilename, waveSet);
    load(loadstr);
end
%% Calculate inter-trial/site phase clustering with morlet wavelet convolution
if calculateIXPC
    tic
    dtypedir = 'analyticSignal';
    waveSet = ixpc.wp.waveSet;
    baseind = ixpc.wp.baseind;
    for ian = 1:length(ixpc.animals)
        animal = ixpc.animals{ian};
        [ixpc.IndTypes{ian}, ixpc.datatypes{ian}, ixpc.datatypesMat{ian}] = getIndsByType(indsSet,eState, ian, 'dataByDay', ixpc.dataByDay{ian});
        IndTypes = ixpc.IndTypes{ian};
        nNTrodes = ixpc.wp.nNTrodes;
        
        switch calcfunction
            case 'ITPC'
                for nt = 1:nNTrodes %use parfor
                    [ph] = loadAS(animal, nt, waveSet, 'PH'); ph = ph.ph; %loads as 'ph.ph'
                    [baseLmeanITPC{nt}, pctbasedITPCout{nt}, logbasedITPCout{nt}] = ...
                        computeIXPC(ph,animal,nt, sprintf('%d', nt),IndTypes,baseind,indsSet, 'ixpcType', calcfunction, 'numNTs', nNTrodes);
                end
                ixpc.baseLmeanITPC{ian} = baseLmeanITPC;
                ixpc.pctbasedITPCout{ian} = pctbasedITPCout;
                ixpc.logbasedITPCout{ian} = logbasedITPCout;
                [ixpc.areas,ixpc.logbasedareaITPCmean{ian}] = combineDataByArea(ixpc.index, animals, ian, ixpc.logbasedITPCout{ian});
                [ixpc.areas,ixpc.pctbasedareaITPCmean{ian}] = combineDataByArea(ixpc.index, animals, ian, ixpc.pctbasedITPCout{ian});
            case 'power'
                for nt = 1:nNTrodes %use parfor
                    [as] = loadAS(animal, nt, waveSet, 'AS'); as = as.as; %loads as 'ph.ph'
                    [baseLmeanpower{nt}, pctbasedpowerout{nt}, logbasedpowerout{nt}] = ...
                        computePower(as,animal,nt, sprintf('%d', nt),IndTypes,baseind,indsSet);
                end
                ixpc.baseLmeanpower{ian} = baseLmeanpower;
                ixpc.pctbasedpowerout{ian} = pctbasedpowerout;
                ixpc.logbasedpowerout{ian} = logbasedpowerout;
                [ixpc.areas,ixpc.logbasedareapowermean{ian}] = combineDataByArea(ixpc.index, animals, ian, ixpc.logbasedpowerout{ian});
                [ixpc.areas,ixpc.pctbasedareapowermean{ian}] = combineDataByArea(ixpc.index, animals, ian, ixpc.pctbasedpowerout{ian});
            case 'ISPC'
                if isstruct(ph)
                    for nt = 1:nNTrodes %use parfor
                        %                     clear ph
                        pause
                        tmpph = loadAS(animal, nt, waveSet, 'PH');
                        ph{nt} = tmpph.ph; %loads as 'ph.ph'
                    end
                end
                ixpc.ntpairs{ian} = nchoosek([1:nNTrodes],2);
                clear baseLmeanISPC pctbasedISPCout logbasedISPCout
                for ipair = 1:size(ntpairs,1);
                    ntA = ntpairs(ipair,1);
                    ntB = ntpairs(ipair,2);
                    phDiff = bsxfun(@minus, ph{ntA}, ph{ntB});
                    %ISPC is running ITPC on the phase difference between sites
                    [baseLmeanISPC{ipair}, pctbasedISPCout{ipair}, logbasedISPCout{ipair}] = ...
                    computeIXPC(phDiff,animal,ipair,sprintf('pair %d - %d ', ntA, ntB),IndTypes,baseind,indsSet, 'ixpcType', calcfunction, 'numNTs', size(ntpairs,1));
                end
                ixpc.baseLmeanISPC{ian} = baseLmeanISPC;
                ixpc.pctbasedISPCout{ian} = pctbasedISPCout;
                ixpc.logbasedISPCout{ian} = logbasedISPCout;
                [ixpc.areas{ian},ixpc.logbasedareaISPCmean{ian}] = combineDataByArea(ixpc.ntpairs{ian}, animals, ian, ixpc.logbasedISPCout{ian}, 'ntindcolumn', [1 2]);
                [ixpc.areas{ian},ixpc.pctbasedareaISPCmean{ian}] = combineDataByArea(ixpc.ntpairs{ian}, animals, ian, ixpc.pctbasedISPCout{ian}, 'ntindcolumn', [1 2]);
        end
    end
    %% ---------------- Save RESULTS Output ---------------------------------------------------
    if saveIXPC == 1;
        dtypedir = calcfunction;
        if ~isdir(sprintf('%s%s/', investInfo{3}, dtypedir));
            mkdir(sprintf('%s%s/', investInfo{3}, dtypedir));
        end
        %     savestr = sprintf('%s%s/%s_%dHz-%dHz_rCyc%.1f-%.1f_%s.mat', investInfo{3}, calcfunction, resultfilename, min_freq, max_freq,range_cycles(1), range_cycles(2), indsSet);
        savestr = sprintf('%s%s/%s_waveSet-%s_indsSet-%s.mat', investInfo{3}, dtypedir, resultfilename, waveSet, indsSet);
        save(savestr, 'ixpc','-v7.3');
        disp(sprintf('SAVED RESULTS ++++++++++ %s',savestr))
    end
    ixpctoc = toc;
end
%% ---------------- Load Results Output ---------------------------------------------------
if loadIXPCresults == 1;
    dtypedir = calcfunction;
    loadstr = sprintf('%s%s/%s_waveSet-%s_indsSet-%s.mat', investInfo{3}, dtypedir, resultfilename, waveSet, indsSet);
    load(loadstr);
end
%% permutation test % need to fix
if runPermutationTesting
    tic
    %% permutation testing... this will take TIME and a shit ton of RAM.. 
    permTestType = 'outboundVSinbound';
    ixpc = permutationTest(ixpc, as, ph, ian, intr, IndTypes, DataTypeFields, n_permutes, permTestType);
    permTestType = 'correctOutvsmistakeOut';
    ixpc = permutationTest(ixpc, as, ph, ian, intr, IndTypes, DataTypeFields, n_permutes, permTestType);
    %% ---------------- Save RESULTS Output ---------------------------------------------------
    if savePermutationOutput == 1;
        dtypedir = sprintf('%s/%s_permutes',calcfunction, calcfunction);
        if ~isdir(sprintf('%s%s/', investInfo{3}, dtypedir));
            mkdir(sprintf('%s%s/', investInfo{3}, dtypedir));
        end
        %     savestr = sprintf('%s%s/%s_%dHz-%dHz_rCyc%.1f-%.1f_%s.mat', investInfo{3}, calcfunction, resultfilename, min_freq, max_freq,range_cycles(1), range_cycles(2), indsSet);
        savestr = sprintf('%s%s/%s_waveSet-%s_indsSet-%s.mat', investInfo{3}, dtypedir, resultfilename, waveSet, indsSet);
        save(savestr, 'ixpc','-v7.3');
        disp(sprintf('SAVED RESULTS ++++++++++ %s',savestr))
    end
    permtoc = toc;
end
%% ---------------- plot ITPC---------------------------------------------------------------------------------------------
if plotIXPC
    if plotByNTrode % need to fix
        error('need to fix')
        warning('off', 'MATLAB:contour:ConstantData')
        fig = getFigspecs(figspecs);
        for ian = 1:length(animals)
            [nt] = getNTinfo(ixpc.index, animals, ian, colorSet);
            animalID = ixpc.animals{ian};
            for iDT = 1:(length(ixpc.datatypes))%+2
                [ifig, sfrows, sfcols] = prepFig(savefigs, pausefigs, nt, fig);
                %% ---- loop across all tets for this day ----
                for nt = 1:nt.nNTrodes;
                    introdeID = nt.ntrodesIndices(nt.numsumSortInds(nt),3);
                    isupareatag = nt.ntrodeTags{nt.numsumSortInds(nt)}.suparea;
                    iareatag = nt.ntrodeTags{nt.numsumSortInds(nt)}.area;
                    isubareatag = nt.ntrodeTags{nt.numsumSortInds(nt)}.subarea;
                    iNTcolor = nt.icolors(nt,:);
                    intfig = subaxis(sfrows,sfcols,nt, 'SpacingVertical', fig.SpacingVertical, 'SpacingHorizontal', fig.SpacingHorizontal, ...
                        'Padding', fig.Padding, 'MarginLeft', fig.MarginLeft, 'MarginRight', fig.MarginRight, 'MarginTop', fig.MarginTop, ...
                        'MarginBottom', fig.MarginBottom);
                    %% ---- Main Plotting ----
                    if plotlogbased
                        eval(sprintf('idata2plot = squeeze(ixpc.logbased%sout{ian}{nt.numsumSortInds(intr)}{iar}(:,iDT,:))'';', calcfunction));
                    else
                        eval(sprintf('idata2plot = squeeze(ixpc.percbased%sout{ian}{nt.numsumSortInds(intr)}{iar}(:,iDT,:))'';', calcfunction));
%                         eval(sprintf('idata2plot = squeeze(ixpc.%sout{ian}{nt.numsumSortInds(intr)}{iar}(:,iDT,:))'';', plotoutputtype));
                    end
                    idata2plot = trim2win(idata2plot, ixpc.wp.srate, ixpc.wp.plotwin);
%                     idata2plotNotBased = trim2win(idata2plotNotBased, wp.srate, wp.plotwin);
                    try
                        [~,bn] = contourf(intfig, plottimeWin,frex(1:end-1),idata2plot,fig.contourRes,'linecolor','none'); %
                        set(gca,'ydir','normal','xlim',[plotwin(1) plotwin(2)], 'ylim',[frex(1) frex(end-1)])
                    catch
                        [~,bn] = contourf(intfig, ixpc.wp.plottimeWin,ixpc.wp.frex,idata2plot,fig.contourRes,'linecolor','none'); %
                        set(gca,'ydir','normal','xlim',[ixpc.wp.plotwin(1) ixpc.wp.plotwin(2)], 'ylim',[ixpc.wp.min_freq ixpc.wp.max_freq])
                    end
                    caxis(coloraxis)
                    colormap(fig.usecolormap)
                    hold on;
                    if (strcmp(ixpc.datatypes{ian}{iDT}, 'outB-inB') || strcmp(ixpc.datatypes{ian}{iDT}, 'corrOut-mistOut')) && plotzmask
%                     if (iDT == 6 || iDT == 7) && plotzmask % plot the significance countours
                        eval(sprintf('zmask2plot = squeeze(ixpc.%szmask{ian}{nt.numsumSortInds(int)}{iDT})'';', calcfunction));
                        eval(sprintf('MCminmax = abs(ixpc.MC_%s_minmax{ian}{nt.numsumSortInds(int)}{iDT})'';',calcfunction));
                        eval(sprintf('irawdata = squeeze(ixpc.%sout{ian}{nt.numsumSortInds(int)}(:,iDT,:))'';', calcfunction));
                        zmask2plot = trim2win(zmask2plot, srate, plotwin);
                        irawdata = trim2win(irawdata, srate, plotwin);
                        zmask = zmask2plot;
                        zmask(abs(zmask2plot)<zval) = 0;
                        %                     [~,h] = contour(intfig,plottimeWin,frex(1:end-1),logical(zmask),1);
                        %                     h.LineColor = [.85 .85 .85];
                        hold on;
                        MCthresh = MCminmax(ceil(length(MCminmax)*(1-pval)));
                        MCthreshmap = irawdata;
                        MCthreshmap(abs(irawdata)<MCthresh) = 0;
                        [~,mc] = contour(intfig,plottimeWin,frex,logical(MCthreshmap),1);
                        mc.LineColor = fig.mcLineColor;%[.85 .85 .85]; %[.6 .6 .6];
                    end
                    %                 set(gca,'ydir','normal','xlim',[plotwin(1) plotwin(2)], 'ylim',[frex(1) frex(end)])
                    if mod(nt, sfcols) ~= 1;
                        set(gca, 'yscale','log','ytick',[])
                        %                                         set(gca, 'yscale','log','ytick',ceil(logspace(log10(min(frex)),log10(max(frex)),6)))
                        set(gca, 'FontSize',8, 'FontName', 'Arial');%
                    else
                        set(gca, 'yscale','log','ytick',ceil(logspace(log10(min(frex)),log10(max(frex)),6)))
                        %                     set(gca, 'YTick',[round(linspace(min(frex),max(frex),6))])
                        set(gca, 'FontSize',8, 'FontName', 'Arial');%
                    end
                    if nt <= (sfrows-1)*sfcols;
                        set(gca, 'XTick', []);
                    else
                        set(gca, 'XTick',[plotwin(1):plotwin(2)/2:plotwin(2)],'FontSize',8, 'FontName', 'Arial');
                    end
                    %% ---- Source ripple line, subplot title ----
                    Xline = [0 0];
                    Yline = [min_freq max_freq];
                    line(Xline, Yline, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                    ititle = title(sprintf('nT%d %s %s',introdeID, iareatag, num2str(isubareatag)));
                    set(ititle,'FontSize',10,'Color', iNTcolor, 'FontName', 'Arial','FontWeight','bold');
                end
                %% ---- super title and colorbar----
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                xlab = 'Time (s)';
                ylab = 'Frequency (Hz)';
                supxlabel = text(.5, .02, xlab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'Parent', sprtitleax,...
                    'Units', 'normalized', 'horizontalAlignment', 'center');
                supylabel = text(.01, .5, ylab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'rotation', 90, ...
                    'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
                if iDT == 6 || iDT == 7
                    %                 caxis([-2 2]);
                else
                    %                 caxis([-2 2]);
                end
                
                clrbar = colorbar('location','eastoutside', 'FontSize',6,'FontName', 'Arial');%, 'FontWeight','bold');
                caxis(coloraxis)
                colormap(fig.usecolormap)
                posx1=get(gca,'position');
                posx=get(clrbar,'Position');
                posx(1)= 1-fig.MarginRight+.01;
                posx(2)= fig.MarginBottom;
                posx(3)= .01;
                posx(4)= 1 - fig.MarginTop - fig.MarginBottom;
                set(clrbar,'Position',posx)
                set(gca,'position',posx1)
                sprtit = sprintf('%s %s %s %s D%s %d-%dHz_rCyc%.1f-%.1f', calcfunction,ixpc.datatypes{ian}{iDT}, strjoin(epochEnvironment,'-'), animalID, strrep(num2str(days,'%-2d'), ' ', '-'), min_freq,max_freq, range_cycles);
                iStitle = text(.5, .95, [{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]], ...
                    'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
                clrbartit = text(posx(1)+posx(3)/2, posx(2)-fig.MarginBottom/2, calcfunction, 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial','horizontalAlignment', 'center');
                
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    if ~isdir(figdirectory);
                        mkdir(figdirectory);
                    end
                    sprtitsave = strrep(sprtit,' ', '_'); %put '_' character in place of whitespace for filename
                    set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
                    currfigfile = sprintf('%s%s',figdirectory, sprtitsave);
                    print(currfigfile,'-dpng', '-r0')
                    disp(sprintf('plot %s saved', sprtit))
                end
                close all
            end
        end
    end
    if plotByArea
        warning('off', 'MATLAB:contour:ConstantData')
        fig = getFigspecs(areafigspecs);
        for ian = 1:length(ixpc.animals)
            [nt] = getAreainfo(ixpc.areas{ian}, ian, colorSet);
            animalID = ixpc.animals{ian};
            nAreas = length(ixpc.areas{ian});
            for iDT = 1:(size(ixpc.datatypesMat{ian},1)-2) %minus 2 because i haven't done the diff for out/in corr/err trials yet for ISPC
                conditions = cell2mat(ixpc.datatypesMat{ian}(iDT,:));
                [ifig, sfrows, sfcols] = prepFig(savefigs, pausefigs, nt, fig, 'figspecs', areafigspecs, 'nDays', length(conditions));
                %% ---- loop across all tets for this day ----
                for iar = 1:nAreas;
                    iID = nt.areas(iar,:);
                    icolor = nt.icolors(iar,:);
                    for icond = conditions;
                        isf = (iar-1)*length(conditions)+icond;
                        intfig = subaxis(sfrows,sfcols,isf, 'SpacingVertical', fig.SpacingVertical, 'SpacingHorizontal', fig.SpacingHorizontal, ...
                            'Padding', fig.Padding, 'MarginLeft', fig.MarginLeft, 'MarginRight', fig.MarginRight, 'MarginTop', fig.MarginTop, ...
                            'MarginBottom', fig.MarginBottom);
                        %% ---- Main Plotting ----
                        if plotlogbased
                            eval(sprintf('idata2plot = squeeze(ixpc.logbasedarea%smean{ian}{iar}{iDT}(:,icond,:))'';', calcfunction));
%                             eval(sprintf('idata2plot = squeeze(ixpc.logbasedarea%smeanUnspec{ian}{iar}{iday}(:,iDT,:))'';', plotoutputtype));
                        else
                            eval(sprintf('idata2plot = squeeze(ixpc.pctbasedarea%smean{ian}{iar}{iDT}(:,icond,:))'';', calcfunction));
%                             eval(sprintf('idata2plot = squeeze(ixpc.pctbasedarea%smeanUnspec{ian}{iar}{iday}(:,iDT,:))'';', plotoutputtype));
                        end
                        idata2plot = trim2win(idata2plot, ixpc.wp.srate, ixpc.wp.plotwin);
                        %                         try
                        %                             [~,bn] = contourf(intfig, wp.plottimeWin,wp.frex(1:end-1),idata2plot,fig.contourRes,'linecolor','none'); %
                        %                             set(gca,'ydir','normal','xlim',[wp.plotwin(1) wp.plotwin(2)], 'ylim',[wp.frex(1) wp.frex(end-1)])
                        %                         catch
                        [~,bn] = contourf(intfig, ixpc.wp.plottimeWin,ixpc.wp.frex,idata2plot,fig.contourRes,'linecolor','none'); %
                        set(gca,'ydir','normal','xlim',[ixpc.wp.plotwin(1) ixpc.wp.plotwin(2)], 'ylim',[ixpc.wp.min_freq ixpc.wp.max_freq])
                        %                         end
                        
                        eval(sprintf('caxis(fig.coloraxis%s);', calcfunction))
                        
                        %                         caxis(coloraxis)
                        colormap(fig.usecolormap)
                        hold on;
                        if (any(strcmp(ixpc.datatypes{ian}{1}{iDT}, 'outB-inB')) || any(strcmp(ixpc.datatypes{ian}{1}{iDT}, 'corrOut-mistOut'))) && plotzmask
                            %                         (iDT == 6 || iDT == 7) && plotzmask % plot the significance countours
                            eval(sprintf('zmask2plot = squeeze(ixpc.area%szmaskmean{ian}{iar}{iDT})'';', calcfunction));
                            %                         eval(sprintf('MCminmax = abs(ixpc.MC_%s_minmax{ian}{iar}{iDT})'';',plotoutputtype));
                            %                         eval(sprintf('irawdata = squeeze(ixpc.%sout{ian}{iar}(:,iDT,:))'';', plotoutputtype));
                            zmask2plot = trim2win(zmask2plot, ixpc.wp.srate, ixpc.wp.plotwin);
                            %                         irawdata = trim2win(irawdata, srate, plotwin);
                            zmask = zmask2plot;
                            zmask(abs(zmask2plot)<zval) = 0;
                            [~,h] = contour(intfig,ixpc.wp.plottimeWin,ixpc.wp.frex,logical(zmask),1);
                            h.LineColor = mcLineColor;
                            hold on;
                            %                         MCthresh = MCminmax(ceil(length(MCminmax)*(1-pval)));
                            %                         MCthreshmap = irawdata;
                            %                         MCthreshmap(abs(irawdata)<MCthresh) = 0;
                            %                         [~,mc] = contour(intfig,plottimeWin,frex(1:end-1),logical(MCthreshmap),1);
                            %                         mc.LineColor = 'm';%[.85 .85 .85]; %[.6 .6 .6];
                        end
                        %                     %                 set(gca,'ydir','normal','xlim',[plotwin(1) plotwin(2)], 'ylim',[frex(1) frex(end)])
                        %                     if mod(int, sfcols) ~= 1;
                        %                         set(gca, 'yscale','log','ytick',[])
                        %                         %                                         set(gca, 'yscale','log','ytick',ceil(logspace(log10(min(frex)),log10(max(frex)),6)))
                        %                         set(gca, 'FontSize',8, 'FontName', 'Arial');%
                        %                     else
                        set(gca, 'yscale','log','ytick',round(logspace(log10(min(ixpc.wp.frex)),log10(max(ixpc.wp.frex)),6)))
                        %                     set(gca, 'YTick',[round(linspace(min(frex),max(frex),6))])
                        set(gca, 'FontSize',fig.yticksize, 'FontName', fig.fontname);%
                        %                     end
                        ititleStr = sprintf('%s',strjoin(iID, '-'));
                        if isf > (sfrows-1)*sfcols; %if the last row
                            set(gca, 'XTick',[ixpc.wp.plotwin(1):ixpc.wp.plotwin(2)/2:ixpc.wp.plotwin(2)],'FontSize',fig.xticksize, 'FontName', fig.fontname);
                            set(gca, 'XTickLabelRotation',fig.xtickrotation);
                        else
                            set(gca, 'XTick', []);
                        end
                        if strcmp(areafigspecs, 'AreasByDays')
                            if isf <= sfcols; %if first row
                                idaytitle = title(sprintf('%d',icond));
                                if isf == 1;
                                    idaytitle = title(sprintf('Day %d',icond));
                                end
                                set(idaytitle,'FontSize',fig.daytitlesize,'Color', fig.daycolor, 'FontName', fig.fontname,'FontWeight','bold');
                            end
                            if ~mod(isf-1,sfcols); %if first column
                                if length(icolor) > 1
                                    ypos = ylabel(' ');
                                    possf=ypos.Position;
%                                     possftmp = bsxfun(@minus, possf([1 3]), .01);
%                                     possfY = possf;
%                                     possfY(1) = possftmp(1);
%                                     possfY(3) = possftmp(2);
%                                     ititle = ylabel(ititleStr);
                                    ititle = text(ypos.Position(1)*(2/3),1/ypos.Position(3), ['\color[rgb]',sprintf('{%d %d %d} %s-', icolor{1}, iID{1}), '\color[rgb]', sprintf('{%d %d %d} %s', icolor{2}, iID{2})],...
                                        'Parent', intfig, 'Units', 'normalized');
                                    set(ititle,'FontSize',fig.areatitlesize,'FontName', fig.fontname,'FontWeight','bold', 'VerticalAlignment', 'top',...
                                        'Rotation', 45, 'HorizontalAlignment', 'center');
                                else
                                    ititle = ylabel(ititleStr);
                                    set(ititle,'FontSize',fig.areatitlesize,'Color', icolor{1}, 'FontName', fig.fontname,'FontWeight','bold');
                                end
                            else
                                set(gca, 'yscale','log','ytick',[])
                            end
                        else
                            ititle = title(ititleStr);
                            set(ititle,'FontSize',fig.areatitlesize,'Color', icolor, 'FontName', fig.fontname,'FontWeight','bold');
                        end
                        %% ---- Source ripple line, subplot title ----
                        Xline = [0 0];
                        Yline = [ixpc.wp.min_freq ixpc.wp.max_freq];
                        line(Xline, Yline, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                    end
                end
                %% ---- super title and colorbar----
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                xlab = 'Time (s)';
%                 ylab = 'Frequency (Hz)';
                supxlabel = text(.5, .02, xlab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName',fig.fontname, 'Parent', sprtitleax,...
                    'Units', 'normalized', 'horizontalAlignment', 'center');
%                 supylabel = text(.01, .5, ylab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', fig.fontname, 'rotation', 90, ...
%                     'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
                
                clrbar = colorbar('location','eastoutside', 'FontSize',6,'FontName', fig.fontname);%, 'FontWeight','bold');
%                 caxis(fig.coloraxis)
                eval(sprintf('caxis(fig.coloraxis%s);', calcfunction))
                colormap(fig.usecolormap)
                posx1=get(gca,'position');
                posx=get(clrbar,'Position');
                posx(1)= 1-fig.MarginRight+.01;
                posx(2)= fig.MarginBottom;
                posx(3)= .01;
                posx(4)= 1 - fig.MarginTop - fig.MarginBottom;
                set(clrbar,'Position',posx)
                set(gca,'position',posx1)
                sprtit = sprintf('%s %s %s %s %s D%s %d-%dHz_rCyc%.1f-%.1f',areafigspecs, calcfunction,ixpc.datatypes{ian}{1}{iDT}, strjoin(epochEnvironment,'-'), animalID, strrep(num2str(days,'%-2d'), ' ', '-'), ixpc.wp.min_freq,ixpc.wp.max_freq, ixpc.wp.range_cycles);
                iStitle = text(.5, .95, [{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]], ...
                    'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', fig.fontname, 'horizontalAlignment', 'center');
                clrbartit = text(posx(1)+posx(3)/2, posx(2)-fig.MarginBottom/2, {calcfunction;'(dB-norm)'}, 'FontSize',10,'FontWeight','bold','Color','k', 'FontName',fig.fontname,'horizontalAlignment', 'center', 'verticalAlignment', 'middle');
                
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    if ~isdir(figdirectory);
                        mkdir(figdirectory);
                    end
                    sprtitsave = strrep(sprtit,' ', '_'); %put '_' character in place of whitespace for filename
                    set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
                    currfigfile = sprintf('%s%s',figdirectory, sprtitsave);
                    print(currfigfile,'-dpng', '-r0')
                    disp(sprintf('plot %s saved', sprtit))
                end
                close all
            end
        end
    end
end





    %% ---------------- plot ISPC---------------------------------------------------------------------------------------------
    %% ---------------- plot ISPC---------------------------------------------------------------------------------------------
    %% ---------------- plot ISPC---------------------------------------------------------------------------------------------
    if 0; %plotISPC
        clear F %save space on memory
        position = [.1 .1 .4 .3];
        SpacingHorizontal = 0.00;
        SpacingVertical = 0.00;
        Spacing = 0.00;
        Padding = 0.00;
        MarginLeft = 0.05;
        MarginRight = 0.05;
        MarginTop = 0.14;
        MarginBottom =  0.08;
        for ian = 1:length(ixpc.animals)
            %% ---- loadtetinfostruct ----
            animalinfo = animaldef(lower(animals{ian}));
            animalID = animalinfo{1,3}; %use anim prefix for name
            FFanimdir =  sprintf('%s',animalinfo{1,2});
            load([FFanimdir, animalID, 'tetinfo']);
            tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
            for icond = days
                ntrodePairIndices = ixpc.index{ian}{icond};
                tetpairs = unique(ntrodePairIndices(:,[3 6]),'rows','stable');
                ntetPairs = size(tetpairs,1);
                %% ---- reorder the LFP traces by suparea|area|subarea tags (in that priority) ----
                ntrodePairIndicesLIN = [ntrodePairIndices(:,[1:3]); ntrodePairIndices(:,[4:6])];
                ixpc.wp.nNTrodes = size(unique(ntrodePairIndicesLIN, 'stable'),1);
                [~, tagIndMap] = ismember(ntrodePairIndicesLIN,tetinfoAll.index, 'rows');
                ntrodeTags = tetinfoAll.values(tagIndMap);
                ntrodeTags = [ntrodeTags(1:length(ntrodeTags)/2) ntrodeTags(length(ntrodeTags)/2+1:end)];
                try
                    numsumSupAreas = cellfun(@(x) sum(uint16(x.suparea)), ntrodeTags, 'UniformOutput', false);
                    numsumAreas = cellfun(@(x) sum(uint16(x.area)), ntrodeTags, 'UniformOutput', false);
                    numsumSubAreas = cellfun(@(x) sum(uint16(x.subarea)), ntrodeTags, 'UniformOutput', false);
                    strSupAreas = cellfun(@(x) x.suparea, ntrodeTags, 'UniformOutput', false);
                    strAreas = cellfun(@(x) x.area, ntrodeTags, 'UniformOutput', false);
                    strSubAreas = cellfun(@(x) x.subarea, ntrodeTags, 'UniformOutput', false);
                catch
                    error('all ntrodes need to have a suparea, subarea, and area tag, even if blank')
                end
                icolors = [colorPicker(colorSet, strAreas(:,1), strSubAreas(:,1)) colorPicker(colorSet, strAreas(:,2), strSubAreas(:,2))];
                %             icolors = [icolors(1:length(icolors)/2,1:3) icolors(length(icolors)/2+1:end,1:3)];
                %             numsumtags = [numsumSupAreas numsumAreas numsumSubAreas];
                numsumallareatags = cell2mat([numsumSupAreas(:,1) numsumAreas(:,1) numsumSubAreas(:,1) numsumSupAreas(:,2) numsumAreas(:,2) numsumSubAreas(:,2)]);
                %             catnumsumallareatags = [numsumallareatags(1:length(numsumallareatags)/2, 1:3) numsumallareatags(length(numsumallareatags)/2+1:end, 1:3)];
                [numsumallSort, numsumSortInds] = sortrows(numsumallareatags);%,[-1 -2 -3]); % -Col to 'descend'
                icolors = icolors(numsumSortInds,:);
                %             sfrows = ntets;
                %             sfcols = ntets;
                %             onemat = logical(tril(ones(ntets,ntets),-1))'; %get indices of non-duplicates (below comb triangle)
                %             linsfpos = find(onemat); %get linear index of ntrode combination indices
                %% ---- loop across all tetpairs for this day ----
                for intrPair = 1:ntetPairs;
                    if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
                        ifig1 = figure('Visible','off','units','normalized','position',position);
                    else
                        ifig1 = figure('units','normalized','position',position);
                    end
                    set(gcf,'color','white')
                    intfig = subaxis(1,1,1,'SpacingVertical', SpacingVertical, 'SpacingHorizontal', SpacingHorizontal, ...
                        'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight, 'MarginTop', MarginTop, ...
                        'MarginBottom', MarginBottom);
                    introdeIDs = ntrodePairIndices(numsumSortInds(intrPair),[3 6]);
                    isupareatag1 = ntrodeTags{numsumSortInds(intrPair),1}.suparea;
                    iareatag1 = ntrodeTags{numsumSortInds(intrPair),1}.area;
                    isubareatag1 = ntrodeTags{numsumSortInds(intrPair),1}.subarea;
                    isupareatag2 = ntrodeTags{numsumSortInds(intrPair),2}.suparea;
                    iareatag2 = ntrodeTags{numsumSortInds(intrPair),2}.area;
                    isubareatag2 = ntrodeTags{numsumSortInds(intrPair),2}.subarea;
                    iNTcolor1 = icolors(intrPair,1:3);
                    iNTcolor2 = icolors(intrPair,4:6);
                    %                 intfig = subaxis(sfrows,sfcols,linsfpos(intrPair), 'SpacingVertical', SpacingVertical, 'SpacingHorizontal', SpacingHorizontal, ...
                    %                 'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight, 'MarginTop', MarginTop, ...
                    %                 'MarginBottom', MarginBottom);
                    %% ---- Main Plotting ----
                    contourf(timeWin,frex,squeeze(ixpc.phOut{ian}{icond}(numsumSortInds(intrPair),:,:))',numfrex,'linecolor','none');
                    %                 patch([-.8, .8, .8, -.8], [-10 -10 80 80], iNTcolor, 'edgecolor','none')
                    hold on;
                    set(gca,'clim',clims,'ydir','normal','xlim',[-.5 .5])
                    set(gca, 'YScale', 'log')
                    
                    %                 if mod(intrPair, sfcols) ~= 1;
                    %                     set(gca, 'YTick', []);
                    %                 else
                    set(gca, 'YTick',[round(linspace(min(frex),max(frex),6))],'FontSize',8, 'FontName', 'Arial');%
                    %                 end
                    %                 if intrPair <= (sfrows-1)*sfcols;
                    %                     set(gca, 'XTick', []);
                    %                 else
                    set(gca, 'XTick',[win(1):win(2)/2:win(2)],'FontSize',8, 'FontName', 'Arial');%
                    %                 end
                    %% ---- Source ripple line, subplot title ----
                    Xline = [0 0];
                    Yline = [min_freq max_freq];
                    line(Xline, Yline, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                    
                    %                 set([{ititle1}; {ititle2}],'FontSize',10,'Color', iNTcolor, 'FontName', 'Arial','FontWeight','bold'); %
                    
                    %% ---- super title and colorbar----
                    sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig1);
                    xlab = 'Time (s)';
                    ylab = 'Frequency (Hz)';
                    supxlabel = text(.5, .02, xlab, 'FontSize',8,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'Parent', sprtitleax,...
                        'Units', 'normalized', 'horizontalAlignment', 'center');
                    supylabel = text(.01, .5, ylab, 'FontSize',8,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'rotation', 90, ...
                        'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
                    caxis(clims);
                    colormap(usecolormap)
                    clrbar = colorbar('location','eastoutside', 'FontSize',6,'FontName', 'Arial');%, 'FontWeight','bold');
                    posx1=get(gca,'position');
                    posx=get(clrbar,'Position');
                    posx(1)= 1-MarginRight+.01;
                    posx(2)= MarginBottom;
                    posx(3)= .01;
                    posx(4)= 1 - MarginTop - MarginBottom;
                    set(clrbar,'Position',posx)
                    set(gca,'position',posx1)
                    %             sprtit = sprintf('%s D%s E%s T%d', ianimalinfo{1}, strjoin(arrayfun(@(x) num2str(x),days','UniformOutput',false),'-'), strjoin(arrayfun(@(x) num2str(x),idayEpTet(:,2)','UniformOutput',false),'-'), idayEpTet(1,3));
                    sprtit = sprintf('%s %s %s D%d', calcfunction, epochEnvironment, animalID, icond);
                    %             if plotNTrodesAcrossDays
                    %                 sprTags = sprintf('%s %s', iarea, isubarea);
                    %                 iclr = icolors(iIndInds(1),:);
                    ititle1 = sprintf('%d %s %s',introdeIDs(1), iareatag1, num2str(isubareatag1));
                    ititle2 = sprintf('%d %s %s',introdeIDs(2), iareatag2, num2str(isubareatag2));
                    iStitle = text(.5, .93, [{filenameTitle};...
                        {['\fontsize{10} \color[rgb]' sprintf('{%d %d %d} %s ', iNTcolor1, ititle1)]};...
                        {['\fontsize{10} \color[rgb]' sprintf('{%d %d %d} %s ', iNTcolor2, ititle2)]}]);
                    %                     {['\fontsize{8} \color[rgb]{.5 .5 .5}', sprintf('{%s}', resultfilename)]}], 'Parent', sprtitleax, 'Units', 'normalized');
                    
                    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
                    %             else
                    %                 iStitle = text(.5, .95, [{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]], ...
                    %                     'Parent', sprtitleax, 'Units', 'normalized');
                    %             end
                    set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
                    clrbartit = text(posx(1)+posx(3)/2, posx(2)-MarginBottom/2, calcfunction, 'FontSize',8,'FontWeight','bold','Color','k', 'FontName', 'Arial','horizontalAlignment', 'center');
                    %% ---- pause, save figs ----
                    %             return
                    if pausefigs
                        pause
                    end
                    if savefigs
                        if ~isdir(figdirectory);
                            mkdir(figdirectory);
                        end
                        if ~isdir([figdirectory resultfilename]);
                            mkdir([figdirectory resultfilename]);
                        end
                        sprtitsave = strrep(sprtit,' ', '_'); %put '_' character in place of whitespace for filename
                        set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
                        currfigfile = sprintf('%s%s/%s_nT%d-%d',figdirectory, resultfilename, sprtitsave, introdeIDs);
                        print(currfigfile,'-dpng', '-r0')
                        disp(sprintf('plot %s nT%d-%d saved (%d of %d)', sprtit, introdeIDs, intrPair, ntetPairs))
                    end
                    close all
                end
            end
        end
    end
    
