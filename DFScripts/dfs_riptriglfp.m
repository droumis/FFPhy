
% clear all
close all
runFilterFramework = 0;
saveFilterOutput = 0;% runFilterFramework;
loadFilterOutput = 0;
% EpochMean = 1;
resaveFilterOutput = 0;
plotLFPtraces = 1;
savefigs = 0;
pausefigs = 1;
% plotSpec_EpochMean = 1;
% plotSpec_allEpochs = 0;
outputDirectory = '/typhoon/droumis/analysis';
%% ---------------- plotting params --------------------------
AreaLabelRotate = 45;
usecolormap = 'colorcube';
win = [.5 .5]; %in seconds
indwin = win*1500;
%% ---------------- Data Filters --------------------------

animals = {'JZ1'};
% animals = {'JZ1', 'D13'};
days = [1:14];
filtfunction = 'riptriglfp';
% LFPtypes = [{'eeg'}, {'ripple'}];
eventtype = 'rippleskons';
% eventarea = 'ca1';
epochEnvironment = 'wtrack';% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
epochType = 'run';
eventSourceArea = 'ca1';
tetAreas = {'ca1', 'mec', 'por'}; %ca1, mec, por

consensus_numtets = 1;   % minimum # of tets for consensus event detection
minthresh = 3;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
minvelocity = 0;
maxvelocity = 4;
filename = sprintf('%s_%s%s_%s_%s.mat', filtfunction, eventSourceArea, eventtype, epochEnvironment, cell2mat(animals));
filenameTitle = strrep(filename,'_', '\_');
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1;
    epochfilter =    sprintf('(isequal($type, ''%s'')) && (isequal($environment, ''%s''))',epochType, epochEnvironment); %'isequal($type, ''run'') && (isequal($environment, ''MultipleW''))'; %%'(isequal($type, ''sleep''))'; %%%&& isequal($descript, ''post-allruns''))';%   %%% %'isequal($type, ''run'') && isequal($environment, ''WTrackA'') && (($exposure>=1) && ($exposure<=10))';  %
    iterator = 'multitetrodeanal'; %multitetrodeanal
    %tetfilter: the ntrodeID's that pass this filter get stored into f.eegdata{day}{epoch}[ntrodeIDs]
    tetfilter = sprintf('(isequal($area,''%s'') || isequal($area,''%s'') || isequal($area,''%s''))', tetAreas{1}, tetAreas{2}, tetAreas{3}); % || isequal($area,''v2l'') || isequal($area,''sub''))';
    % timefilter{1} = {'get2dstate','($velocity<4)'};
    % timefilter{2} = {'kk_getriptimes','($nripples>=1)',[],'tetfilter',tetfilter,'minthresh',5};
    timefilter{1} = {'getconstimes', '($cons == 1)', [eventSourceArea eventtype],1,'consensus_numtets',consensus_numtets,...
        'minthresh',minthresh,'exclusion_dur',exclusion_dur,'minvelocity',minvelocity,'maxvelocity',maxvelocity};
    %----------F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);--------
    F = createfilter('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);
    %----------f = setfilteriterator(f, funcname, loadvariables, options)--------
    F = setfilterfunction(F, ['dfa_' filtfunction], {'eeg', 'ripple', [eventSourceArea eventtype]},'eventtype',[eventSourceArea eventtype]);
    tic
    F = runfilter(F);
    F(1).filterTimer = toc; F(1).filterTimer
    F(1).worldDateTime = clock;
    F(1).dataFilter = struct('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator, 'filename', filename);
end
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    if ~isdir(sprintf('%s/filter_output/%s/', outputDirectory, filtfunction));
        mkdir(sprintf('%s/filter_output/%s/', outputDirectory, filtfunction));
    end
    %save the entire workspace for filter provenance
    save(sprintf('%s/filter_output/%s/%s',outputDirectory, filtfunction, filename), 'F','-v7.3');
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s/filter_output/%s/%s',outputDirectory,filtfunction, filename))
end

if plotLFPtraces %to do
    
    %% STEP 2: For each EEG + Ripple Type... select nTrodes and Lookup Rip time windows
    
    for iLFPtype = 1%:length(LFPtypes); % For each LFP type (wideband EEG, ripple band, etc), load all of the regions LFP files into eegstruct
        
        %% plot LFP traces for all areas for each ripple time window
        if plotLFPtraces
            for ianimal = 1:length(F)
                %loadtetinfostruct
                animalinfo = animaldef(lower(animals{ianimal}));
                animalID = animalinfo{1,3}; %use anim prefix for name
                FFanimdir =  sprintf('%s',animalinfo{1,2});
                load([FFanimdir, animalID, 'tetinfo']);
                tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
                for iday = 1:length(F(ianimal).output)
                    %load the ripplekons file for each area for this day
                    for iareakons = 1:length(tetAreas);
                        load(sprintf('%s%s%s%s%02d.mat', FFanimdir, animalID, tetAreas{iareakons}, eventtype, iday));
                        eval([sprintf('eventKons.%s%s = %s%s;', tetAreas{iareakons}, eventtype, tetAreas{iareakons}, eventtype)])
                        eval([sprintf('clear %s%s', tetAreas{iareakons}, eventtype)])
                    end
                    for iepoch = 1:length(F(ianimal).output{iday})
                        if isempty(F(ianimal).output{iday}(iepoch).data)
                            continue %if this an day epoch is empty, skip to the next
                        end
                        %get data info for rips within this day epoch
                        iepochLFPtimes = []; ripStartIndices = []; ripEndIndices = []; ntrodesIndices = []; win = [];
                        iepochLFPtimes = F(ianimal).output{iday}(iepoch).LFPtimes;
                        ripStartIndices = F(ianimal).output{iday}(iepoch).eventStartIndices;
                        ripEndIndices = F(ianimal).output{iday}(iepoch).eventEndIndices;
                        win = F(ianimal).output{iday}(iepoch).win;
                        indwin = win*1500;
                        ntrodesIndices = F(ianimal).output{iday}(iepoch).index;
                        %reorder the LFP traces by suparea|area|subarea tags (in that priority)
                        [~, tagIndMap] = ismember(ntrodesIndices,tetinfoAll.index, 'rows');
                        ntrodeTags = tetinfoAll.values(tagIndMap);
                        try
                            numsumSupAreas = cellfun(@(x) sum(uint8(x.suparea)), ntrodeTags, 'UniformOutput', false);
                            numsumAreas = cellfun(@(x) sum(uint8(x.area)), ntrodeTags, 'UniformOutput', false);
                            numsumSubAreas = cellfun(@(x) sum(uint8(x.subarea)), ntrodeTags, 'UniformOutput', false);
                            strSupAreas = cellfun(@(x) x.suparea, ntrodeTags, 'UniformOutput', false);
                            strAreas = cellfun(@(x) x.area, ntrodeTags, 'UniformOutput', false);
                            strSubAreas = cellfun(@(x) x.subarea, ntrodeTags, 'UniformOutput', false);
                        catch
                            error('all ntrodes need to have a suparea, subarea, and area tag, even if blank')
                        end
                        numsumallareatags = cell2mat([numsumSupAreas numsumAreas numsumSubAreas]);
                        [numsumallSort, numsumSortInds] = sortrows(numsumallareatags);%,[-1 -2 -3]); % -Col to 'descend'
                        for irip = 1:length(F(ianimal).output{iday}(iepoch).data{1})
                            if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
                                ifig = figure('Visible','off','units','normalized','position',[.1 .1 .6 .8]);
                            else
                                ifig = figure('units','normalized','position',[.1 .1 .6 .8]);
                            end
                            clrmat = (colormap(usecolormap)); %flipud
                            colorslength = length(clrmat(:,1));
                            iripWinTimes = iepochLFPtimes(ripStartIndices(irip)-indwin: ripStartIndices((irip))+indwin);
                            %                             icolor = mod(allareatagsSorted(:,:), colorslength);
                            %                             icolor = 1-(icolor-min(min(icolor)))./max(max(icolor)); %normalize
                            icolors = clrmat(mod(prod(numsumallSort,2), colorslength), :);
                            [uniqColors, uniqColorsInds, uniqColorsInds2] = unique(icolors, 'rows', 'stable');
                            [uniqNumsum, uniqNumSumInds, uniqNumSumInds2] = unique(numsumallSort, 'rows', 'stable');
                            titleColors = icolors(uniqColorsInds,:);
                            titleSupAreas = strSupAreas(numsumSortInds(uniqColorsInds));
                            titleAreas = strAreas(numsumSortInds(uniqColorsInds));
                            titleSubAreas = strSubAreas(numsumSortInds(uniqColorsInds));
                            iYlim = [];
                            for iLFPtype = 1:length(F(ianimal).output{iday}(iepoch).LFPtypes)
                                for introde = 1:length(F(ianimal).output{iday}(iepoch).data{iLFPtype}{irip}(:,1))
                                    introdeID = F(ianimal).output{iday}(iepoch).index(numsumSortInds(introde), 3);
                                    isupareatag = ntrodeTags{numsumSortInds(introde)}.suparea;
                                    iareatag = ntrodeTags{numsumSortInds(introde)}.area;
                                    isubareatag = ntrodeTags{numsumSortInds(introde)}.subarea;
                                    %                                     iTags = {{isupareatag},{iareatag}, {isubareatag}};
                                    iNTcolor = icolors(introde,:);
                                    figure(ifig)
                                    subaxis(1,length(F(ianimal).output{iday}(iepoch).data),iLFPtype, 'Spacing', 0.02, 'Padding', 0.0, 'Margin', 0.09);
                                    %                                     subplot(1,length(F(ianimal).output{iday}(iepoch).data),iLFPtype)
                                    introdeiripLFP = F(ianimal).output{iday}(iepoch).data{iLFPtype}{irip}(numsumSortInds(introde),:);
                                    
                                    
                                    if introde == 1;
                                        plot(iripWinTimes, introdeiripLFP, 'Color', iNTcolor, 'Linewidth', 1)
                                        hold on;
                                        traceVertOffset = 0;
                                        previousLFP = introdeiripLFP;
                                        iYlim(1,1) = max(introdeiripLFP); %first trace max
                                    else
                                        traceVertOffset(introde) = traceVertOffset(introde-1) + abs(min(previousLFP)) + abs(max(introdeiripLFP));
                                        plot(iripWinTimes, introdeiripLFP - traceVertOffset(introde), 'Color', iNTcolor, 'LineWidth',1);
                                        hold on;
                                        previousLFP = introdeiripLFP;
                                    end
                                    
                                    %get index of matching region iareakons
                                    %set patch for this ntrode
                                    eval([sprintf(' ripsStartTimes = eventKons.%s%s{%d}{%d}{1}.starttime;', iareatag, eventtype, iday, iepoch)])
                                    eval([sprintf(' ripsEndTimes = eventKons.%s%s{%d}{%d}{1}.endtime;', iareatag, eventtype, iday, iepoch)])
                                    %                                     ripsStartTimes = eventKons{1}.ca1rippleskons{iday}{iepoch}{1}.starttime;
                                    %                                     ripsEndTimes = eventKons{1}.ca1rippleskons{iday}{iepoch}{1}.endtime;
                                    ripsinwinInds = find((ripsStartTimes>iripWinTimes(1) & ripsStartTimes<iripWinTimes(end))); %if the rip start time is within the current window
                                    ripsinwinTimes = [ripsStartTimes(ripsinwinInds) ripsEndTimes(ripsinwinInds)];
                                    for iripinwin = 1:length(ripsinwinTimes(:,1))
                                        ilmax = max(introdeiripLFP - traceVertOffset(introde));
                                        ilmin = min(introdeiripLFP - traceVertOffset(introde));
                                        iripstarttime = ripsStartTimes(ripsinwinInds(iripinwin));
                                        iripendtime = ripsEndTimes(ripsinwinInds(iripinwin));
                                        %                                     line([iripstarttime iripstarttime], [ilmin ilmax],'Color',[.8 .8 .8],'LineWidth',1.5);
                                        Xpatch = [iripstarttime iripendtime iripendtime iripstarttime];
                                        Ypatch = [ilmin ilmin ilmax ilmax];
                                        patch(Xpatch, Ypatch, iNTcolor, 'FaceAlpha', .25, 'edgecolor','none'); %triggering-ripple patch
                                    end
                                    NTlabel = text(iripWinTimes(1)-.01, -traceVertOffset(introde), num2str(introdeID));
                                    set(NTlabel, 'FontName', 'Arial', 'FontSize',8, 'FontWeight', 'normal', 'horizontalAlignment', 'right');
                                    %overlay event window patch for this
                                    %ntrode
                                    ilmax = max(introdeiripLFP - traceVertOffset(introde));
                                    ilmin = min(introdeiripLFP - traceVertOffset(introde));
                                    iripstarttime = iepochLFPtimes(ripStartIndices(irip));
                                    iripendtime = iepochLFPtimes(ripEndIndices(irip));
                                    Xpatch = [iripstarttime iripstarttime];
                                    Ypatch = [ilmin ilmax];
                                    %                                     patch(Xpatch, Ypatch, 'k', 'edgecolor','none', 'FaceAlpha', .1);
                                    line(Xpatch, Ypatch, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1.5);
                                    %                                     patch(Xpatch, Ypatch, Xpatch, 'FaceAlpha', .1, 'edgecolor','none',); %triggering-ripple patch
                                end
                                ititle = title(sprintf('%s',F(ianimal).output{iday}(iepoch).LFPtypes{iLFPtype}));
                                set(ititle,'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial');
                                iYlim(1,2) = min(introdeiripLFP - traceVertOffset(introde)); %last trace min
                                set(gca,'children',flipud(get(gca,'children'))) %send the patch behind the LFP traces
                                ylim([iYlim(1,2) iYlim(1,1)])
                                xlim([iripWinTimes(1) iripWinTimes(end)])
                                set(gca, 'YTick', []);
                                set(gca,'XTick',[iripWinTimes(1):(iripWinTimes(end) - iripWinTimes(1))/10:iripWinTimes(end)], 'FontSize',8,'FontWeight','normal')
                                set(gca, 'XTickLabel', round([iripWinTimes(1):(iripWinTimes(end) - iripWinTimes(1))/10:iripWinTimes(end)],2), 'XTickLabelRotation',30);
                                %                             set(gca, 'XTickLabel', [-win(1):win(1)/5:win(2)])
                                xlabel('time(s)','FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial')
                                if iLFPtype == 1
                                    for ititleTag = 1:length(titleColors(:,1))
                                        try %place in the middle of this area's traces
                                            titleY = -traceVertOffset(uniqColorsInds(ititleTag)) - (traceVertOffset(uniqColorsInds(ititleTag+1)) - traceVertOffset(uniqColorsInds(ititleTag)))/2;
                                        catch %will catch when on the last title, at this point, place between the current vert offset and total minimum
                                            titleY = -traceVertOffset(uniqColorsInds(ititleTag)) - (-iYlim(1,2) - traceVertOffset(uniqColorsInds(ititleTag)))/2;
                                        end
                                        itxt = text(iripWinTimes(1)-.05, titleY, ['{\color[rgb]' sprintf('{%d %d %d} %s %s}',titleColors(ititleTag,:), titleAreas{ititleTag}, num2str(titleSubAreas{ititleTag}))]);
                                        %                                         set(itxt, 'rotation', AreaLabelRotate, 'FontName', 'Arial', 'FontSize',12, 'FontWeight', 'Bold');
                                        set(itxt, 'rotation', 0, 'FontName', 'Arial', 'FontSize',12, 'FontWeight', 'Bold', 'horizontalAlignment', 'right');
                                    end
                                end
                                %                                 if iLFPtype == length(F(ianimal).output{iday}(iepoch).LFPtypes)
                                %                                     iStitle = supertitle([{sprintf('%s D%d E%d R%d', animalID, iday, iepoch, irip)}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]]);
                                %                                     set(iStitle, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial');
                                %                                 end
                            end
                            iStitle = supertitle([{sprintf('%s D%d E%d R%d', animalID, iday, iepoch, irip)}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]]);
                            set(iStitle, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial');
                            %                             return
                            if pausefigs
                                pause
                                close all
                            end
                        end
                    end
                end
            end
        end
    end
end

if 0
    for currrip=1:length(ripsStartTimes{srcRegionind}) %for each ripple from source region..
        close all
        clear ripsinwin
        if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
            fig = figure('Visible','off','units','normalized','position',[.1 .1 .2 .8]);
        else
            fig = figure('units','normalized','position',[.1 .1 .2 .8]);
        end
        for p = 1:length(YripLFPdataMAT{iLFPtype}{currrip}(1,:)); %for each lfp trace
            if p == 1; %if the first trace
                plot(Xwindowtimes{iLFPtype}{currrip},YripLFPdataMAT{iLFPtype}{currrip}(:,p), 'Color',regionclr(lfptraceLUTregion(p),:), 'LineWidth',1); hold on;
                traceVertOffset = 0;
            else %compute offset from the absMax+absMin to place the current trace under the last trace
                traceVertOffset(p) = traceVertOffset(p-1) + abs(min(YripLFPdataMAT{iLFPtype}{currrip}(:,p-1))) + abs(max(YripLFPdataMAT{iLFPtype}{currrip}(:,p))); %use e.g. abs(max(YripLFPdataMAT{LFP}{currrip}(:,p)))/3 to overlay traces more
                plot(Xwindowtimes{iLFPtype}{currrip},YripLFPdataMAT{iLFPtype}{currrip}(:,p) - traceVertOffset(p), 'Color',regionclr(lfptraceLUTregion(p),:), 'LineWidth',1); hold on;
            end
        end
        %% Plot Ripple times as patchs and other accessory ploting stuff
        for iArea = 1:length(regions); %for each region
            ripsinwinInds = (ripsStartTimes{iArea}>WindowStartEndTimes(currrip,1) & ripsStartTimes{iArea}<WindowStartEndTimes(currrip,2)); %if the start time is within the current window
            ripsinwin{currrip}{iArea} = [ripsStartTimes{iArea}(ripsinwinInds) ripsEndTimes{iArea}(ripsinwinInds)];
            if ~isempty(ripsinwin{currrip}{iArea}) %if there are any ripples from this region in this window
                for m = 1:length(ripsinwin{currrip}{iArea}(:,1)) %for each ripple within the current window
                    Ylfpranges4region{iArea} = [-traceVertOffset(find(lfptraceLUTregion == iArea,1,'last'))-(abs(min(YripLFPdataMAT{iLFPtype}{currrip}(:,find(lfptraceLUTregion == iArea,1,'last'))))) -traceVertOffset(find(lfptraceLUTregion == iArea,1,'first'))+max(YripLFPdataMAT{iLFPtype}{currrip}(:,find(lfptraceLUTregion == iArea,1,'first')))];
                    %                 Ylfpranges4region{i} = [-lfpoffset(find(lfptraceLUTregion == i,1,'last'))-(abs(min(YripLFPdataMAT{LFP}{currrip}(find(lfptraceLUTregion == i,1,'last'))))) abs(max(YripLFPdataMAT{LFP}{currrip}(:,find(lfptraceLUTregion == i,1,'first'))))];
                    line([ripsinwin{currrip}{iArea}(m,1) ripsinwin{currrip}{iArea}(m,1)], Ylfpranges4region{iArea},'Color',regionclr(iArea,:),'LineWidth',1.5);
                    Xpatch = [ripsinwin{currrip}{iArea}(m,1) ripsinwin{currrip}{iArea}(m,2) ripsinwin{currrip}{iArea}(m,2) ripsinwin{currrip}{iArea}(m,1)];
                    %                 Ypatch = [Ylfpranges4region{i}(1) Ylfpranges4region{i}(1)+diff(Ylfpranges4region{i})/4 Ylfpranges4region{i}(2)-diff(Ylfpranges4region{i})/4 Ylfpranges4region{i}(2)]; %trapezoidal patch that decays toward rip end
                    Ypatch = [Ylfpranges4region{iArea}(1) Ylfpranges4region{iArea}(1) Ylfpranges4region{iArea}(2) Ylfpranges4region{iArea}(2)];
                    patch(Xpatch, Ypatch, patchclr(iArea,:), 'edgecolor','none'); %triggering-ripple patch
                end
            end
        end
        centerripStartTime=ripsStartTimes{srcRegionind}(currrip);
        Yextent = [-traceVertOffset(end)-(abs(min(YripLFPdataMAT{iLFPtype}{currrip}(:,end)))) max(max(YripLFPdataMAT{iLFPtype}{currrip}))];
        line([centerripStartTime centerripStartTime], Yextent ,'Color',regionclr(srcRegionind,:), 'LineStyle', '--','LineWidth',1.5) %line for the center trigger-ripple
        set(gca,'children',flipud(get(gca,'children'))) %send the patch behind the LFP traces
        ylim([Yextent(1) Yextent(2)])
        xlim([WindowStartEndTimes(currrip,1) WindowStartEndTimes(currrip,2)])
        xl = xlim;
        for k=1:length(regions)
            text(xl(1)-diff(WindowStartEndTimes(currrip,:))/10, -traceVertOffset(find(lfptraceLUTregion == k,1,'first')), regions{k}, 'Color', regionclr(k,:),'FontSize',13,'FontWeight','bold')
        end
        set(gca, 'YTick', []);
        set(gca,'XTick',[WindowStartEndTimes(currrip,1):diff(WindowStartEndTimes(currrip,:))/10:WindowStartEndTimes(currrip,2)], 'FontSize',10,'FontWeight','bold')
        set(gca, 'XTickLabel', [-windowsize:windowsize/5:windowsize])
        xlabel('seconds from rip start','FontSize',12,'FontWeight','bold','Color','k')
        title({[sprintf('%s d%de%d %sRip(%d) Triggered %s-LFP',animal, day, epoch, rippleregion, currrip, LFPtypes{iLFPtype})];[sprintf('timestamp: %16.f', centerripStartTime)]; [sprintf('peak nstds: %d', ripout{srcRegionind}{day}{epoch}.maxthresh(currrip))]},'FontSize',12,'FontWeight','bold')
        
        %% Pause and/or Save Figs
        if pausefigs
            pause
        end
        if savefigs
            if iLFPtype == 1 && currrip == 1; %make a new directory with datetime string if it's the first rip of the first LFPtype
                currfigdirectory = sprintf('%s/%s_d%02d-e%02d_%s-ripTriggeredLFP/',figdirectory, dateprintstr, day, epoch, regions{srcRegionind});
                mkdir(currfigdirectory)
            end
            set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
            currfigfile = sprintf('%s%s_%s%d',currfigdirectory, regions{srcRegionind}, LFPtypes{iLFPtype}, currrip);
            print(currfigfile,'-dpng', '-r0')
            disp(sprintf('%s plot %d of %d saved', LFPtypes{iLFPtype}, currrip, length(ripsStartTimes{srcRegionind})))
        end
        %         end
        %     end
        % end
        
        %% Combine Ripple + EEG plots in same directory and create merged plots to facilitate viewing
        
        if  createmergedplots && plotLFPtraces
            %copy figdirectory
            %     combdir = sprintf('/mnt/data19/droumis/D10/Figs/eegripple_combined/%s/',dateprintstr);
            %     mkdir(combdir)
            %     cd(combdir)
            cd(currfigdirectory);
            for iArea = 1:length(ripsStartTimes{srcRegionind})
                system(sprintf('convert +append *ripple%d.png *eeg%d.png ca1_eegripple%d.png',iArea,iArea,iArea)); %horizontally append eeg and ripple plots and save
            end
            disp('done combining images')
        end
        %% Do ripple band plotting, then population LFP fig, then population ripple timing fig tonight
        
        %plot each other regions' population source-rip triggered riptimes histogrammed
        
        %plt each other regions' LFP and rip band LFP envelope source-rip triggered average trace.. also try to square then squareroot version of this to account for polarity
        
        %take all the ripples in a given std and compute (1) the likelihood of rips in both MEC/crtx within some time window
        %(2) the correlation with peak std in rip band within mec, cortex... is the scatterplot of ca1 vs cortex it bimodal, while ca1-mec unimodal?? aka is mec a ripple gate
        %(3) the correlation with peak std in widband within mec, cortex
        %(4) all these measures seperated by whether the animal is on incorrect or correct trial .. use ripsinstate
        
        
        %% Tomorrow/this weekend do propagation for every stf category of ripple
        
        %save each other regions ripple time distance from source rip start to plot/save population histograms later
        
        %save a cell array with a cell for every sourcerip with all the snippets and other region rip times
    end
end