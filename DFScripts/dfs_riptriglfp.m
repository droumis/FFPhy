
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
days = [1:14];
filtfunction = 'riptriglfp';
LFPtype1 = 'eeg';
LFPtype2 = 'ripple';
eventtype = 'ca1rippleskons';
% eventarea = 'ca1';
epochEnvironment = 'sleep';% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
epochType = 'sleep';
% tetAreas = ['ca1', 'mec', 'por']; %ca1, mec, por

consensus_numtets = 1;   % minimum # of tets for consensus event detection
minthresh = 3;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
minvelocity = 0;
maxvelocity = 4;
% filenameTitle = sprintf('%s\\_%s\\_%s\\_%s', filtfunction, eventtype, epochEnvironment, cell2mat(animals));
filename = sprintf('%s_%s_%s_%s.mat', filtfunction, eventtype, epochEnvironment, cell2mat(animals));
stitlename = strrep(filename,'_', ' ');
stitlename  = strrep(stitlename,'.mat','');
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1;
    epochfilter =    sprintf('(isequal($type, ''%s'')) && (isequal($environment, ''%s''))',epochType, epochEnvironment); %'isequal($type, ''run'') && (isequal($environment, ''MultipleW''))'; %%'(isequal($type, ''sleep''))'; %%%&& isequal($descript, ''post-allruns''))';%   %%% %'isequal($type, ''run'') && isequal($environment, ''WTrackA'') && (($exposure>=1) && ($exposure<=10))';  %
    iterator = 'multitetrodeanal'; %multitetrodeanal
    %tetfilter: the ntrodeID's that pass this filter get stored into f.eegdata{day}{epoch}[ntrodeIDs]
    tetfilter = '(isequal($area,''ca1'') || isequal($area,''mec'') || isequal($area,''por''))'; % || isequal($area,''v2l'') || isequal($area,''sub''))';
    % timefilter{1} = {'get2dstate','($velocity<4)'};
    % timefilter{2} = {'kk_getriptimes','($nripples>=1)',[],'tetfilter',tetfilter,'minthresh',5};
    timefilter{1} = {'getconstimes', '($cons == 1)', eventtype,1,'consensus_numtets',consensus_numtets,...
        'minthresh',minthresh,'exclusion_dur',exclusion_dur,'minvelocity',minvelocity,'maxvelocity',maxvelocity};
    %----------F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);--------
    F = createfilter('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);
    %----------f = setfilteriterator(f, funcname, loadvariables, options)--------
    F = setfilterfunction(F, ['dfa_' filtfunction], {LFPtype1, LFPtype2, eventtype},'eventtype',eventtype);
    tic
    F = runfilter(F);
    F.filterTimer = toc; F.filterTimer
    F.worldDateTime = clock;
    F.dataFilter = struct('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator, 'filename', filename);
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

            for ianimal = 1:length(animals)
                %loadtetinfostruct
                ianimal = 1;
                animalinfo = animaldef(lower(animals{ianimal}));
                animalID = animalinfo{1,3}; %use anim prefix for name
                FFanimdir =  sprintf('%s',animalinfo{1,2});
                load([FFanimdir, animalID, 'tetinfo']);
                tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
                for iday = 1%:length(F.output)
                    for iepoch = 1%:length(F.output{iday})
                        %get data info for rips within this day epoch
                        iepochLFPtimes = []; ripStartIndices = []; ripEndIndices = []; ntrodesIndices = []; win = [];
                        iepochLFPtimes = F.output{iday}(iepoch).LFPtimes;
                        ripStartIndices = F.output{iday}(iepoch).eventStartIndices;
                        ripEndIndices = F.output{iday}(iepoch).eventEndIndices;
                        win = F.output{iday}(iepoch).win;
                        indwin = win*1500;
                        ntrodesIndices = F.output{iday}(iepoch).index;
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
                        for irip = 1%:length(F.output{iday}(iepoch).eegdata)
                            if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
                                ifig = figure('Visible','off','units','normalized','position',[.1 .1 .6 .8]);
                            else
                                ifig = figure('units','normalized','position',[.1 .1 .6 .8]);
                            end
                            clrmat = (colormap(usecolormap));
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
                            for iLFPtype = 1:length(LFPtypes)
                                for introde = 1:length(F.output{iday}(iepoch).eegdata{irip}(:,1))
                                    isupareatag = ntrodeTags{numsumSortInds(introde)}.suparea;
                                    iareatag = ntrodeTags{numsumSortInds(introde)}.area;
                                    isubareatag = ntrodeTags{numsumSortInds(introde)}.subarea;
                                    iTags = {{isupareatag},{iareatag}, {isubareatag}};
                                    %                                     iTags = {{strSupAreas(areatagsSortedInds(introde))},{strSupAreas(areatagsSortedInds(introde))}, {strSupAreas(areatagsSortedInds(introde))}};
                                    iNTcolor = icolors(introde,:);
                                    %                                     ripNTcolors(introde) = iNTcolor;
                                    figure(ifig)
                                    subplot(1,2,1)
                                    introdeiripLFP = F.output{iday}(iepoch).eegdata{irip}(numsumSortInds(introde),:);
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
                                    ititle = title('wideband(1-400Hz)');
                                    set(ititle,'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial');
                                    %overlay event window patches
                                    ilmax = max(introdeiripLFP - traceVertOffset(introde));
                                    ilmin = min(introdeiripLFP - traceVertOffset(introde));
                                    iripstarttime = iepochLFPtimes(ripStartIndices(irip));
                                    iripendtime = iepochLFPtimes(ripEndIndices(irip));
                                    line([iripstarttime iripstarttime], [ilmin ilmax],'Color',[.8 .8 .8],'LineWidth',1.5);
                                    %                                 line([traceVertOffset(introde) traceVertOffset(introde)], [],'Color',iNTcolor,'LineWidth',1.5);
                                    Xpatch = [iripstarttime iripendtime iripendtime iripstarttime];
                                    Ypatch = [ilmin ilmin ilmax ilmax];
                                    patch(Xpatch, Ypatch, [.8 .8 .8], 'edgecolor','none'); %triggering-ripple patch
                                end
                            end
                            iYlim(1,2) = min(introdeiripLFP - traceVertOffset(introde)); %last trace min
                            for ititleTag = 1:length(titleColors(:,1))
                                try %place in the middle of this area's traces
                                    titleY = -traceVertOffset(uniqColorsInds(ititleTag)) - (traceVertOffset(uniqColorsInds(ititleTag+1)) - traceVertOffset(uniqColorsInds(ititleTag)))/2;
                                catch %will catch when on the last title, at this point, place between the current vert offset and total minimum
                                    titleY = -traceVertOffset(uniqColorsInds(ititleTag)) - (-iYlim(1,2) - traceVertOffset(uniqColorsInds(ititleTag)))/2;
                                end
                                itxt = text(iripWinTimes(1)-.15, titleY, ['{\color[rgb]' sprintf('{%d %d %d} %s %s}',titleColors(ititleTag,:), titleAreas{ititleTag}, num2str(titleSubAreas{ititleTag}))]);
                                set(itxt, 'rotation', AreaLabelRotate, 'FontName', 'Arial', 'FontSize',12, 'FontWeight', 'Bold');
                            end
%                             supertitle([{filenameTitle};{['\fontsize{16}black {\color{magenta}magenta ' '\color[rgb]{0 .5 .5}teal \color{red}red} black again']}])
                            iStitle = supertitle(stitlename);
                            set(iStitle, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial');
                            set(gca,'children',flipud(get(gca,'children'))) %send the patch behind the LFP traces
                            ylim([iYlim(1,2) iYlim(1,1)])
                            xlim([iripWinTimes(1) iripWinTimes(end)])
                            set(gca, 'YTick', []);
                            set(gca,'XTick',[iripWinTimes(1):(iripWinTimes(end) - iripWinTimes(1))/10:iripWinTimes(end)], 'FontSize',10,'FontWeight','normal')
                            set(gca, 'XTickLabel', round([iripWinTimes(1):(iripWinTimes(end) - iripWinTimes(1))/10:iripWinTimes(end)],2), 'XTickLabelRotation',45);
%                             set(gca, 'XTickLabel', [-win(1):win(1)/5:win(2)])
                            xlabel('time(s)','FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial')
                            return
%                             if pausefigs
%                                 pause
%                             end
                        end
                    end
                end
            end
        end
    end
end


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