
%% Dashboard
close all
runFilterFramework = 1;
    ; saveFilterOutput = runFilterFramework;
    ; loadFilterOutput = 0;
processnTrodeMean = 0;
processAreaMean = 0;
plotfigs = 0;
    ; plotLFPtraces = 0;
    ; plotnTrodeMeanLFPtraces = 0;
    ; plotAreaMeanLFPtraces = 0;
    ; savefigs = 0;
    ; pausefigs = 0;
%% ---------------- plotting params --------------------------
% colorsMEC = cbrewer('seq', 'Blues', 10, 'PCHIP');
% colorsCTX = cbrewer('seq', 'Reds', 10, 'PCHIP');
% allthecolors = {[0 0 0], colorsMEC, colorsCTX};
% usecolormap = 'jet'; %colorcube %lines %jet winter
colorSet = 'DR1';
win = [1 1]; %in seconds
% indwin = win*1500;
plotwin = [.5 .5]; %in seconds
plotHilbertPower = 0;

% indwin = plotwin*1500;
%% ---------------- Data Filters --------------------------
% animals = {'D10'};
% days = [1:12]; % 1:12 
% days = [1:12]; % 1:12
% animals = {'D12'};
% days = [1:7]; % 1:7
% animals = {'D13'};
% days = [1:10];
animals = {'JZ1'};
days = [1]; %1:14 1:9
%  animals = {'JZ2'};
%  days = [1:5];
% animals = {'JZ3'};
% days = [1 2 5:10]; %[1 2 5:10]
% animals = {'JZ4'};
% days = [1:10];%1:10
filtfunction = 'riptriglfp';
LFPtypes = {'eeg', 'ripple', 'theta', 'lowgamma', 'fastgamma'}; %
LFPrangesHz = {'1-400', '140-250', '6-9', '20-50', '65 - 140'}; %need to make this a lookup from the filter mat
eventtype = 'rippleskons';
% eventtype = 'noisekons';
% eventarea = 'ca1';
epochType = {'sleep'};
epochEnvironment = {'sleep'};% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
eventSourceArea = 'ca1';
ripAreas = {'ca1', 'mec', 'por'};
ntAreas = {'ca1','sub','mec', 'por', 'v2l', 'ref'};  %'sub', 
plotenv = {'sleep'};

consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 7;        % STD. how big your ripples are
exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
minvelocity = 0;
maxvelocity = 4;
%% ---------------- Paths and Title strings ---------------------------------------------------
investInfo = animaldef(lower('Demetris'));
filtOutputDirectory = sprintf('%s%s/', investInfo{2}, filtfunction);
figdirectory = sprintf('%s%s/', investInfo{4}, filtfunction);
filenamesave = sprintf('%s%sSTD%d_%s_%s_%s_D%s', eventSourceArea, eventtype, minstdthresh, strjoin(epochEnvironment,'-'), cell2mat(animals), cell2mat(LFPtypes),strjoin(arrayfun(@(x) num2str(x),days,'UniformOutput',false),'-'));
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
filenameTitle = strrep(filename,'_', '\_');
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1;
    %     epochfilter =    sprintf('(isequal($type, ''%s'')) && (isequal($environment, ''%s''))',epochType, epochEnvironment);
    eptypeEnv = [epochType; epochEnvironment];
    epochfilter = sprintf('((isequal($type, ''%s'')) && (isequal($environment, ''%s''))) || ',eptypeEnv{:});
    yuck = strfind(epochfilter, '||');
    epochfilter = epochfilter(1:yuck(end)-1);  %cut off the trailing '||'
    iterator = 'multitetrodeanal'; %multitetrodeanal
    %tetfilter: the ntrodeID's that pass this filter get stored into f.eegdata{day}{epoch}[ntrodeIDs]
    %     tetfilter = sprintf('(isequal($area,''%s'') || isequal($area,''%s'') || isequal($area,''%s'') || isequal($area,''%s'') || isequal($area,''%s'') || isequal($area,''%s''))', ntAreas{1}, ntAreas{2}, ntAreas{3}, ntAreas{4}, ntAreas{5}, ntAreas{6});
    tetfilter = sprintf('(isequal($area,''%s'')) || ', ntAreas{:});
    yuck = strfind(tetfilter, '||');
    tetfilter = tetfilter(1:yuck(end)-1); %cut off the trailing '||'
    % timefilter{1} = {'get2dstate','($velocity<4)'};
    % timefilter{2} = {'kk_getriptimes','($nripples>=1)',[],'tetfilter',tetfilter,'minthresh',5};
    timefilter{1} = {'getconstimes', '($cons == 1)', [eventSourceArea eventtype],1,'consensus_numtets',consensus_numtets,...
        'minstdthresh',minstdthresh,'exclusion_dur',exclusion_dur,'minvelocity',minvelocity,'maxvelocity',maxvelocity};
%     if strcmp(animals{1}, 'D12')
        timefilter{2} = {'excludenoiseevents', '($noise == 0)', [eventSourceArea,'noisekons'], 1, };
%     end
    %---------- save all filtering parameters in workspace into struct
    %     iF.datafilter = whos;
    %----------F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);--------
    F = createfilter('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);
    %----------f = setfilteriterator(f, funcname, loadvariables, options)--------
    eval(['F = setfilterfunction(F, [''dfa_'' filtfunction], {[eventSourceArea eventtype],' ...
        strjoin(arrayfun(@(x) sprintf('''%s'',', cell2mat(x)), LFPtypes,'UniformOutput', ...
        false)) '},''eventtype'', [eventSourceArea eventtype], ''LFPtypes'', LFPtypes, ...''win'', win);']);
    %     eval(['F = setfilterfunction(F, [''dfa_'' filtfunction], {[eventSourceArea eventtype],' strjoin(reshape(repmat(arrayfun(@(x) sprintf('''%s'',', cell2mat(x)), LFPtypes,'UniformOutput',false), 2, 1), 1,length(LFPtypes)*2)) '},''eventtype'', [eventSourceArea eventtype], ''LFPtypes'', LFPtypes);']);
    tic
    F = runfilter(F);
    F(1).filterTimer = toc; F(1).filterTimer
    F(1).worldDateTime = clock;
    F(1).dataFilter = struct('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator, 'filename', filename);
end
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    if ~isdir(filtOutputDirectory);
        mkdir(filtOutputDirectory);
    end
    %save the entire workspace for filter provenance
    save(sprintf('%s/%s',filtOutputDirectory, filename), 'F','-v7.3');
    fprintf('filteroutput saved to %s/%s\n',filtOutputDirectory, filename)
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s/%s',filtOutputDirectory, filename))
    fprintf('filteroutput loaded: %s/%s\n',filtOutputDirectory, filename)
end
%% plot
if plotfigs
    %% plot LFP traces per ntrode, per ripple
    if plotLFPtraces
        %% STEP 2: For each EEG + Ripple Type... select nTrodes and Lookup Rip time windows
        % For each LFP type (wideband EEG, ripple band, etc), load all of the regions LFP files into eegstruct
%         for iLFPtype = 1:length(LFPtypes); 
            for ianimal = 1:length(F)
                    %% ---- loadtetinfostruct ----
                    animalinfo = animaldef(lower(animals{ianimal}));
                    animalID = animalinfo{1,3}; %use anim prefix for name
                    FFanimdir =  sprintf('%s',animalinfo{1,2});
                    load([FFanimdir, animalID, 'tetinfo']);
                    tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
                    days = find(cell2mat(cellfun(@(x) ~isempty(x), F.output, 'un', 0)));
                    for iday = 1:length(days);
                        day = days(iday);
                        %% ---- load the ripplekons file for each area for this day ----
                        for iareakons = 1:length(ripAreas);
                            load(sprintf('%s%s%s%s%02d.mat', FFanimdir, animalID, ripAreas{iareakons}, eventtype, day));
                            eval([sprintf('eventKons.%s%s = %s%s;', ripAreas{iareakons}, eventtype, ripAreas{iareakons}, eventtype)])
                            eval([sprintf('clear %s%s', ripAreas{iareakons}, eventtype)])
                        end
                        %add filtering via getconstimes to get ripples of
                        %certain strength
                        epochs = find(cell2mat(cellfun(@(x) ~isempty(x), {F(ianimal).output{day}(:).data}, 'un', 0)));
                        for iepoch = 1:length(epochs)
                            epoch = epochs(iepoch);
                            if isempty(F(ianimal).output{day}(epoch).data)
                                continue %if this an day epoch is empty, skip to the next
                            end
                            %% ---- get data info for rips within this day epoch ----
                            iepochLFPtimes = []; ripStartIndices = []; ripEndIndices = []; ntrodesIndices = []; %plotwin = [];
                            iepochLFPtimes = F(ianimal).output{day}(epoch).LFPtimes;
                            ripStartIndices = F(ianimal).output{day}(epoch).eventStartIndices;
                            ripEndIndices = F(ianimal).output{day}(epoch).eventEndIndices;
%                             plotwin = F(ianimal).output{day}(epoch).win;
                            indwin = plotwin*1500;
                            ntrodesIndices = F(ianimal).output{day}(epoch).index;
                            %% ---- reorder the LFP traces by suparea|area|subarea tags (in that priority) ----
                            [~, tagIndMap] = ismember(ntrodesIndices,tetinfoAll.index, 'rows');
                            ntrodeTags = tetinfoAll.values(tagIndMap);
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
                            icolors = colorPicker(colorSet, strAreas, strSubAreas);
                            numsumallareatags = cell2mat([numsumSupAreas numsumAreas numsumSubAreas]);
                            [numsumallSort, numsumSortInds] = sortrows(numsumallareatags);%,[-1 -2 -3]); % -Col to 'descend'
                            icolors = icolors(numsumSortInds,:);
                            %% ---- loop across all ripples for this day epoch sourcearea ----
                            for irip = 1:length(F(ianimal).output{day}(epoch).data{1})
                                %if you want to save figs but not pause them to take a look.. just keep them invisible. 
                                % i think this makes it faster and returns keyboard/mouse control during saving
                                if savefigs && ~pausefigs; 
                                    ifig = figure('Visible','off','units','normalized','position',[.1 .1 .9 .8]);
                                else
                                    ifig = figure('units','normalized','position',[.1 .1 .9 .8]);
                                end
                                %% ---- Colors ----
                                set(gcf,'color','white')
                                %                             clrmat = unique(colormap(usecolormap), 'rows', 'stable'); %flipud
                                %                             clrmat = clrmat(all(clrmat < 1,2),:); % only use dark/bold colors
                                %                             colorslength = length(clrmat(:,1));
                                iripWinTimes = iepochLFPtimes(ripStartIndices(irip)-indwin: ripStartIndices((irip))+indwin);
                                %                             icolor = mod(allareatagsSorted(:,:), colorslength);
                                %                             icolor = 1-(icolor-min(min(icolor)))./max(max(icolor)); %normalize
                                
                                
                                %                             icolors = clrmat(mod(prod(numsumallSort,2), colorslength)+1,:);   % the +1 is to assure all indices > 0 (bc matlab is 1-based ..)
                                %                             [uniqColors, uniqColorsInds, uniqColorsInds2] = unique(icolors, 'rows', 'stable');
                                %                             [uniqNumsum, uniqNumSumInds, uniqNumSumInds2] = unique(numsumallSort, 'rows', 'stable');
                                %                             titleColors = icolors(uniqColorsInds,:);
                                %                             titleSupAreas = strSupAreas(numsumSortInds(uniqColorsInds));
                                %                             titleAreas = strAreas(numsumSortInds(uniqColorsInds));
                                %                             titleSubAreas = strSubAreas(numsumSortInds(uniqColorsInds));
                                iYlim = [];
                                %% ---- PLOT THE POWERTRACE OF THE CONSCENSUS RIPPLES FILTERED LFP ----
                                %                             eval([sprintf(' iconsensusTimes = eventKons.%s%s{%d}{%d}{1}.eegtimesvec_ref;', eventSourceArea, eventtype, iday, iepoch)])
                                %                             iconsensusInds = find((iconsensusTimes == iripWinTimes(1) | iconsensusTimes == iripWinTimes(end))); %if the rip start time or end time is within the current window
                                %                             eval([sprintf(' iconsensusTrace = eventKons.%s%s{%d}{%d}{1}.powertrace(%d:%d);', eventSourceArea, eventtype, iday, iepoch, iconsensusInds(1), iconsensusInds(2))])
                                %                             conTraceAx = subaxis(2,2,2, 'Spacing', 0.02, 'Padding', 0.0, 'Margin', 0.09);
                                %                             plot(iripWinTimes, iconsensusTrace, 'Color', icolors(end,:), 'Linewidth', 2)
                                %                             hold on;
                                %% ---- loop over each LFP type ----
                                for iLFPtype = 1:length(F(ianimal).output{day}(epoch).LFPtypes)
                                    for introde = 1:length(F(ianimal).output{day}(epoch).data{iLFPtype}{irip}(:,1))
                                        introdeID = F(ianimal).output{day}(epoch).index(numsumSortInds(introde), 3);
                                        isupareatag = ntrodeTags{numsumSortInds(introde)}.suparea;
                                        iareatag = ntrodeTags{numsumSortInds(introde)}.area;
                                        isubareatag = ntrodeTags{numsumSortInds(introde)}.subarea;
                                        %                                     clrpicker = find(strcmp(iareatag, ntAreas));
                                        %                                     allthecolors{clrpicker};
                                        iNTcolor = icolors(introde,:);
                                        subaxis(1,length(F(ianimal).output{day}(epoch).data),iLFPtype, 'Spacing', 0.02, 'Padding', 0.0, 'MarginLeft', 0.04, 'MarginRight', 0.02, 'MarginTop', 0.09, 'MarginBottom', 0.09);
                                        %                                     subaxis(2,2,iLFPtype+2, 'Spacing', 0.02, 'Padding', 0.0, 'Margin', 0.09);
                                        %                                     subplot(1,length(F(ianimal).output{iday}(iepoch).data),iLFPtype)
                                        introdeiripLFP = double(F(ianimal).output{day}(epoch).data{iLFPtype}{irip}(numsumSortInds(introde),:));
                                        middleind = ceil(length(introdeiripLFP)/2);
                                        introdeiripLFP = introdeiripLFP(middleind-(plotwin(1)*1500):middleind+(plotwin(2)*1500));
                                        % ---- plot ntrode traces ----
                                        if introde == 1; %if it's the first ntrode
                                            plot(iripWinTimes, introdeiripLFP, 'Color', iNTcolor, 'Linewidth', 1)
                                            hold on;
                                            traceVertOffset(introde) = 0;
                                            previousTrace = introdeiripLFP;
                                            iYlim(1,1) = max(introdeiripLFP); %first trace max
                                        else
                                            traceVertOffset(introde) = traceVertOffset(introde-1) + abs(min(previousTrace)) + abs(max(introdeiripLFP));
                                            plot(iripWinTimes, introdeiripLFP - traceVertOffset(introde), 'Color', iNTcolor, 'LineWidth',1);
                                            hold on;
                                            previousTrace = introdeiripLFP;
                                        end
                                        %% ---- set ripple patches for this ntrode ----
                                        if any(strcmp(iareatag,ripAreas));
                                            %get index of matching region iareakons
                                            eval([sprintf(' ripsStartTimes = eventKons.%s%s{%d}{%d}{1}.starttime;', iareatag, eventtype, day, epoch)])
                                            eval([sprintf(' ripsEndTimes = eventKons.%s%s{%d}{%d}{1}.endtime;', iareatag, eventtype, day, epoch)])
                                            
                                            %                                     ripsStartTimes = eventKons{1}.ca1rippleskons{iday}{iepoch}{1}.starttime;
                                            %                                     ripsEndTimes = eventKons{1}.ca1rippleskons{iday}{iepoch}{1}.endtime;
                                            ripsinwinInds = find(ripsStartTimes>iripWinTimes(1) & ripsEndTimes<iripWinTimes(end)); %if the rip start time and end time is within the current window
                                            ripsinwinTimes = [ripsStartTimes(ripsinwinInds) ripsEndTimes(ripsinwinInds)];
                                            
                                            
                                            
                                            %                                     if (introde == 1) && (iareatag == eventSourceArea);
                                            %                                         %PLOT THE POWERTRACE OF THE CONSCENSUS RIPPLES FILTERED LFP
                                            % %                                         srcarea find(strcmp(eventSourceArea, tetAreas));
                                            %                                         eval([sprintf(' iconsensusTimes = eventKons.%s%s{%d}{%d}{1}.eegtimesvec_ref;', iareatag, eventtype, iday, iepoch)])
                                            %                                         iconsensusInds = find((iconsensusTimes == iripWinTimes(1) | iconsensusTimes == iripWinTimes(end))); %if the rip start time or end time is within the current window
                                            %                                         eval([sprintf(' iconsensusTrace = eventKons.%s%s{%d}{%d}{1}.powertrace(%d:%d);', iareatag, eventtype, iday, iepoch, iconsensusInds(1), iconsensusInds(2))])
                                            %                                     end
                                            for iripinwin = 1:length(ripsinwinTimes(:,1))
                                                ilmax = max(introdeiripLFP - traceVertOffset(introde));
                                                ilmin = min(introdeiripLFP - traceVertOffset(introde));
                                                iripstarttime = ripsStartTimes(ripsinwinInds(iripinwin));
                                                iripendtime = ripsEndTimes(ripsinwinInds(iripinwin));
                                                %                                     line([iripstarttime iripstarttime], [ilmin ilmax],'Color',[.8 .8 .8],'LineWidth',1.5);
                                                Xpatch = [iripstarttime iripendtime iripendtime iripstarttime];
                                                Ypatch = [ilmin ilmin ilmax ilmax];
                                                patch(Xpatch, Ypatch, iNTcolor, 'FaceAlpha', .25, 'edgecolor','none'); %triggering-ripple patch
                                                %% ---- plot the maxthresh scores for all visible ripples ----
                                                if ~isempty(strfind(eventtype, LFPtypes{iLFPtype})) %iLFPtype == length(F(ianimal).output{iday}(iepoch).LFPtypes) %if it's the rightmost subplot, plot the maxthresh scores for all visible ripples
                                                    imaxthresh = [];
                                                    eval([sprintf(' imaxthresh = eventKons.%s%s{%d}{%d}{1}.maxthresh(%d);', iareatag, eventtype, day, epoch, ripsinwinInds(iripinwin))])
                                                    iriplabel = text(double(iripendtime),double(ilmax-(ilmax-ilmin)/10),num2str(round(imaxthresh,1)));
                                                    set(iriplabel, 'Color', iNTcolor, 'FontName', 'Arial', 'FontSize',6, 'FontWeight', 'normal', 'horizontalAlignment', 'left', 'verticalAlignment', 'bottom');
                                                end
                                            end
                                        end
                                        NTlabel = text(iripWinTimes(1)-.01, -traceVertOffset(introde), num2str([introdeID]));
                                        set(NTlabel, 'Color', [.5 .5 .5], 'FontName', 'Arial', 'FontSize',8, 'FontWeight', 'normal', 'horizontalAlignment', 'right');
                                        if iLFPtype == 1
                                            NTlabel = text(iripWinTimes(1)-.07, -traceVertOffset(introde), num2str([iareatag, ' ', num2str(isubareatag)]));
                                            set(NTlabel, 'Color', iNTcolor, 'FontName', 'Arial', 'FontSize',8, 'FontWeight', 'normal', 'horizontalAlignment', 'right');
                                        end
                                        %                                     if iareatag == eventSourceArea;
                                        %                                         sourceAreaColor = iNTcolor;
                                        %                                     end
                                    end
                                    %% ---- Source ripple line, subplot title ----
                                    %                                 set(conTraceAx,'Color', sourceAreaColor) %set(subplot(2,2,1),'Color','Red')
                                    iYlim(1,2) = min(introdeiripLFP - traceVertOffset(introde)); %last trace min
                                    iripstarttime = iepochLFPtimes(ripStartIndices(irip));
                                    %                                 iripendtime = iepochLFPtimes(ripEndIndices(irip));
                                    Xline = [iripstarttime iripstarttime];
                                    Yline = [[iYlim(1,2) iYlim(1,1)]];
                                    line(Xline, Yline, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1.5);
                                    
                                    ititle = title(sprintf('%s %s Hz',F(ianimal).output{day}(epoch).LFPtypes{iLFPtype}, LFPrangesHz{iLFPtype}));
                                    set(ititle,'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial');
                                    set(gca,'children',flipud(get(gca,'children'))) %send the patch behind the LFP traces
                                    ylim([iYlim(1,2) iYlim(1,1)])
                                    xlim([iripWinTimes(1) iripWinTimes(end)])
                                    set(gca, 'YTick', []);
                                    set(gca,'XTick',[iripWinTimes(1):(iripWinTimes(end) - iripWinTimes(1))/10:iripWinTimes(end)], 'FontSize',8,'FontWeight','normal')
                                    set(gca, 'XTickLabel', round([iripWinTimes(1):(iripWinTimes(end) - iripWinTimes(1))/10:iripWinTimes(end)],2), 'XTickLabelRotation',30);
                                    %                             set(gca, 'XTickLabel', [-win(1):win(1)/5:win(2)])
                                    xlabel('time(s)','FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial')
                                    %% ---- area axes label ----
                                    %                                 if iLFPtype == 1
                                    %                                     for ititleTag = 1:length(titleColors(:,1))
                                    %                                         try %place in the middle of this area's traces
                                    %                                             titleY = -traceVertOffset(uniqColorsInds(ititleTag)) - (traceVertOffset(uniqColorsInds(ititleTag+1)) - traceVertOffset(uniqColorsInds(ititleTag)))/2;
                                    %                                         catch %will catch when on the last title, at this point, place between the current vert offset and total minimum
                                    %                                             titleY = -traceVertOffset(uniqColorsInds(ititleTag)) - (-iYlim(1,2) - traceVertOffset(uniqColorsInds(ititleTag)))/2;
                                    %                                         end
                                    %                                         itxt = text(iripWinTimes(1)-.05, titleY, ['{\color[rgb]' sprintf('{%d %d %d} %s %s}',titleColors(ititleTag,:), titleAreas{ititleTag}, num2str(titleSubAreas{ititleTag}))]);
                                    %                                         %                                         set(itxt, 'rotation', AreaLabelRotate, 'FontName', 'Arial', 'FontSize',12, 'FontWeight', 'Bold');
                                    %                                         set(itxt, 'rotation', 0, 'FontName', 'Arial', 'FontSize',12, 'FontWeight', 'Bold', 'horizontalAlignment', 'right', 'verticalAlignment', 'bottom');
                                    %                                     end
                                    %                                 end
                                    %% ---- super title ----
                                end
                                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                                sprtit = sprintf('%s D%d E%d %s%d', animalID, day, epoch, eventtype, irip);
                                iStitle = text(.5, .95, [{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]], 'Parent', sprtitleax, 'Units', 'normalized');
                                %                             iStitle = supertitle([{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]]);
                                set(iStitle, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
                                %% ---- pause, save figs ----
                                if pausefigs
                                    pause
                                end
                                if savefigs
                                    if ~isdir(figdirectory);
                                        mkdir(figdirectory);
                                    end
                                    if ~isdir([figdirectory filenamesave]);
                                        mkdir([figdirectory filenamesave]);
                                    end
                                    sprtitsave = strrep(sprtit,' ', '_'); %put '_' character in place of whitespace for filename
                                    set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
                                    currfigfile = sprintf('%s%s/%s',figdirectory, filenamesave, sprtitsave);
                                    print(currfigfile,'-dpng', '-r0')
                                    disp(sprintf('plot %s saved', sprtit))
                                end
                                close all
                            end
                        end
                    end
                end
%         end
    end
    %% plot lfp traces per ntrode, all ripples
    if plotnTrodeMeanLFPtraces
        for ianimal = 1:length(F)
            %% ---- loadtetinfostruct ----
            animalinfo = animaldef(lower(animals{ianimal}));
            animalID = animalinfo{1,3}; %use anim prefix for name
            FFanimdir =  sprintf('%s',animalinfo{1,2});
            load([FFanimdir, animalID, 'tetinfo']);
            tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
            days = find(~cellfun(@isempty, F.output));
            for day = days; % 1:length(F(ianimal).output)
                %% ---------- combine all rip trig traces across epochs for each env type --------------------------------
                load([FFanimdir, animalID, 'task', sprintf('%02d.mat', day)]);
                epochs = find(~arrayfun(@isempty, F.output{day}));
                %% ---------- get unique epoch environment types  --------------------------------
                eptypes = [];
                for epoch = epochs;
                    eptypes{epoch,1} = task{day}{epoch}.environment;
                end
                [envtypes,IndA,IndC] = unique(eptypes, 'stable');
                
                for iepochEnv = 1:size(envtypes,1); %for each envtype. cat3 all rip trig traces across epochs
                    if isempty(strfind(envtypes{iepochEnv}, plotenv))
                        continue
                    end
                    ienvTypeInds = find(iepochEnv == IndC);
                    numLFPtypes = length(F(ianimal).output{1, day}(ienvTypeInds(1)).LFPtypes);
                    if processnTrodeMean
                        for iLFPtype = 1:numLFPtypes;
                            try
                            epienvrips = arrayfun(@(x) x.data{iLFPtype}, F(ianimal).output{1, day}(ienvTypeInds),'un', 0);
                            catch
                                continue
                            end
                            ienvTallrips = [epienvrips{:}];
                            ienvTallripsCat = double(cat(3,ienvTallrips{:}));
                            if plotHilbertPower
                                %for each event, calculate the hilbert transform for each ntrode.
                                % make sure the time domain is the first dimension... as this is the only way the hilbert works with a matrix input
                                ienvTallHilbCat = nan(size(ienvTallripsCat));
                                for ievent = 1:size(ienvTallripsCat,3)
                                    tmplfp = squeeze(ienvTallripsCat(:,:,ievent))';
                                    tmphilb = abs(hilbert(tmplfp));
                                    ienvTallHilbCat(:,:,ievent) = tmphilb';
                                end
                                %                             temppow  = nanmean(ienvTallHilbCat.^2,3);
                                %                             pow = 10*log10(temppow./nanmean(temppow(:,1:750),2));
                                %                             pow = 10*log10(inv(diag(nanmean(temppow(:,1:750),2)))*temppow);
                                
                                ienvTallripsMean{iLFPtype} = nanmean(ienvTallHilbCat,3);
                                ienvTallripsStd{iLFPtype} = std(ienvTallHilbCat,0,3,'omitnan');
                            else
                                ienvTallripsMean{iLFPtype} = nanmean(ienvTallripsCat,3);
                                ienvTallripsStd{iLFPtype} = std(ienvTallripsCat,0,3,'omitnan');
                            end
                            %                     ienvTallripsSem{iLFPtype} = std(ienvTallripsCat,0,3,'omitnan') / sqrt(size(ienvTallripsCat,3));
                            %                 imagesc(ienvTallripsMean{iLFPtype})
                            %                 imagesc(ienvTallripsStd{iLFPtype})
                        end
                    end
%                     plotwin = F(ianimal).output{1, day}(ienvTypeInds(1)).win;
                    iepochLFPtimes = [-plotwin:1/1500:plotwin]; % lfp time for average plots is just the window length
                    ntrodesIndices = [];
                    useep = 1;
                    while isempty(ntrodesIndices)
                        ntrodesIndices = F(ianimal).output{day}(ienvTypeInds(useep)).index;
                        useep = useep + 1;
                        if useep > length(ienvTypeInds)
                            continue
                        end
                    end
                    %% ---- reorder the LFP traces by suparea|area|subarea tags (in that priority) ----
                    [~, tagIndMap] = ismember(ntrodesIndices,tetinfoAll.index, 'rows');
                    ntrodeTags = tetinfoAll.values(tagIndMap);
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
                    icolors = colorPicker(colorSet, strAreas, strSubAreas);
                    numsumallareatags = cell2mat([numsumSupAreas numsumAreas numsumSubAreas]);
                    [numsumallSort, numsumSortInds] = sortrows(numsumallareatags);%,[-1 -2 -3]); % -Col to 'descend'
                    icolors = icolors(numsumSortInds,:);
                    %% ---- PLOT ----
                    for irip = 1;%:length(F(ianimal).output{iday}(iepoch).data{1})
                        if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
                            ifig = figure('Visible','off','units','normalized','position',[.1 .1 .9 .8]);
                        else
                            ifig = figure('units','normalized','position',[.1 .1 .9 .8]);
                        end
                        %% ---- Colors ----
                        set(gcf,'color','white')
                        %                     iripWinTimes = iepochLFPtimes(ripStartIndices(irip)-indwin: ripStartIndices((irip))+indwin);
                        iripWinTimes = iepochLFPtimes; % lfp time for average plots is just the window length
                        iYlim = [];
                        %% ---- loop over each LFP type ----
                        for iLFPtype = 1:numLFPtypes;
                            for introde = 1:size(ienvTallripsMean{iLFPtype},1)
                                introdeID = F(ianimal).output{day}(ienvTypeInds(1)).index(numsumSortInds(introde), 3);
                                isupareatag = ntrodeTags{numsumSortInds(introde)}.suparea;
                                iareatag = ntrodeTags{numsumSortInds(introde)}.area;
                                isubareatag = ntrodeTags{numsumSortInds(introde)}.subarea;
                                iNTcolor = icolors(introde,:);
                                subaxis(1,numLFPtypes,iLFPtype, 'Spacing', 0.02, 'Padding', 0.0, 'MarginLeft', 0.04, 'MarginRight', 0.02, 'MarginTop', 0.09, 'MarginBottom', 0.09);
                                introdeiripLFP = ienvTallripsMean{iLFPtype}(numsumSortInds(introde),:);
                                middleind = ceil(length(introdeiripLFP)/2);
                                introdeiripLFP = introdeiripLFP(middleind-(plotwin(1)*1500):middleind+(plotwin(2)*1500));
                                introdeiripLFPStd = ienvTallripsStd{iLFPtype}(numsumSortInds(introde),:);
                                %                             introdeiripLFP = double(F(ianimal).output{iday}(iepoch).data{iLFPtype}{irip}(numsumSortInds(introde),:));
                                % ---- plot ntrode traces ----
                                if introde == 1; %if it's the first ntrode
                                    %                                 plot(iripWinTimes, introdeiripLFP, 'Color', iNTcolor, 'Linewidth', 1);
                                    %                                 hold on;
                                    traceVertOffset(introde) = 0;
                                    previousTrace = introdeiripLFP;
                                    iYlim(1,1) = max(introdeiripLFP); %first trace max
                                else
                                    traceVertOffset(introde) = traceVertOffset(introde-1) + abs(min(previousTrace)) + abs(max(introdeiripLFP));
                                    
                                    %                                 plot(iripWinTimes, introdeiripLFP - traceVertOffset(introde), 'Color', iNTcolor, 'LineWidth',1);
                                    %                                 hold on;
                                    previousTrace = introdeiripLFP;
                                    introdeiripLFP = introdeiripLFP - traceVertOffset(introde);
                                end
                                %                             shadedErrorBar(iripWinTimes,introdeiripLFP,introdeiripLFPStd,{'Color', iNTcolor, 'Linewidth', 1},1);
                                lh = plot(iripWinTimes, introdeiripLFP, 'Linewidth', 1.2);
                                lh.Color = [iNTcolor];
                                
                                hold on;                          %
                                
                                NTlabel = text(iripWinTimes(1)-.01, -traceVertOffset(introde), num2str([introdeID]));
                                set(NTlabel, 'Color', [.5 .5 .5], 'FontName', 'Arial', 'FontSize',8, 'FontWeight', 'normal', 'horizontalAlignment', 'right');
                                if iLFPtype == 1
                                    NTlabel = text(iripWinTimes(1)-.07, -traceVertOffset(introde), num2str([iareatag, ' ', num2str(isubareatag)]));
                                    set(NTlabel, 'Color', iNTcolor, 'FontName', 'Arial', 'FontSize',8, 'FontWeight', 'normal', 'horizontalAlignment', 'right');
                                end
                            end
                            %% ---- Source ripple line, subplot title ----
                            %                                 set(conTraceAx,'Color', sourceAreaColor) %set(subplot(2,2,1),'Color','Red')
                            iYlim(1,2) = min(introdeiripLFP); % - traceVertOffset(introde)); %last trace min
                            iripstarttime = 0; %iepochLFPtimes(ripStartIndices(irip));
                            %                                 iripendtime = iepochLFPtimes(ripEndIndices(irip));
                            Xline = [iripstarttime iripstarttime];
                            Yline = [[iYlim(1,2) iYlim(1,1)]];
                            line(Xline, Yline, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1.5);
                            
                            ititle = title(sprintf('%s %s Hz',F(ianimal).output{day}(epoch).LFPtypes{iLFPtype}, LFPrangesHz{iLFPtype}));
                            set(ititle,'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial');
                            set(gca,'children',flipud(get(gca,'children'))) %send the patch behind the LFP traces
                            ylim([iYlim(1,2) iYlim(1,1)])
                            %                         ylim auto
                            xlim([iripWinTimes(1) iripWinTimes(end)])
                            set(gca, 'YTick', []);
                            set(gca,'XTick',[iripWinTimes(1):(iripWinTimes(end) - iripWinTimes(1))/10:iripWinTimes(end)], 'FontSize',8,'FontWeight','normal')
                            set(gca, 'XTickLabel', round([iripWinTimes(1):(iripWinTimes(end) - iripWinTimes(1))/10:iripWinTimes(end)],2), 'XTickLabelRotation',30);
                            %                             set(gca, 'XTickLabel', [-win(1):win(1)/5:win(2)])
                            xlabel('time(s)','FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial')
                            %% ---- super title ----
                        end
                        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                        sprtit = sprintf('%s D%d %s %s STD>%d', animalID, day, eventtype, envtypes{iepochEnv}, minstdthresh);
                        iStitle = text(.5, .95, [{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]], 'Parent', sprtitleax, 'Units', 'normalized');
                        %                             iStitle = supertitle([{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]]);
                        set(iStitle, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
                        %% ---- pause, save figs ----
                        if pausefigs
                            pause
                        end
                        if savefigs
                            if ~isdir(figdirectory);
                                mkdir(figdirectory);
                            end
                            if plotHilbertPower
                                if ~isdir([figdirectory filenamesave '/ntrode_mean/hilbert']);
                                    mkdir([figdirectory filenamesave '/ntrode_mean/hilbert']);
                                end
                            elseif ~isdir([figdirectory filenamesave '/ntrode_mean/']);
                                mkdir([figdirectory filenamesave '/ntrode_mean/']);
                            end
                            sprtitsave = strrep(sprtit,' ', '_'); %put '_' character in place of whitespace for filename
                            set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
                            if plotHilbertPower
                                currfigfile = sprintf('%s%s/ntrode_mean/hilbert/%s',figdirectory, filenamesave, sprtitsave);
                            else
                                currfigfile = sprintf('%s%s/ntrode_mean/%s',figdirectory, filenamesave, sprtitsave);
                            end
                            print(currfigfile,'-dpng', '-r0')
                            disp(sprintf('plot %s saved', sprtit))
                        end
                        close all
                    end
                end
            end
        end
    end
    %% plot lfp traces per Area/subarea, all ripples
    if plotAreaMeanLFPtraces
        for ianimal = 1:length(F)
            %% ---- loadtetinfostruct ----
            animalinfo = animaldef(lower(animals{ianimal}));
            animalID = animalinfo{1,3}; %use anim prefix for name
            FFanimdir =  sprintf('%s',animalinfo{1,2});
            load([FFanimdir, animalID, 'tetinfo']);
            tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
            days = find(~cellfun(@isempty, F.output));
            for day = days; % 1:length(F(ianimal).output)
                %% ---------- combine all rip trig traces across epochs for each env type --------------------------------
                load([FFanimdir, animalID, 'task', sprintf('%02d.mat', day)]);
                epochs = find(~arrayfun(@isempty, F.output{day}));
                %% ---------- get unique epoch environment types  --------------------------------
                eptypes = [];
                for epoch = epochs;
                    eptypes{epoch,1} = task{day}{epoch}.environment;
                end
                [envtypes,IndA,IndC] = unique(eptypes, 'stable');
                for iepochEnv = 1:size(envtypes,1); %for each envtype. cat3 all rip trig traces across epochs
                    ienvTypeInds = find(iepochEnv == IndC);
                    numLFPtypes = length(F(ianimal).output{1, day}(ienvTypeInds(1)).LFPtypes);
                    
                    ntrodesIndices = F(ianimal).output{day}(ienvTypeInds(1)).index;
                    %% ---- reorder the LFP traces by suparea|area|subarea tags (in that priority) ----
                    [~, tagIndMap] = ismember(ntrodesIndices,tetinfoAll.index, 'rows');
                    ntrodeTags = tetinfoAll.values(tagIndMap);
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
                    icolors = colorPicker(colorSet, strAreas, strSubAreas);
                    numsumallareatags = cell2mat([numsumSupAreas numsumAreas numsumSubAreas]);
                    [numsumallSort, numsumSortInds] = sortrows(numsumallareatags);%,[-1 -2 -3]); % -Col to 'descend'
                    %                 icolors = icolors(numsumSortInds,:);
                    numsumAreas = cell2mat([numsumSupAreas numsumAreas]);
                    [subAreas,sAindA,sAindC] = unique(numsumAreas,'stable','rows');
                    if processAreaMean; %toggle to compute
                        ienvTallripsMean = [];
                        for iLFPtype = 1:numLFPtypes;
                            epienvrips = arrayfun(@(x) x.data{iLFPtype}, F(ianimal).output{1, day}(ienvTypeInds),'un', 0);
                            ienvTallrips = [epienvrips{:}];
                            ienvTallripsCat = double(cat(3,ienvTallrips{:}));
                            for isA = 1:max(sAindC); %for each unique area
                                isAInds{iLFPtype}{isA} = find(sAindC == isA);
                                ienvTallripsAreaCat = reshape(permute(ienvTallripsCat(isAInds{iLFPtype}{isA},:,:), [2 1 3]),1,size(ienvTallripsCat,2),[]);
                                ienvTallripsMean{iLFPtype}(isA,:) = nanmean(ienvTallripsAreaCat,3);
                                ienvTallripsStd{iLFPtype}(isA,:) = std(ienvTallripsAreaCat,0,3,'omitnan');
                            end
                        end
                    end
                    
%                     plotwin = F(ianimal).output{1, day}(ienvTypeInds(1)).win;
                    iepochLFPtimes = [-plotwin:1/1500:plotwin]; % lfp time for average plots is just the window length
                    %% ---- PLOT ----
                    for irip = 1;%:length(F(ianimal).output{iday}(iepoch).data{1})
                        if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
                            ifig = figure('Visible','off','units','normalized','position',[.1 .1 .9 .8]);
                        else
                            ifig = figure('units','normalized','position',[.1 .1 .9 .8]);
                        end
                        %% ---- Colors ----
                        set(gcf,'color','white')
                        %                     iripWinTimes = iepochLFPtimes(ripStartIndices(irip)-indwin: ripStartIndices((irip))+indwin);
                        iripWinTimes = iepochLFPtimes; % lfp time for average plots is just the window length
                        iYlim = [];
                        %% ---- loop over each LFP type ----
                        for iLFPtype = 1:numLFPtypes;
                            for introde = 1:size(ienvTallripsMean{iLFPtype},1)
                                introdeID = isAInds{iLFPtype}{introde};
                                unsortntrodeInd = find(isAInds{iLFPtype}{introde}(1) == numsumSortInds, 1);
                                isupareatag = ntrodeTags{isAInds{iLFPtype}{introde}(1)}.suparea;
                                iareatag = ntrodeTags{isAInds{iLFPtype}{introde}(1)}.area;
                                %                             isubareatag = ntrodeTags{numsumSortInds(introde)}.subarea;
                                iNTcolor = icolors(isAInds{iLFPtype}{introde}(1),:);
                                subaxis(1,numLFPtypes,iLFPtype, 'Spacing', 0.02, 'Padding', 0.0, 'MarginLeft', 0.04, 'MarginRight', 0.02, 'MarginTop', 0.09, 'MarginBottom', 0.09);
                                introdeiripLFP = ienvTallripsMean{iLFPtype}(introde,:);
                                introdeiripLFPStd = ienvTallripsStd{iLFPtype}(introde,:);
                                %                             introdeiripLFP = double(F(ianimal).output{iday}(iepoch).data{iLFPtype}{irip}(numsumSortInds(introde),:));
                                % ---- plot ntrode traces ----
                                if introde == 1; %if it's the first ntrode
                                    %                                 plot(iripWinTimes, introdeiripLFP, 'Color', iNTcolor, 'Linewidth', 1);
                                    %                                 hold on;
                                    traceVertOffset(introde) = 0;
                                    previousTrace = introdeiripLFP;
                                    iYlim(1,1) = max(introdeiripLFP); %first trace max
                                else
                                    traceVertOffset(introde) = traceVertOffset(introde-1) + abs(min(previousTrace)) + abs(max(introdeiripLFP));
                                    
                                    %                                 plot(iripWinTimes, introdeiripLFP - traceVertOffset(introde), 'Color', iNTcolor, 'LineWidth',1);
                                    %                                 hold on;
                                    previousTrace = introdeiripLFP;
                                    introdeiripLFP = introdeiripLFP - traceVertOffset(introde);
                                end
                                %                             shadedErrorBar(iripWinTimes,introdeiripLFP,introdeiripLFPStd,{'Color', iNTcolor, 'Linewidth', 1},1);
                                %                             shadedErrorBar(iripWinTimes,introdeiripLFP,introdeiripLFPStd,{'Linewidth', 1},1);
                                lh = plot(iripWinTimes, introdeiripLFP, 'Linewidth', 1.2);
                                lh.Color = [iNTcolor, .7];
                                
                                hold on;                          %
                                
                                NTlabel = text(iripWinTimes(1)-.01, -traceVertOffset(introde), num2str([introdeID]));
                                set(NTlabel, 'Color', [.5 .5 .5], 'FontName', 'Arial', 'FontSize',8, 'FontWeight', 'normal', 'horizontalAlignment', 'right');
                                if iLFPtype == 1
                                    NTlabel = text(iripWinTimes(1)-.07, -traceVertOffset(introde), num2str([iareatag, ' ']));
                                    set(NTlabel, 'Color', iNTcolor, 'FontName', 'Arial', 'FontSize',8, 'FontWeight', 'normal', 'horizontalAlignment', 'right');
                                end
                            end
                            %% ---- Source ripple line, subplot title ----
                            %                                 set(conTraceAx,'Color', sourceAreaColor) %set(subplot(2,2,1),'Color','Red')
                            iYlim(1,2) = min(introdeiripLFP); % - traceVertOffset(introde)); %last trace min
                            iripstarttime = 0; %iepochLFPtimes(ripStartIndices(irip));
                            %                                 iripendtime = iepochLFPtimes(ripEndIndices(irip));
                            Xline = [iripstarttime iripstarttime];
                            Yline = [[iYlim(1,2) iYlim(1,1)]];
                            line(Xline, Yline, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1.5);
                            
                            ititle = title(sprintf('%s %s Hz',F(ianimal).output{day}(epoch).LFPtypes{iLFPtype}, LFPrangesHz{iLFPtype}));
                            set(ititle,'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial');
                            set(gca,'children',flipud(get(gca,'children'))) %send the patch behind the LFP traces
                            ylim([iYlim(1,2) iYlim(1,1)])
                            %                         ylim auto
                            xlim([iripWinTimes(1) iripWinTimes(end)])
                            set(gca, 'YTick', []);
                            set(gca,'XTick',[iripWinTimes(1):(iripWinTimes(end) - iripWinTimes(1))/10:iripWinTimes(end)], 'FontSize',8,'FontWeight','normal')
                            set(gca, 'XTickLabel', round([iripWinTimes(1):(iripWinTimes(end) - iripWinTimes(1))/10:iripWinTimes(end)],2), 'XTickLabelRotation',30);
                            %                             set(gca, 'XTickLabel', [-win(1):win(1)/5:win(2)])
                            xlabel('time(s)','FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial')
                            %% ---- super title ----
                        end
                        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                        %                     sprtit = sprintf('%s D%d E%d %s %s', animalID, iday, iepoch, eventtype, envtypes{iepochEnv});
                        sprtit = sprintf('%s D%d %s %s STD>%d', animalID, day, eventtype, envtypes{iepochEnv}, minstdthresh);
                        iStitle = text(.5, .95, [{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]], 'Parent', sprtitleax, 'Units', 'normalized');
                        %                             iStitle = supertitle([{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]]);
                        set(iStitle, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
                        %% ---- pause, save figs ----
                        if pausefigs
                            pause
                        end
                        if savefigs
                            if ~isdir(figdirectory);
                                mkdir(figdirectory);
                            end
                            if ~isdir([figdirectory filenamesave '/area_mean']);
                                mkdir([figdirectory filenamesave '/area_mean']);
                            end
                            sprtitsave = strrep(sprtit,' ', '_'); %put '_' character in place of whitespace for filename
                            set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
                            currfigfile = sprintf('%s%s/area_mean/%s',figdirectory, filenamesave, sprtitsave);
                            print(currfigfile,'-dpng', '-r0')
                            disp(sprintf('plot %s saved', sprtit))
                        end
                        close all
                    end
                end
            end
        end
    end
    
    %% DEPRICATED Combine Ripple + EEG plots in same directory and create merged plots to facilitate viewing
    
    if 0 % createmergedplots && plotLFPtraces
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
    %% Dopopulation LFP fig, then population ripple timing fig
    
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