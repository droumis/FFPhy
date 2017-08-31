
%% To do
%make this compatibl with multi-antimal F output.. currently animals has to
%be a single animal
%% Dashboard
close all
runFilterFramework = 1;
; saveFilterOutput = runFilterFramework;
; loadFilterOutput = 0;
gatherResults = 1;
; saveResults = 0;
loadResults = 0;
plotfigs = 1;
; plotEnvironments = 1;
; plotStates = 0;
% ; plotSpec_EpochMean = 0;
; savefigs = 1;
; pausefigs = 0;
%% ---------------- plotting params --------------------------

%the total height and length of the figure dictates the sizing of the subplots
figLeft = 1;
figBottom = 1;
figHeight = 400;
figlength = 1800;
sfSpacing = .05;

sfStartLeft = sfSpacing;
sfEndRight = 1 - sfSpacing;
sfAllLength = sfEndRight - sfStartLeft;

sfAllBottom = .2;
sfAllTop = .8 - sfAllBottom;
clims = [-.25 3];
colorSet = 'DR1';
position = [.1 .1 .9 .8];
SpacingHorizontal = 0.01;
SpacingVertical = 0.02;
Spacing = 0.00;
Padding = 0.0;
MarginLeft = 0.04;
MarginRight = 0.04;
MarginTop = 0.09;
MarginBottom =  0.08;
usecolormap = 'jet';
%% ---------------- Data Filters --------------------------
filtfunction = 'riptrigspectrogram';
behavestruct = 'BehaveState';
animals = {'JZ1'};
days = [6];
eventtype = 'rippleskons';
eventSourceArea = 'ca1';
% epochEnvironment = 'openfield'; %wtrack, wtrackrotated, openfield, sleep
epochTypes = {'run'}; %'sleep', 'run', 
epochEnvironments = {'wtrack'};% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
eventSourceArea = 'ca1'; %ca1, mec, por
% ntAreas = {'ca1', 'sub', 'mec', 'por', 'v2l', 'ref'};
ntAreas = {'mec', 'por', 'v2l'};
fpass = [2 300];
timeHalfbwProduct = 2;% higher numbers increase time res, reduce freq res
numSlepianTapers = 3; % higher numbers reduce freq res (should be less than or equal to 2TW-1)
params = {};
params.Fs = 1500;
params.fpass = fpass;
params.trialave = 0;
params.tapers = [timeHalfbwProduct numSlepianTapers];
win = [.5 .5];
cwin = [0.1 0.01];

consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 3;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
minvelocity = 0;
maxvelocity = 4;
%% ---------------- Paths and Title strings ---------------------------------------------------
investInfo = animaldef(lower('Demetris'));
filtOutputDirectory = sprintf('%s%s/', investInfo{2}, filtfunction);
ResultsOutDirectory = sprintf('%s%s/', investInfo{3}, filtfunction);
figdirectory = sprintf('%s%s/', investInfo{4}, filtfunction);
filenamesave = sprintf('%s%sSTD%d_%s_%s_D%s_%dHB_%dTap', eventSourceArea, eventtype, minstdthresh, ...
    strjoin(epochEnvironments,'-'), cell2mat(animals),strjoin(arrayfun(@(x) num2str(x),days,'un',0),'-'),...
    timeHalfbwProduct, numSlepianTapers);
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
% filename = sprintf('riptrigspec_%s%s_%s_%sTets_%s.mat', eventSourceArea, eventtype, epochEnvironment, eventSourceArea, cell2mat(animals));
filenameTitle = strrep(filename,'_', '\_');
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1;
    eptypeEnv = [epochTypes; epochEnvironments];
    epochfilter = sprintf('((isequal($type, ''%s'')) && (isequal($environment, ''%s''))) || ',eptypeEnv{:});
    yuck = strfind(epochfilter, '||');
    epochfilter = epochfilter(1:yuck(end)-1);  %cut off the trailing '||'
    %     iterator = 'multitetrodeanal'; %multitetrodeanal
    %     epochfilter = sprintf('(isequal($type, ''%s'')) && (isequal($environment, ''%s''))',epochType, epochEnvironment); %'isequal($type, ''run'') && (isequal($environment, ''MultipleW''))'; %%'(isequal($type, ''sleep''))'; %%%&& isequal($descript, ''post-allruns''))';%   %%% %'isequal($type, ''run'') && isequal($environment, ''WTrackA'') && (($exposure>=1) && ($exposure<=10))';  %
    iterator = 'epocheeganal';
    %     tetfilter = sprintf('(isequal($area,''%s''))', eventSourceArea);
    tetfilter = sprintf('(isequal($area,''%s'')) || ', ntAreas{:});
    yuck = strfind(tetfilter, '||');
    tetfilter = tetfilter(1:yuck(end)-1); %cut off the trailing '||'
    timefilter{1} = {'getconstimes', '($cons == 1)',[eventSourceArea,eventtype],1,'consensus_numtets',consensus_numtets,...
        'minstdthresh',minstdthresh,'exclusion_dur',exclusion_dur,'minvelocity',minvelocity,'maxvelocity',maxvelocity};
    %----------F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);--------
    F = createfilter('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);
    %----------f = setfilteriterator(f, funcname, loadvariabl   es, options)--------
    F = setfilterfunction(F, ['dfa_' filtfunction], {'eeg', [eventSourceArea,eventtype]},'eventtype',eventtype, 'params', params, 'win', win, 'cwin', cwin);
    tic
    F = runfilter(F, 'useparpool', 1);
    F.filterTimer = toc;
    F.worldDateTime = clock;
    F.dataFilter = struct('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);
end
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    if ~isdir(filtOutputDirectory);
        mkdir(filtOutputDirectory);
    end
    save(sprintf('%s/%s',filtOutputDirectory, filename), 'F','-v7.3');
    disp(sprintf('filteroutput saved to %s/%s',filtOutputDirectory, filename))
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s/%s',filtOutputDirectory, filename))
    disp(sprintf('filteroutput loaded: %s/%s',filtOutputDirectory, filename))
end
%% ---------------- Collect across days and epochs for each site (ntrode) --------------------------------------------
if gatherResults == 1;
    for iAn = 1:length(F)
        animalinfo = animaldef(lower(animals{iAn}));
        animalID = animalinfo{1,3}; %use anim prefix for name
        FFanimdir =  sprintf('%s',animalinfo{1,2});
        matInds = cell2mat({F(iAn).output{1}.index}');
        [days, daysInds, daysInds2] = unique(matInds(:,[1]), 'rows', 'stable'); %get all matches across epochs and days
        days = days';
        dayepsenvset = cell(1,(length(epochEnvironments)));
        %% ---------- get dayepoch/environment mapping  --------------------------------
        for iday = 1:length(days(:,1))
            day = days(iday);
            load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'task', iday));
%             [eps, epsInds, epsInds2 ] = unique(matInds(matInds(:,1) == iday,[2]), 'rows', 'stable'); %get all matches across epochs and days
            eptypes = [];
            
            for iep = 1:length(task{iday}); %eps(:,1))
%                 ep = eps(iep);
                eptypes{iep,1} = task{iday}{iep}.environment;
            end
            [envtypes,~,ep2env] = unique(eptypes, 'stable');
            for ienv = 1:length(epochEnvironments)
                envienv = find(strcmp(epochEnvironments{ienv},envtypes));
                if ~envienv
                    continue
                end
                epsienv = find(ep2env == envienv);
                dayepsenvset{ienv} = [dayepsenvset{ienv}; [repmat(day,length(epsienv(:,1)),1) epsienv]];
            end
        end
        
        %% ---------- for each ntrode/environment for this animal, average all the event spectrograms -----------
        [nts, ntsInds, ntsInds2] = unique(matInds(:,[3]), 'rows', 'stable'); %get all matches across epochs and days
        Fg(iAn).epochEnvironments = epochEnvironments;
        Fg(iAn).t = [];
        Fg(iAn).f = [];
        for int = 1:length(nts(:,1))
            nt = nts(int);
            for ienv = 1:length(dayepsenvset)
                ienvNtInds = find(ismember(matInds(:,[1:2]),dayepsenvset{ienv}, 'rows') & matInds(:,3) == nt);
                Fg(iAn).meanspecs{ienv}{nt} = nanmean(cat(3,F(iAn).output{1}(ienvNtInds).S),3);
                if isempty(Fg(iAn).t)
                    Fg(iAn).ntIDs = nts;
                    Fg(iAn).t = F(iAn).output{1}(ienvNtInds).t;
                    Fg(iAn).f = F(iAn).output{1}(ienvNtInds).f;
                end
            end
        end
        
        %% ---------- Gather by wtrack -- inbound outbount correctoutbound mistakeoutbound--------------------------------
        wNTData = cell(max(nts),1);
        load(sprintf('%s%s%s.mat',animalinfo{1,2}, animalID, behavestruct));
        eventstate = [];
        eventstate.eventStartTimes = [];
        eventstate.state = [];
        eventstate.fields = [];
        wtrackInd = find(strcmp(epochEnvironments, 'wtrack'));
        [wdays, ~, wdaysInds2eps] = unique(dayepsenvset{wtrackInd}(:,1), 'stable');
        %get the animal state for every event
        for iday = 1:length(wdays(:,1))
            wday = wdays(iday);
            [weps, ~, wepsInds2dayepsmap] = unique(dayepsenvset{wtrackInd}(dayepsenvset{wtrackInd}(wdaysInds2eps(:,1)) == wday,2), 'rows', 'stable'); %get all matches across epochs and days
            %load linpos per day
            load(sprintf('%s%s%s%02d.mat',animalinfo{1,2}, animalID, 'linpos', wday));
            % get the data output indices of the [wday weps]
            % concat z-stack the wep events
            % get the eventstarttimes and the LFP-times
            % get a list of the eventtimes with event state
            % filter by eventstate types
            
            for iep = 1:length(weps(:,1))
                wep = weps(iep);
                wenvNtInds = find(ismember(matInds(:,[1:2]),[wday wep], 'rows'));
                load(sprintf('%sEEG/%s%s%02d-%d-01.mat',animalinfo{1,2}, animalID, 'eeg', wday, wep));
                numdata = length(eeg{wday}{wep}{1}.data);
                samprate = eeg{wday}{wep}{1}.samprate;
                starttime = eeg{wday}{wep}{1}.starttime;
                LFPtimes = [starttime:1/samprate:starttime+(numdata-1)/samprate]';
                eventStartTimes = F(iAn).output{1}(wenvNtInds(1)).triggers + starttime;
                %cat the event spects across days,eps for each ntrode
                wntrodes =  matInds(wenvNtInds,3);
                for iwntrode = 1:length(wntrodes)
                    wntrode = wntrodes(iwntrode);
                    wNTData{wntrode} = cat(3,wNTData{wntrode},F(iAn).output{1}(wenvNtInds(iwntrode)).S);
                end
                % get trajectory, segment, and linddist
                %                 eventStartIndices = F(iAn).output{wday}(wep).eventStartIndices;
                %                 LFPtimes = F(iAn).output{wday}(wep).LFPtimes;
                %                 eventStartTimes = LFPtimes(eventStartIndices);
                statematrix = linpos{wday}{wep}.statematrix;
                statemat.data = [statematrix.time statematrix.traj statematrix.segmentIndex statematrix.lindist];
                statemat.fields = 'postime traj segment lindist';
                stateindex = knnsearch(statemat.data(:,1), eventStartTimes);
                eventTrajs = statemat.data(stateindex,:);
                
                %get performance state
                trialIO = BehaveState.statechanges{wday}{wep}.statechangeseq;
                trialIOfields = BehaveState.statechanges{wday}{wep}.fields;
                corrcol = find(cell2mat(cellfun(@(x) strcmp(x,'correct'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
                timeportoutcol = find(cell2mat(cellfun(@(x) strcmp(x,'timeportout'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
                outBcol = find(cell2mat(cellfun(@(x) strcmp(x,'outbound'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
                inBcol = find(cell2mat(cellfun(@(x) strcmp(x,'inbound'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
                lastTimecol = find(cell2mat(cellfun(@(x) strcmp(x,'lasttime'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
                currTimecol = find(cell2mat(cellfun(@(x) strcmp(x,'currenttime'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
                outBCorr = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==1),[lastTimecol currTimecol corrcol timeportoutcol]);
                outBMist = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==0),[lastTimecol currTimecol corrcol timeportoutcol]);
                corrvec = list2vec(outBCorr(:,[1:2]),LFPtimes);
                mistvec = list2vec(outBMist(:,[1:2]),LFPtimes);
                
                eventindex = knnsearch(LFPtimes, eventStartTimes);
                tmpeventstate = [];
                tmpeventstate = [LFPtimes(eventindex) eventTrajs corrvec(eventindex) mistvec(eventindex)];
                eventstate.state = [eventstate.state; tmpeventstate];
                eventstate.fields = 'LFPtime postime traj segment lindist correct mistake';
                eventstate.eventStartTimes = [eventstate.eventStartTimes; eventStartTimes];
            end
        end
        if max(abs(diff(eventstate.state(:,1:2),[],2))) > 0.033;
            error('max event-time offset between pos and lfp times is more than 33ms (1 cam frame)')
        end
        if length(eventstate.state(:,1)) ~= size(wNTData{1},3)
            error('mismatch between number of state info and lfp sanples')
        end
        % split the events out into groups based on correct out / mistake out
        corrOutEventInd = find(eventstate.state(:,6) == 1);
        mistOutEventInd = find(eventstate.state(:,7) == 1);
        outBEventInd = [corrOutEventInd; mistOutEventInd];
        inBEventInd = setdiff([1:length(eventstate.state(:,1))], outBEventInd);
        wNTData = wNTData(~cell2mat(cellfun(@(x) isempty(x), wNTData, 'un', 0))); %get rid of empty cells caused by skipping ntrodes
        wDataType{1} = cellfun(@(x) nanmean(x(:,:,:),3), wNTData, 'un', 0);
        wDataType{2} = cellfun(@(x) nanmean(x(:,:,corrOutEventInd),3), wNTData, 'un', 0);
        wDataType{3} = cellfun(@(x) nanmean(x(:,:,mistOutEventInd),3), wNTData, 'un', 0);
        wDataType{4} = cellfun(@(x) nanmean(x(:,:,outBEventInd),3), wNTData, 'un', 0);
        wDataType{5} = cellfun(@(x) nanmean(x(:,:,inBEventInd),3), wNTData, 'un', 0);
        %         wDataTypeFields = {'all', 'corrOut', 'mistOut', 'outB', 'inB'};
        Fg(iAn).wStateData = wDataType;
        Fg(iAn).wStateDataFields = {'all', 'corrOut', 'mistOut', 'outB', 'inB'};
        %
        %     % add tags from tetinfo struct as a field in the F.tetout
        %     %load the tetinfostruct
        %     animalinfo = animaldef(animals{1});
        %     load([animalinfo{2},animalinfo{1} 'tetinfo']);
        %     supareainfo = cellfetch(tetinfo, 'suparea');
        %     areainfo = cellfetch(tetinfo, 'area');
        %     subareainfo = cellfetch(tetinfo, 'subarea');
        %
        %     %scavange the tags from each tetrode
        %     for ci = 1:length(F.tetout)
        %         tmpInds = F.tetout(ci).indices(1,:); %just grab the first Day Ep Tet ind... each tetrode should only be in one area for the duration of the exp
        %
        %         targetInd = find(ismember(supareainfo.index, tmpInds, 'rows'));
        %         F.tetout(ci).suptag = supareainfo.values{targetInd};
        %
        %         targetInd = find(ismember(areainfo.index, tmpInds, 'rows'));
        %         F.tetout(ci).areatag = areainfo.values{targetInd};
        %
        %         targetInd = find(ismember(subareainfo.index, tmpInds, 'rows'));
        %         F.tetout(ci).subtag = subareainfo.values{targetInd};
        %     end
    end
    clear F
    %% ---------------- Save Results Output ---------------------------------------------------
    if saveResults == 1;
        if ~isdir(ResultsOutDirectory);
            mkdir(ResultsOutDirectory);
        end
        save(sprintf('%s/%s',ResultsOutDirectory, filename), 'Fg','-v7.3');
        disp(sprintf('results output saved to %s/%s',ResultsOutDirectory, filename))
    end
end
%% ---------------- Load Results Output ---------------------------------------------------
if loadResults == 1;
    load(sprintf('%s/%s',ResultsOutDirectory, filename));
    disp(sprintf('loaded results output %s/%s',ResultsOutDirectory, filename))
end
%% plot
if plotfigs
    %% ---------------- Plot Averaged Across Epochs, Days for each ntrode/task -------------------------------------------
    if plotEnvironments == 1
        for iAn = 1:length(Fg)
            animalinfo = animaldef(lower(animals{iAn}));
            animalID = animalinfo{1,3}; %use anim prefix for name
            FFanimdir =  sprintf('%s',animalinfo{1,2});
            %% ---- loadtetinfostruct ----
            load([FFanimdir, animalID, 'tetinfo']);
            tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
            ntrodeIDs = Fg(iAn).ntIDs;
            ntets = size(unique(ntrodeIDs(:,1),'stable'),1);
            %         tets = unique(ntrodesIndices(:,1));
            %         if ntets ~= length(tets);
            %             error('data doesnt match indices')
            %         end
            %% ---- reorder the LFP traces by suparea|area|subarea tags (in that priority) ----
            [~, tagIndMap] = ismember(ntrodeIDs,tetinfoAll.index(:,3), 'rows');
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
            sfrows = floor(sqrt(ntets));
            sfcols = ceil(sqrt(ntets));
            
            for ienv = 1:length(Fg(iAn).epochEnvironments)
                if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
                    ifig = figure('Visible','off','units','normalized','position',position);
                else
                    ifig = figure('units','normalized','position',position);
                end
                set(gcf,'color','white')
                
                %% ---- loop across all tets for this day ----
                for introde = 1:ntets;
                    introdeID = ntrodeIDs(numsumSortInds(introde),1);
                    isupareatag = ntrodeTags{numsumSortInds(introde)}.suparea;
                    iareatag = ntrodeTags{numsumSortInds(introde)}.area;
                    isubareatag = ntrodeTags{numsumSortInds(introde)}.subarea;
                    iNTcolor = icolors(introde,:);
                    intfig = subaxis(sfrows,sfcols,introde, 'SpacingVertical', SpacingVertical, 'SpacingHorizontal', SpacingHorizontal, ...
                        'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight, 'MarginTop', MarginTop, ...
                        'MarginBottom', MarginBottom);
                    %% ---- Main Plotting ----
                    
                    %                 contourf(intfig,timeWin,frex,squeeze(ixpc.output{ianimal}{iDT}(numsumSortInds(introde),:,:))',num_frex,'linecolor','none');%, )
                    imagesc(Fg(iAn).t,Fg(iAn).f,Fg(iAn).meanspecs{ienv}{numsumSortInds(introde)}',[-0.5,.5])
                    %                 set(gca,'YDir','normal')
                    
%                     set(gca,'clim',clims,'ydir','normal','xlim',[-win(1) win(2)])
                    set(gca,'ydir','normal','xlim',[-win(1) win(2)])
                    hold on;
                    %                 set(gca, 'YScale', 'log')
                    
                    if mod(introde, sfcols) ~= 1;
                        set(gca, 'YTick', []);
                    else
                        set(gca, 'YTick',[round(linspace(min(Fg(iAn).f),max(Fg(iAn).f),8))],'FontSize',8, 'FontName', 'Arial');
                    end
                    if introde <= (sfrows-1)*sfcols;
                        set(gca, 'XTick', []);
                    else
                        set(gca, 'XTick',[win(1):win(2)/2:win(2)],'FontSize',8, 'FontName', 'Arial');
                    end
                    %% ---- Source ripple line, subplot title ----
                    Xline = [0 0];
                    Yline = [min(Fg(iAn).f) max(Fg(iAn).f)];
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
%                 caxis(clims);
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
                sprtit = sprintf('%s %s %s D%s %.1f-%.1fHz %dHB %dTap', filtfunction, epochEnvironments{ienv}, animalID, strrep(num2str(days), '  ', '-'), min(Fg(iAn).f),max(Fg(iAn).f), timeHalfbwProduct, numSlepianTapers);
                %             if plotNTrodesAcrossDays
                %                 sprTags = sprintf('%s %s', iarea, isubarea);
                %                 iclr = icolors(iIndInds(1),:);
                %                 iStitle = text(.5, .95, [{sprtit};...
                %                     {['\fontsize{12} \color[rgb]' sprintf('{%d %d %d} %s ', iclr, sprTags) '\color[rgb]{0 0 0} \fontsize{10}' sprintf('(phase %s %d)', phaseTetArea,selection_phasetet)]};...
                %                     {['\color[rgb]{.5 .5 .5} \fontsize{8}', sprintf(' {%s}', filenameTitle)]}], 'Parent', sprtitleax, 'Units', 'normalized');
                %                 set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
                %             else
                iStitle = text(.5, .95, [{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]], ...
                    'Parent', sprtitleax, 'Units', 'normalized');
                %             end
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
                clrbartit = text(posx(1)+posx(3)/2, posx(2)-MarginBottom/2, 'Z-power', 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial','horizontalAlignment', 'center');
                
                %% ---- pause, save figs ----
                %             return
                if pausefigs
                    pause
                end
                if savefigs
                    if ~isdir(figdirectory);
                        mkdir(figdirectory);
                    end
                    %                 if ~isdir([currfigdirectory filenamesave]);
                    %                     mkdir([currfigdirectory filenamesave]);
                    %                 end
                    sprtitsave = strrep(sprtit,' ', '_'); %put '_' character in place of whitespace for filename
                    set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
                    currfigfile = sprintf('%s%s',figdirectory, sprtitsave);
                    print(currfigfile,'-dpng', '-r0')
                    disp(sprintf('plot %s saved', sprtit))
                    
                    %             imagesc(F.tetout(sortInds(t)).t,F.tetout(sortInds(t)).f,F.meanspecs{sortInds(t)}',[-0.2,2.5])
                    %             set(gca,'YDir','normal')
                    %             %         plot(rand(3,4))
                    %             %         axis xy
                    %             %             xlabel('seconds from rip start','FontSize',12,'FontWeight','bold','Color','k')
                    %                 set(gca, 'YTick', []);
                    %                 set(gca,'XTick', [])
                    %                 set(gca, 'XTickLabel', [])
                end
                
                close all
                
            end
        end
    end
        %% ---------------- Plot w track conditions -------------------------------------------
    if plotStates == 1
        for iAn = 1:length(Fg)
            animalinfo = animaldef(lower(animals{iAn}));
            animalID = animalinfo{1,3}; %use anim prefix for name
            FFanimdir =  sprintf('%s',animalinfo{1,2});
            %% ---- loadtetinfostruct ----
            load([FFanimdir, animalID, 'tetinfo']);
            tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
            ntrodeIDs = Fg(iAn).ntIDs;
            ntets = size(unique(ntrodeIDs(:,1),'stable'),1);
            %         tets = unique(ntrodesIndices(:,1));
            %         if ntets ~= length(tets);
            %             error('data doesnt match indices')
            %         end
            %% ---- reorder the LFP traces by suparea|area|subarea tags (in that priority) ----
            [~, tagIndMap] = ismember(ntrodeIDs,tetinfoAll.index(:,3), 'rows');
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
            sfrows = floor(sqrt(ntets));
            sfcols = ceil(sqrt(ntets));
            
            for istate = 1:length(Fg(iAn).wStateData)
                if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
                    ifig = figure('Visible','off','units','normalized','position',position);
                else
                    ifig = figure('units','normalized','position',position);
                end
                set(gcf,'color','white')
                
                %% ---- loop across all tets for this day ----
                for introde = 1:ntets;
                    introdeID = ntrodeIDs(numsumSortInds(introde),1);
                    isupareatag = ntrodeTags{numsumSortInds(introde)}.suparea;
                    iareatag = ntrodeTags{numsumSortInds(introde)}.area;
                    isubareatag = ntrodeTags{numsumSortInds(introde)}.subarea;
                    iNTcolor = icolors(introde,:);
                    intfig = subaxis(sfrows,sfcols,introde, 'SpacingVertical', SpacingVertical, 'SpacingHorizontal', SpacingHorizontal, ...
                        'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight, 'MarginTop', MarginTop, ...
                        'MarginBottom', MarginBottom);
                    %% ---- Main Plotting ----
                    
                    %                 contourf(intfig,timeWin,frex,squeeze(ixpc.output{ianimal}{iDT}(numsumSortInds(introde),:,:))',num_frex,'linecolor','none');%, )
                    if istate == 6;
                        spect2plot = (Fg(iAn).wStateData{4}{numsumSortInds(introde)}' - Fg(iAn).wStateData{5}{numsumSortInds(introde)}');
                    else
                        spect2plot = Fg(iAn).wStateData{istate}{numsumSortInds(introde)}';
                    end
                    imagesc(Fg(iAn).t,Fg(iAn).f,spect2plot,[-0.5,.5])
                    %                 set(gca,'YDir','normal')
                    
                    set(gca,'clim',clims,'ydir','normal','xlim',[-.5 .5])
                    hold on;
                    %                 set(gca, 'YScale', 'log')
                    
                    if mod(introde, sfcols) ~= 1;
                        set(gca, 'YTick', []);
                    else
                        set(gca, 'YTick',[round(linspace(min(Fg(iAn).f),max(Fg(iAn).f),8))],'FontSize',8, 'FontName', 'Arial');
                    end
                    if introde <= (sfrows-1)*sfcols;
                        set(gca, 'XTick', []);
                    else
                        set(gca, 'XTick',[win(1):win(2)/2:win(2)],'FontSize',8, 'FontName', 'Arial');
                    end
                    %% ---- Source ripple line, subplot title ----
                    Xline = [0 0];
                    Yline = [min(Fg(iAn).f) max(Fg(iAn).f)];
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
                if istate == 6
                    sprtit = sprintf('%s %s %s D%s %.1f-%.1fHz', filtfunction, 'outB-inB', animalID, strrep(num2str(days), '  ', '-'), min(Fg(iAn).f),max(Fg(iAn).f));
                else
                    sprtit = sprintf('%s %s %s D%s %.1f-%.1fHz', filtfunction, Fg(iAn).wStateDataFields{istate}, animalID, strrep(num2str(days), '  ', '-'), min(Fg(iAn).f),max(Fg(iAn).f));
                end
                %             if plotNTrodesAcrossDays
                %                 sprTags = sprintf('%s %s', iarea, isubarea);
                %                 iclr = icolors(iIndInds(1),:);
                %                 iStitle = text(.5, .95, [{sprtit};...
                %                     {['\fontsize{12} \color[rgb]' sprintf('{%d %d %d} %s ', iclr, sprTags) '\color[rgb]{0 0 0} \fontsize{10}' sprintf('(phase %s %d)', phaseTetArea,selection_phasetet)]};...
                %                     {['\color[rgb]{.5 .5 .5} \fontsize{8}', sprintf(' {%s}', filenameTitle)]}], 'Parent', sprtitleax, 'Units', 'normalized');
                %                 set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
                %             else
                iStitle = text(.5, .95, [{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]], ...
                    'Parent', sprtitleax, 'Units', 'normalized');
                %             end
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
                clrbartit = text(posx(1)+posx(3)/2, posx(2)-MarginBottom/2, 'Z-power', 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial','horizontalAlignment', 'center');
                
                %% ---- pause, save figs ----
                %             return
                if pausefigs
                    pause
                end
                if savefigs
                    if ~isdir(figdirectory);
                        mkdir(figdirectory);
                    end
                    %                 if ~isdir([currfigdirectory filenamesave]);
                    %                     mkdir([currfigdirectory filenamesave]);
                    %                 end
                    sprtitsave = strrep(sprtit,' ', '_'); %put '_' character in place of whitespace for filename
                    set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
                    currfigfile = sprintf('%s%s',figdirectory, sprtitsave);
                    print(currfigfile,'-dpng', '-r0')
                    disp(sprintf('plot %s saved', sprtit))
                    
                    %             imagesc(F.tetout(sortInds(t)).t,F.tetout(sortInds(t)).f,F.meanspecs{sortInds(t)}',[-0.2,2.5])
                    %             set(gca,'YDir','normal')
                    %             %         plot(rand(3,4))
                    %             %         axis xy
                    %             %             xlabel('seconds from rip start','FontSize',12,'FontWeight','bold','Color','k')
                    %                 set(gca, 'YTick', []);
                    %                 set(gca,'XTick', [])
                    %                 set(gca, 'XTickLabel', [])
                end
                
                close all
                
            end
        end
    end
end