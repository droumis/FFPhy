%to do:

%run other areas aside from ca1.. use tet 28 as the phasereftet

%seperate by epoch to double check consistency across a day

%plot all tets, states, days on same figure

%run sleep, openfield


% clear all
close all
runFilterFramework =0;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 1;
gatherXepochs = 1;
resaveEpochsMean = gatherXepochs;
plotXFC = 1;
savefigs = 1;
pausefigs = 0;
%% ---------------- plotting params --------------------------
colorSet = 'DR1';
position = [.1 .1 .8 .6];
SpacingVert = 0.00;
Spacing = 0.00;
Padding = 0.0;
MarginLeft = 0.06;
MarginRight = 0.04;
MarginTop = 0.12;
MarginBottom =  0.07;
usecolormap = 'jet';
%% ---------------- Data Filters --------------------------
consensus_numtets = 3;   % minimum # of tets for consensus event detection
minstdthresh = 2;        % STD. how big your ripples are
exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
minvelocity = 0;
maxvelocity = 4;
eventtype = 'rippleskons';
eventSourceArea = 'ca1';

% target_det =   [];
% stareftet_only = 0;
PhaseFreqVector=0:0.5:50;
AmpFreqVector=10:5:250;
PhaseFreq_BandWidth=2;
AmpFreq_BandWidth=5;
selection_phasetet = 20;%  0: Phase tet is same tetrode %  1, 2, 3, 4 : other options refer to STAreftet (cc, CA2, CA3, DG)
norip =  {'getconstimes', '($cons == 0)', [eventSourceArea eventtype],1,'consensus_numtets',consensus_numtets,...
    'minstdthresh',minstdthresh,'exclusion_dur',exclusion_dur,'maxvelocity',maxvelocity,'minvelocity',minvelocity};
rip =  {'getconstimes', '($cons == 1)', [eventSourceArea eventtype],1,'consensus_numtets',consensus_numtets,...
    'minstdthresh',minstdthresh,'exclusion_dur',exclusion_dur,'maxvelocity',maxvelocity,'minvelocity',minvelocity};
vel4 = {'kk_get2dstate', '(abs($velocity) > 4)'};
vel20 = {'kk_get2dstate', '(abs($velocity) > 20)'};
velunder20 = {'kk_get2dstate', '((abs($velocity) < 20))'};
immobile2 = {'kk_get2dstate', '($immobility == 1)','immobility_velocity',4,'immobility_buffer',1}; % notice 1 s buffer
clear statespec
%open timefilterscript.m for options
% statespec_descript = 'state1: vel4 // state2: immobile2 '; %format: state<n>: <statenlabel> // ...
% statespec{1} = { vel4, norip };
% statespec{2} = { immobile2 } ;

% statespec_descript = 'state1: rip  // state2: immobile2 '; %format: state<n>: <statenlabel> // ...
% statespec{1} = { immobile2, rip};
% statespec{2} = { immobile2 } ;

statespec_descript = 'state1: rip  // state2: norip '; %format: state<n>: <statenlabel> // ...
statespec{1} = { immobile2, rip};
statespec{2} = { immobile2, norip};


%     statespec{1} = { vel20 };
%     statespec{2} = { vel4, velunder20, norip };
%     statespec_descript = 'state1: sleepc NREM // state2: sleepc, REM  ';
%     statespec{1} = { sleepc, noREM } ;
%     statespec{2} = { sleepc, REM } ;
% end
ntAreas = {'ca1', 'sub', 'mec', 'por', 'v2l', 'ref'};
% ntAreas = {'v2l'};
% ntAreas = {'mec'};
animals = {'JZ1'};
days = [1:6];
filtfunction = 'XfreqCoupling';
epochType = 'run'; %run or sleep
epochEnvironment = 'wtrack';% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
investInfo = animaldef(lower('Demetris'));

%% ---------------- Paths and Title strings ---------------------------------------------------
statessave = strsplit(statespec_descript, '//'); %split up state labels
statessave = cellfun(@(x) strtrim(x(strfind(x, ':')+1:end)),statessave, 'UniformOutput', false); %grabstate label and trim whitespace
filtOutputDirectory = sprintf('%s%s/', investInfo{2}, filtfunction);
figdirectory = sprintf('%s%s/', investInfo{4}, filtfunction);
filenamesave = sprintf('%s_%s_PhaseTet%d_%s_D%s', epochEnvironment, strjoin(statessave,'-'), selection_phasetet, strjoin(animals,'-'), strjoin(arrayfun(@(x) num2str(x),days,'UniformOutput',false),'-'));
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
filenameTitle = strrep(filename,'_', '\_');
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1;
    %     parpool(24); %use 96 workers on typhoon. use 6 on haboob
    %     iterator = 'eeganal_stableOutputStruct';
    iterator = 'eeganal';
    %tetfilter: the ntrodeID's that pass this filter get stored into f.eegdata{day}{epoch}[ntrodeIDs]
    %         tetfilter = '(isequal($area, ''ca1''))';
    tetfilter = sprintf('(isequal($area,''%s'')) || ', ntAreas{:});
    yuck = strfind(tetfilter, '||');
    tetfilter = tetfilter(1:yuck(end)-1); %cut off the trailing '||'
    timefilter = {};  % done via statespec, within analysis function
    epochfilter =    sprintf('(isequal($type, ''%s'')) && (isequal($environment, ''%s''))',epochType, epochEnvironment);
    for iAn = 1:length(animals)
        clear iF
        animal = animals{iAn}; %loop over animals to get access to the currnt animal in the iterator in order to run multiple internal timefillters via statespec
        %----------F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);--------
        iF = createfilter('animal', animal, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,...
            'iterator', iterator);
        %----------f = setfilteriterator(f, funcname, loadvariables, options)--------
        iF = setfilterfunction(iF,'dfa_XfreqCoupling',{'eeg'},'animal', animal, 'statespec', statespec,'statespec_descript', statespec_descript,...
            'PhaseFreqVector',PhaseFreqVector,'AmpFreqVector',AmpFreqVector,...
            'PhaseFreq_BandWidth',PhaseFreq_BandWidth,'AmpFreq_BandWidth',AmpFreq_BandWidth, 'selection_phasetet',selection_phasetet);
        iF.datafilter = whos; %save info about filtering parameters in workspace into struct
        save(sprintf('%s/_filterparams_%s',filtOutputDirectory, filename)); %save filtering params to a '.filterparams' file
        tic
        iF = runfilter(iF, 'outputDayEpTetCells', 1); %parpool = 0, outputDayEpTetCells = 0 by default
        iF.filterTimer = toc; iF.filterTimer
        iF.worldDateTime = clock;
        F{iAn} = iF;
    end
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
end
%% ----------------  Gather ntrode data across epochs ---------------------------------------------------
if gatherXepochs
    %Collect data on individual tetrodes from across epochs
    for iAn = 1:length(animals)
        F{iAn}.tetout = [];
        disp(['@@@@@@@@@@@ Animal ' animals{iAn}])
        days = find(~cellfun(@isempty,F{iAn}.output));
        for iday = 1:length(days)
            day = days(iday);
            indices = cellfetch(F{iAn}.output, 'index');
            indices = indices.index;
            %             indices = cell2mat(arrayfun(@(x) x.index, F{iAn}.output{iday}, 'UniformOutput', false)'); % collect the ntrode's epochs
            [tets ia ic] = unique(indices(:,3), 'stable'); %find the indices of the epochs for each ntrode
            statespecind = find(cell2mat(arrayfun(@(x) strcmp(x.name, 'statespec'), F{iAn}.datafilter, 'UniformOutput', false))); %find the ind of statespec
            numstates = F{iAn}.datafilter(statespecind).size(2);
            %             F{iAn}.tetout{iday} = struct;
            eps = find(~cellfun(@isempty,F{iAn}.output{day}));
            for istate = 1:numstates
                disp(['$$$$$ state ' num2str(istate)])
                tets = find(~cellfun(@isempty,F{iAn}.output{day}{eps(1)}));
                for itet = 1:length(tets)
                    tet = tets(itet);
                    disp(['ntrode' num2str(tet)])
                    %                     xEps = find(ic == iuniq); %across epoch matches for day, ntrode
                    F{iAn}.tetout{day}{tet}.XFC{istate}=[];
                    F{iAn}.tetout{day}{tet}.indices = [];
                    F{iAn}.tetout{day}{tet}.duration = [];
                    useEps = find(~cellfun(@isempty,F{iAn}.output{day}));
                    for iEp = 1:length(useEps) % concat in Z dim
                        ep = useEps(iEp);
                        F{iAn}.tetout{day}{tet}.XFC{istate} = cat(3, F{iAn}.tetout{day}{tet}.XFC{istate}, F{iAn}.output{day}{ep}{tet}.XFC{istate});
                        F{iAn}.tetout{day}{tet}.indices = cat(1,F{iAn}.tetout{day}{tet}.indices, F{iAn}.output{day}{ep}{tet}.index);
                        F{iAn}.tetout{day}{tet}.duration = cat(1,F{iAn}.tetout{day}{tet}.duration, F{iAn}.output{day}{ep}{tet}.duration);
                    end
                    F{iAn}.tetout{day}{tet}.XFCcomb{istate} = zeros(size(F{iAn}.tetout{day}{tet}.XFC{istate},1), size(F{iAn}.tetout{day}{tet}.XFC{istate},2)); %initialize
                    F{iAn}.tetout{day}{tet}.XFCcomb{istate} = nanmean(F{iAn}.tetout{day}{tet}.XFC{istate},3); %mean across epochs
                    F{iAn}.tetout{day}{tet}.frequency_phase = F{iAn}.output{day}{ep}{tet}.frequency_phase;
                    F{iAn}.tetout{day}{tet}.phasefreq_bandwidth = F{iAn}.output{day}{ep}{tet}.phasefreq_bandwidth;
                    F{iAn}.tetout{day}{tet}.frequency_amplitude = F{iAn}.output{day}{ep}{tet}.frequency_amplitude;
                    F{iAn}.tetout{day}{tet}.ampfreq_bandwidth = F{iAn}.output{day}{ep}{tet}.ampfreq_bandwidth;
                    F{iAn}.tetout{day}{tet}.statespec_descript = F{iAn}.output{day}{ep}{tet}.statespec_descript;
                end
            end
        end
    end
end
%% ----------------  resave gathered data across epochs ---------------------------------------------------
if resaveEpochsMean == 1
    save(sprintf('%s/%s',filtOutputDirectory, filename), 'F','-v7.3');
    disp(sprintf('F saved to %s/%s',filtOutputDirectory, filename))
end

%% ----------------   Plot results by tetrode across all days and epochs. ----------------
if plotXFC
    for iAn = 1:length(animals)
        %% ---------- Get animal info--------------------------------
        ianimalinfo = animaldef(lower(animals{iAn}));
        animalID = ianimalinfo{1,3}; %use anim prefix for name
        FFanimdir =  sprintf('%s',ianimalinfo{1,2});
        load([FFanimdir, animalID, 'tetinfo']);
        tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
        %         indices = cell2mat(arrayfun(@(x) x.index, F{iAn}.output{1}, 'UniformOutput', false)'); % collect the ntrode's epochs
        indices = cellfetch(F{iAn}.output, 'index');
        indices = indices.index;
        [tets ia ic] = unique(indices(:,3), 'stable'); %find the indices of the epochs for each ntrode
        %         [DETindices ia ic] = unique(indices,'rows', 'stable'); %find the indices of the epochs for each ntrode
        statespecind = find(cell2mat(arrayfun(@(x) strcmp(x.name, 'statespec'), F{iAn}.datafilter, 'UniformOutput', false))); %find the ind of statespec
        numstates = F{iAn}.datafilter(statespecind).size(2);
        
        %% ---------- Get ntrode tags --------------------------------
        [~, tagIndMap] = ismember(indices,tetinfoAll.index, 'rows');
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
        %         numsumallareatags = cell2mat([numsumSupAreas numsumAreas numsumSubAreas]);
        %         [numsumallSort, numsumSortInds] = sortrows(numsumallareatags);%,[-1 -2 -3]); % -Col to 'descend'
        %         icolors = icolors(numsumSortInds,:);
        days = find(~cellfun(@isempty,F{iAn}.tetout));
        tets = find(~cellfun(@isempty,F{iAn}.tetout{days(1)})); %look into first day for the num of tetrodes
        for itet = 1:length(tets)   
            tet = tets(itet);
            if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
                ifig = figure('Visible','off','units','normalized','position',position);
            else
                ifig = figure('units','normalized','position',position);
            end
            set(gcf,'color','white')
            for iday = 1:length(days)
                day = days(iday);
                idayEpTet = F{iAn}.tetout{day}{tet}.indices;
                [~, iIndInds] = ismember(idayEpTet,indices,'rows');
                iarea = char(cell2mat(strAreas(iIndInds(1)))); %grab the area tags
                isubarea = num2str(cell2mat(strSubAreas(iIndInds(1))));
                idurs = sum(F{iAn}.tetout{day}{tet}.duration,2);
                subtitles = strsplit(F{iAn}.tetout{day}{tet}.statespec_descript, '//');
                for istate = 1:numstates
                    subaxis(numstates,length(days),iday+length(days)*(istate-1), 'Spacing', Spacing, 'SpacingVert', SpacingVert, 'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight,...
                        'MarginTop', MarginTop, 'MarginBottom', MarginBottom);
                    ixfc = F{iAn}.tetout{day}{tet}.XFCcomb{istate};
                    iVecPhase = F{iAn}.tetout{day}{tet}.frequency_phase;
                    iVecAmp = F{iAn}.tetout{day}{tet}.frequency_amplitude;
                    iPhasebin = F{iAn}.tetout{day}{tet}.phasefreq_bandwidth;
                    iAmpbin = F{iAn}.tetout{day}{tet}.ampfreq_bandwidth;
                    contourf(iVecPhase+iPhasebin/2, iVecAmp+iAmpbin/2, ixfc', 30, 'lines', 'none')
                    if iday == 1
                        yt = get(gca,'YTick');
                        set(gca,'YTick',yt,'FontName','Arial','fontsize',10)
                        ylab = ylabel(sprintf('%s',subtitles{istate}));
                        set(ylab, 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial')
                    else
                        set(gca, 'YTick', []);
                    end
                    if istate == numstates;
                        xt = get(gca,'XTick');
                        set(gca,'XTick',xt,'FontName','Arial','fontsize',10)
                    else
                        isubtit = title(sprintf('Day %d',days(iday)));
                        set(isubtit,'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial')
                        set(gca, 'XTick', []);
                    end
                    idur = sprintf('%.02f',idurs(istate)/60);
                    text(.98, .96, idur, 'FontSize',7,'FontWeight','bold','Color',[.5 .5 .5], 'FontName', 'Arial', 'Units', 'normalized', 'horizontalAlignment', 'right');
                end
%                 xtic = get(gca,'XTick');
            end
            %% ---- super title and colorbar----
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
%             posfig=get(gcf,'position');
            xlab = 'Phase Frequency (Hz)';
            ylab = 'Amplitude Frequency (Hz)';
            supxlabel = text(.5, .02, xlab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
            supylabel = text(.01, .5, ylab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'rotation', 90, 'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
            colormap(usecolormap)
            clrbar = colorbar('location','eastoutside', 'FontSize',6,'FontName', 'Arial');%, 'FontWeight','bold');
            %             ylabel(clrbar,'comodulation', 'FontSize', 10, 'FontWeight','bold')
            posx1=get(gca,'position');
            posx=get(clrbar,'Position');
            posx(1)= 1-MarginRight+.01;
            posx(2)= MarginBottom;
            posx(3)= .01;
            posx(4)= 1 - MarginTop - MarginBottom;
            set(clrbar,'Position',posx)
            set(gca,'position',posx1)
            %             clrbar.Label.String = 'comodulation';
            %             clrbar.AxisLocation = 'in';
            if selection_phasetet == 0
                phaseTetArea = 'self';
            else
                [~, phaseTetInd] = ismember(selection_phasetet,indices(:,3));
                phaseTetArea = strAreas{phaseTetInd};
            end
            sprtit = sprintf('%s D%s E%s T%d', ianimalinfo{1}, strjoin(arrayfun(@(x) num2str(x),days','UniformOutput',false),'-'), strjoin(arrayfun(@(x) num2str(x),idayEpTet(:,2)','UniformOutput',false),'-'), idayEpTet(1,3));
            sprTags = sprintf('%s %s', iarea, isubarea);
            %             yuck = strfind(sprtit, ';');
            %             sprtit = sprtit(1:yuck(end)-1); %cut off the trailing ';'
            iclr = icolors(iIndInds(1),:);
            iStitle = text(.5, .95, [{sprtit};...
                {['\fontsize{12} \color[rgb]' sprintf('{%d %d %d} %s ', iclr, sprTags) '\color[rgb]{0 0 0} \fontsize{10}' sprintf('(phase %s %d)', phaseTetArea,selection_phasetet)]};...
                {['\color[rgb]{.5 .5 .5} \fontsize{8}', sprintf(' {%s}', filenameTitle)]}], 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle, 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
            clrbartit = text(posx(1)+posx(3)/2, posx(2)-MarginBottom/2, 'comod', 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial','horizontalAlignment', 'center');
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
                savetit = [sprintf('%s ph%d', sprtit, selection_phasetet)]; %save with phase tet label
                sprtitsave = strrep(savetit,' ', '_'); %put '_' character in place of whitespace for filename
                set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
                currfigfile = sprintf('%s%s/%s', figdirectory, filenamesave, sprtitsave);
                print(currfigfile,'-dpng', '-r0')
                disp(sprintf('plot %s saved', savetit))
            end
            close all
            
        end
    end
end
