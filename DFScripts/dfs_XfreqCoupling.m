

%to do: run other areas aside from ca1
% parrelelioze by tetrode in the iterator.. maybe use multieeg instead and
% parallelize in the filtfunction
%add the subplot code from the spect script to this one
%add duration, state
%seperate by epoch to double check consistency across a day

% clear all
% close all
runFilterFramework = 0;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 0;
gatherXFC = 0;
resaveEpochsMean = 0;
plotXFC = 1;
%% ---------------- plotting params --------------------------

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
selection_phasetet = 0;%  0: Phase tet is same tetrode %  1, 2, 3, 4 : other options refer to STAreftet (cc, CA2, CA3, DG)
norip =  {'getconstimes', '($cons == 0)', [eventSourceArea eventtype],1,'consensus_numtets',consensus_numtets,...
    'minstdthresh',minstdthresh,'exclusion_dur',exclusion_dur,'maxvelocity',maxvelocity,'minvelocity',minvelocity};
vel4 = {'kk_get2dstate', '(abs($velocity) > 4)'};
vel20 = {'kk_get2dstate', '(abs($velocity) > 20)'};
velunder20 = {'kk_get2dstate', '((abs($velocity) < 20))'};
immobile2 = {'kk_get2dstate', '($immobility == 1)','immobility_velocity',4,'immobility_buffer',1}; % notice 1 s buffer
clear statespec
if 1
    statespec_descript = 'state1: vel4 // state2: immobile2   ';
    %     statespec{1} = { vel20 };
    %     statespec{2} = { vel4, velunder20, norip };
    statespec{1} = { vel4, norip };
    statespec{2} = { immobile2 } ;
elseif 1
    statespec_descript = 'state1: sleepc NREM // state2: sleepc, REM  ';
    statespec{1} = { sleepc, noREM } ;
    statespec{2} = { sleepc, REM } ;
end
ntAreas = {'ca1', 'sub', 'mec', 'por', 'v2l', 'ref'};
animals = {'JZ1'};
days = [1];
filtfunction = 'XfreqCoupling';
epochType = 'run'; %run or sleep
epochEnvironment = 'wtrack';% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
investInfo = animaldef(lower('Demetris'));
%% ---------------- Paths and Title strings ---------------------------------------------------
filenamesave = sprintf('%s', epochEnvironment); %add params here
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1;
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
        iF.datafilter = whos; %save all filtering parameters in workspace into struct
        animal = animals{iAn}; %loop over animals to get access to the currnt animal in the iterator in order to run multiple internal timefillters via statespec
        %----------F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);--------
        iF = createfilter('animal', animal, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,...
            'iterator', iterator);
        %----------f = setfilteriterator(f, funcname, loadvariables, options)--------
        iF = setfilterfunction(iF,'dfa_XfreqCoupling',{'eeg'},'animal', animal, 'statespec', statespec,'statespec_descript', statespec_descript,...
            'PhaseFreqVector',PhaseFreqVector,'AmpFreqVector',AmpFreqVector,...
            'PhaseFreq_BandWidth',PhaseFreq_BandWidth,'AmpFreq_BandWidth',AmpFreq_BandWidth, 'selection_phasetet',selection_phasetet);
        tic
        iF = runfilter(iF);
        iF.filterTimer = toc; iF.filterTimer
        iF.worldDateTime = clock;
        F{iAn} = iF;
    end
end
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    if ~isdir(sprintf('%s%s', investInfo{2}, filtfunction));
        mkdir(sprintf('%s%s', investInfo{2}, filtfunction));
    end
    save(sprintf('%s%s/%s',investInfo{2}, filtfunction, filename), 'F','-v7.3');
    disp(sprintf('filteroutput saved to %s%s/%s',investInfo{2}, filtfunction, filename))
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s%s/%s',investInfo{2},filtfunction, filename))
end
%% ----------------  Gather ntrode data across epochs ---------------------------------------------------
if gatherXFC
    %Collect data on individual tetrodes from across epochs
    for iAn = 1:length(animals)
        indices = cell2mat(arrayfun(@(x) x.index, F{1}.output{1}, 'UniformOutput', false)'); % collect the ntrode's epochs
        [tets ia ic] = unique(indices(:,3), 'stable'); %find the indices of the epochs for each ntrode
        statespecind = find(cell2mat(arrayfun(@(x) strcmp(x.name, 'statespec'), F{1}.datafilter, 'UniformOutput', false))); %find the ind of statespec
        numstates = F{iAn}.datafilter(statespecind).size(2);
        F{iAn}.tetout = struct;
        for istate = 1:numstates
            for iuniq = 1:length(tets)
                xEps = find(ic == iuniq); %across epoch matches for day, ntrode
                F{iAn}.tetout(iuniq).XFC{istate}=[];
                F{iAn}.tetout(iuniq).indices = [];
                for iEp = 1:length(xEps) % concat in Z dim
                    F{iAn}.tetout(iuniq).XFC{istate} = cat(3, F{iAn}.tetout(iuniq).XFC{istate}, F{iAn}.output{1}(xEps(iEp)).XFC{istate});
                    F{iAn}.tetout(iuniq).indices = cat(1,F{iAn}.tetout(iuniq).indices, F{iAn}.output{1}(xEps(iEp)).index);
                end
                F{iAn}.tetout(iuniq).XFCcomb{istate} = zeros(size(F{iAn}.tetout(iuniq).XFC{istate},1), size(F{iAn}.tetout(iuniq).XFC{istate},2)); %initialize
                F{iAn}.tetout(iuniq).XFCcomb{istate} = nanmean(F{iAn}.tetout(iuniq).XFC{istate},3); %mean across epochs
                F{iAn}.tetout(iuniq).frequency_phase = F{iAn}.output{1}(xEps(iEp)).frequency_phase;
                F{iAn}.tetout(iuniq).phasefreq_bandwidth = F{iAn}.output{1}(xEps(iEp)).phasefreq_bandwidth;
                F{iAn}.tetout(iuniq).frequency_amplitude = F{iAn}.output{1}(xEps(iEp)).frequency_amplitude;
                F{iAn}.tetout(iuniq).ampfreq_bandwidth = F{iAn}.output{1}(xEps(iEp)).ampfreq_bandwidth;
            end
        end
    end
end
%% ----------------  resave gathered data across epochs ---------------------------------------------------
if resaveEpochsMean == 1
    save(sprintf('%s%s/%s',investInfo{2}, filtfunction, filename), 'F','-v7.3');
    disp(sprintf('F saved to %s%s/%s',investInfo{2}, filtfunction, filename))
end

%% ----------------   Plot results by tetrode across all days and epochs. ----------------  
if plotXFC
    for iAn = 1:length(animals)
        indices = cell2mat(arrayfun(@(x) x.index, F{1}.output{1}, 'UniformOutput', false)'); % collect the ntrode's epochs
        [tets ia ic] = unique(indices(:,3), 'stable'); %find the indices of the epochs for each ntrode
        statespecind = find(cell2mat(arrayfun(@(x) strcmp(x.name, 'statespec'), F{1}.datafilter, 'UniformOutput', false))); %find the ind of statespec
        numstates = F{iAn}.datafilter(statespecind).size(2);
        
        for istate = 1:numstates
            for itet = 1:length(tets)
                ixfc = F{iAn}.tetout(itet).XFCcomb{istate};
                iVecPhase = F{iAn}.tetout(itet).frequency_phase;
                iVecAmp = F{iAn}.tetout(itet).frequency_amplitude;
                iPhasebin = F{iAn}.tetout(itet).phasefreq_bandwidth;
                iAmpbin = F{iAn}.tetout(itet).ampfreq_bandwidth;
                contourf(iVecPhase+iPhasebin/2, iVecAmp+iAmpbin/2, ixfc', 30, 'lines', 'none')
%                 caxis([0 .001])
                set(gca,'fontsize',14)
                ylabel('Amplitude Frequency (Hz)')
                xlabel('Phase Frequency (Hz)')
                colormap('jet')
%                 colorbar
            end
        end
    end
end
