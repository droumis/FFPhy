


%% Dashboard
%%% First create the processed data structure similar to dfs_riptriglfp.m
close all
runFilterFramework = 1;
    ; saveFilterOutput = runFilterFramework;
    ; loadFilterOutput = 0;

%%% Then run phasecoherence on the LFP traces
calcEventState = 0;
saveEventState = 0;%calcEventState;
calculateITPC = 0;
calculateISPC = 0;
calculateIXPC = calculateITPC || calculateISPC;
loadFilterOutput = calculateIXPC || calcEventState;
saveResultsOutput = calculateIXPC;
runPermTest = 0; %calculateIXPC;
loadResultsOutput = 0;
plotITPC = 0;
plotISPC = 0;
savefigs = 0;
pausefigs = 0;
calcfunction = 'ITPC'; %ITPC
behavestruct = 'BehaveState';
plotoutputtype = 'phase'; % 'power' or 'phase'
%% ---------------- plotting params --------------------------
colorSet = 'DR1';
% if strcmp(plotoutputtype, 'phase')
%     clims = [0 .7]; % mean ITPC goes from 0-1
%     usecolormap = 'jet'; %cbrewer('seq', 'Greens', 255, 'linear');
    %     usecolormap = buildcmap('kryww');
% else
    clims = [-8 8]; %zscore power
    %     usecolormap = 'jet';
    usecolormap = flipud(cbrewer('div', 'RdBu', 255, 'linear')); %red high white neutral blue low
    %         usecolormap = flipud(usecolormap);
    % %     calcfunction = 'power';
% end
% clims = [0 1]; %[0 .7]
% usecolormap = 'jet';

% usecolormap = buildcmap('rwb');
srate = 1500;
win = [-2 2]; %in seconds. The morlet wavelet needs to be centered at time zero.. so always make this window symmetric about zero
plotwin = [-.5 .5];
basewin = plotwin(2); %baseline window to use for 1/f normalization
plottimeWin = plotwin(1):1/srate:plotwin(2);
%baseline indices. length of 1/2 plot window, centered on first plot ind
% i.e. 0.5 sec basewindow starting at .25 sec before the start of the
% plotwindow
baselineind = ([(win(2)-plotwin(2))-basewin/2 (win(2)-plotwin(2))+basewin/2]).*srate;
% basedur = abs(win(1))+abs(win(2));
% baselineind = ([1 (srate*basedur)+1]);
% p-value
pval = 0.05;
% convert p-value to Z value
zval = abs(norminv(pval));
% number of permutations
n_permutes = 50; %1000 takes 24 hours for JZ1 wtrack days 1:6
OvsIind = 6;
outboundInd = 4;
inboundInd = 5;
CvsMind = 7;
correctOutInd = 2;
mistakeOutInd = 3;
%% ---------------- Data Filters --------------------------
animals = {'D12', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
days = [];
% animals = {'D10'};
% days = 1:12; %1:12
% animals = {'D12'};
% days = 1:7; %1:7
% animals = {'D13'};
% days = 1:10; %1:10
% animals = {'JZ2'};
% days = 1:5; %1:5
% animals = {'JZ3'};
% days = [1 2 5:10]; %[1 2 5:10]
% animals = {'JZ4'};
% days = [1:10]; %[1:10]
% animals = {'JZ1'};
% days = [1:6];
filtfunction = 'riptriglfp';
LFPtypes = {'eeg', 'ripple', 'theta', 'lowgamma', 'fastgamma'};
% LFPrangesHz = {'1-400', '140-250', '6-9', '20-50', '65 - 95'}; %need to make this a lookup from the filter mat
eventtype = 'rippleskons';
epochEnvironment = {'wtrack'};% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
% epochType = 'run';
eventSourceArea = 'ca1';
% ripAreas = {'ca1', 'mec', 'por'};
% ntAreas = {'ca1', 'sub', 'mec', 'por', 'v2l', 'ref'};
% consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 5;        % STD. how big your ripples are
% exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
% minvelocity = 0;
% maxvelocity = 4;
% outputDirectory = '/typhoon/droumis/analysis';
%% ----------------- wavelet parameters ------------------------------------
min_freq =  2;
max_freq = 30;
frexres = 1; %frequency resolution default 4
num_frex = floor((max_freq - min_freq)/frexres);
% set range for variable number of wavelet cycles
% range_cycles = [min_freq*win(2) max_freq*win(2)]; % is there a principled way to do this?
% more wavelets for the same frequency will reduce the temporal precision.
% more wavelets for the same frequency will increase the frequency precision
% range_cycles = [4 floor(max_freq/10)];
range_cycles = [2 10];
% win = [[1,-1*win(1)]; %seconds
timeWin = win(1):1/srate:win(2);
% frex = logspace(log10(min_freq),log10(max_freq),num_frex);
frex = linspace(min_freq,max_freq,num_frex+1);
nWavecycles = linspace(range_cycles(1),range_cycles(end),num_frex); %number of wavelet cycles per freq
% nWavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex); %number of wavelet cycles per freq
% nWavecycles = repmat(range_cycles,1,num_frex);
% figure
% plot(nWavecycles,frex, 'o-')
% ylabel('frequency')
% xlabel('# cycles')
half_wave_size = (length(timeWin)-1)/2;
%% ---------------- Paths and Title strings ---------------------------------------------------
investInfo = animaldef(lower('Demetris'));
figdirectory = sprintf('%s%s/', investInfo{4}, plotoutputtype);
% filenamesave = sprintf('%s%s_%s_%s_%s', eventSourceArea, eventtype, epochEnvironment, cell2mat(animals), cell2mat(LFPtypes));
filenamesave = sprintf('%s%sSTD%d_%s_%s_%s_D%s', eventSourceArea, eventtype, minstdthresh, strjoin(epochEnvironment,'-'),...
    cell2mat(animals), cell2mat(LFPtypes),strjoin(arrayfun(@(x) num2str(x,'%-2d'),days,'UniformOutput',false),'-'));
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
resultfilename = sprintf('%s_%s%s_%s_%s_D%s', calcfunction, eventSourceArea, eventtype, strjoin(epochEnvironment,'-'), cell2mat(animals), strrep(num2str(days, '%-2d'),' ', '-'));
filenameTitle = strrep(resultfilename,'_', ' ');
iEpTetDataTypeFields = {'all', 'corrOut', 'mistOut', 'outB', 'inB','outB-inB', 'corrOut-mistOut'};

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
    eval(['F = setfilterfunction(F, [''dfa_'' filtfunction], {[eventSourceArea eventtype],' strjoin(arrayfun(@(x) sprintf('''%s'',', cell2mat(x)), LFPtypes,'UniformOutput',false)) '},''eventtype'', [eventSourceArea eventtype], ''LFPtypes'', LFPtypes, ''win'', win);']);
    %     eval(['F = setfilterfunction(F, [''dfa_'' filtfunction], {[eventSourceArea eventtype],' strjoin(reshape(repmat(arrayfun(@(x) sprintf('''%s'',', cell2mat(x)), LFPtypes,'UniformOutput',false), 2, 1), 1,length(LFPtypes)*2)) '},''eventtype'', [eventSourceArea eventtype], ''LFPtypes'', LFPtypes);']);
    tic
    F = runfilter(F);
    F(1).filterTimer = toc; F(1).filterTimer
    F(1).worldDateTime = clock;
    F(1).dataFilter = struct('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator, 'filename', filename);
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s%s/%s',investInfo{2},filtfunction, filename))
end
%% ---------------- Load Results Output ---------------------------------------------------
if loadResultsOutput == 1;
    load(sprintf('%s%s/%s_%dHz-%dHz_rCyc%d-%d.mat', investInfo{3}, calcfunction, resultfilename, min_freq, max_freq,range_cycles(1), range_cycles(2)));
end
%% ---------------- Calc Event State ---------------------------------------------------
if calcEventState
    ixpc.eventstate = [];
    ixpc.animals = [];
    ixpc.index = [];
    for ianimal = 1:length(F)
        animalinfo = animaldef(lower(animals{ianimal}));
        animalID = animalinfo{1,3}; %use anim prefix for name
        ixpc.animals = [ixpc.animals; {animalID}];
        iEpTetData = [];
        load(sprintf('%s%s%s.mat',animalinfo{1,2}, animalID, behavestruct));
        eventstate = [];
        eventstate.eventStartTimes = [];
        eventstate.state = [];
        eventstate.fields = [];
        ixpc.index{ianimal} = [];
        %get the animal state for every event
        for iday = 1:length(days)
            day = days(iday);
            eps = find(~cellfun(@isempty,{F(ianimal).output{day}.index})); %get nonempty eps
            ixpc.index{ianimal} = [ixpc.index{ianimal}; F(ianimal).output{day}(eps(1)).index];
            %load linpos per day
            load(sprintf('%s%s%s%02d.mat',animalinfo{1,2}, animalID, 'linpos', day));
            for iep = eps;
                iEpTetData = cat(2,iEpTetData,F(ianimal).output{day}(iep).data{1,1}); %cat the event snips across eps, days
                % get trajectory, segment, and linddist
                eventStartIndices = F.output{day}(iep).eventStartIndices;
                LFPtimes = F.output{day}(iep).LFPtimes;
                eventStartTimes = LFPtimes(eventStartIndices);
                statematrix = linpos{day}{iep}.statematrix;
                statemat.data = [statematrix.time statematrix.traj statematrix.segmentIndex statematrix.lindist];
                statemat.fields = 'postime traj segment lindist';
                stateindex = knnsearch(statemat.data(:,1), eventStartTimes);
                eventTrajs = statemat.data(stateindex,:);
                
                %get performance state
                trialIO = BehaveState.statechanges{day}{iep}.statechangeseq;
                trialIOfields = BehaveState.statechanges{day}{iep}.fields;
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
            disp(sprintf('%d max event-time offset between pos and lfp times is more than 33ms (1 cam frame).. removing offset events', max(abs(diff(eventstate.state(:,1:2),[],2)))))
            removeevents = find(abs(diff(eventstate.state(:,1:2),[],2)) > 0.033);
            removevec = ones(length(eventstate.state(:,1)),1);
            removevec(removeevents) = 0;
            eventstate.state = eventstate.state(logical(removevec),:);
            iEpTetData = iEpTetData(logical(removevec));
        end
        if length(eventstate.state(:,1)) ~= length(iEpTetData)
            error('mismatch between number of state info and lfp sanples')
        end
        ixpc.eventstate = eventstate;
        ixpc.allNTDataCat{ianimal} = cat(3, iEpTetData{:});
        %         iEpTetData = iEpTetData;
    end
    %% ---------------- Save eventState Output ---------------------------------------------------
    if saveEventState == 1;
        if ~isdir(sprintf('%s%s/', investInfo{3}, calcfunction));
            mkdir(sprintf('%s%s/', investInfo{3}, calcfunction));
        end
        save(sprintf('%s%s/%s_%dHz-%dHz_rCyc%d-%d.mat', investInfo{3}, calcfunction, resultfilename, min_freq, ...
            max_freq,range_cycles(1), range_cycles(2)), 'ixpc','-v7.3');
        disp(sprintf('%s saved', resultfilename))
    end
end

%% Calculate inter-trial phase clustering with morlet wavelet convolution


if calculateIXPC
    tic
    %     datafields = {'phaseoutput','poweroutput'};
    ixpc.phaseoutput = [];
    ixpc.poweroutput = [];
    ixpc.type = calcfunction;
    ixpc.range_cycles = range_cycles;
    ixpc.max_freq = max_freq;
    ixpc.min_freq = min_freq;
    ixpc.num_frex = num_frex;
    ixpc.win = win;
    ixpc.srate = srate;
    % split the events out into groups based on correct out / mistake out
    corrOutEventInd = find(ixpc.eventstate.state(1:end,6) == 1);
    mistOutEventInd = find(ixpc.eventstate.state(1:end,7) == 1);
    outBEventInd = [corrOutEventInd; mistOutEventInd];
    inBEventInd = setdiff([1:length(ixpc.eventstate.state(1:end,1))], outBEventInd)';
    iEpTetIndsType = {[1:length(ixpc.eventstate.state(1:end,1))], corrOutEventInd, mistOutEventInd, outBEventInd, inBEventInd};
    
    for ianimal = 1:length(F)
        nevents = size(ixpc.allNTDataCat{ianimal},3);
        nWave = length(timeWin);
        nNTrodes = 1;
        nsamps = size(ixpc.allNTDataCat{ianimal},2);
        nData = nsamps*nevents;
        nConv = nWave+nData-1; % length of the result of convolution.
        nConv2pow = 2^nextpow2(nConv);
        zpad2pow = nConv2pow - nConv;
        ntets = size(ixpc.allNTDataCat{ianimal},1);
        for introde = 1:ntets
            intDataCat = squeeze(ixpc.allNTDataCat{ianimal}(introde,:,:)); %squeeze will make nsampes x nevents..
            % concat reshape all the events (peri-rip snips) for speed and to reduce edge artifacts, then take the fft.
            % fft with the arg of nConv ensures that matlab will zero-pad the half-wavelet duration on either side of the data series
            %             dataY = fft(reshape(iEpTetDataCat,nNTrodes,nData),nConv,2);
            dataY = fft(reshape(intDataCat,nNTrodes,nData),nConv2pow,2); % reshape reshapes row-wise into a 1 x d vector
%             dataY{introde} = fft(reshape(squeeze(ixpc.allNTDataCat{ianimal}(introde,:,:)),nNTrodes,nData),nConv2pow,2); % reshape reshapes row-wise into a 1 x d vector
            %         %% fix ispc
            %         if calculateISPC
            %             % initialize output time-frequency data
            %             %                 ixpc.output{ianimal}{iday} = zeros(nchoosek(size(iEpTetDataCat,1),2), nsamps, num_frex);
            %             ixpc.phaseoutput{ianimal}{idatatype} = zeros(nchoosek(size(iEpTetDataCat,1),2), nsamps, num_frex);
            
            %% loop over frequencies
            for fi=1:num_frex
                % create wavelet and get its FFT
                % creating the morlet wavelet by combingin the complex sine wave and the gaussian
                sine_wave = exp(2*1i*pi*frex(fi).*timeWin); %make a complex sine wave
                bn = nWavecycles(fi)/(2*pi*frex(fi)); % std of the gaussian of the morlet wavelet. dependent on the freq and #cycles
                gaus_win = exp(-timeWin.^2./(2*bn^2)); %the gaussian
                wavelet = sine_wave .* gaus_win;
                waveletFFT = fft(wavelet,nConv2pow);
%                 waveletFFT{introde}{fi} = fft(exp(2*1i*pi*frex(fi).*timeWin) .* exp(-timeWin.^2./(2*nWavecycles(fi)/(2*pi*frex(fi))^2)),nConv2pow);
                %normalize wavelet to a maximum of 1. this will ensure that the units of convolution are the same as in the original data.
                waveletFFT = waveletFFT ./ max(waveletFFT);
%                 waveletFFT{introde}{fi} = waveletFFT{introde}{fi} ./ max(waveletFFT{introde}{fi});
                % run convolution (filtering) : Time-domain convolution in the frequency domain... because itz so much faster
                astmp = bsxfun(@times,dataY,waveletFFT); %multiply the power spectrum of the data and the wavelet
%                 astmp{introde}{fi} = bsxfun(@times,dataY{introde},waveletFFT); %multiply the power spectrum of the data and the wavelet
                astmp = ifft(astmp,nConv2pow,2); % take the inverse transform
                %trim off the length of half the wavelet at the beginning
                %and at the end. also trim the zero padding at the tail
                astmp = astmp(1,half_wave_size+1:end-half_wave_size-zpad2pow);
                as{introde}(:,:,fi) = reshape(astmp,nsamps,nevents); %reshape it back to ntrodes by nsamples by nevents
                phdata{introde}(:,:,fi) = angle(as{introde}(:,:,fi)); %get phase component of the analytic signal
                %             %% ISPC -- NEED TO FIX :(
                %             if calculateISPC % time series diff of all possible ntrode pairs.
                %                 %the idea here is to treat rows as pairs instead of individual ntrodes, as in ITPC
                %                 ispc = zeros(nchoosek(size(phdata{introde}{fi},1),2), nsamps, nevents);
                %                                 m = logical(tril(ones(size(phdata{introde}{fi}(:,:,1),1)),-1)); %get indices of non-duplicates (below comb triangle)
                %                                 [brow, bcol] = find(m); %get linear index of ntrode combination indices
                %                                 indices = [F(ianimal).output{day}(eps(1)).index(bcol,:) F(ianimal).output{day}(eps(1)).index(brow,:)]; %convert to ntrodeID
                %                 for w = 1:size(phdata{introde}{fi},3); %loop through each event plane (ntrode x samples X event#)
                %                     iprm = permute(phdata{introde}{fi}(:,:,w),[3 2 1]); %transpose
                %                     B = bsxfun(@minus,phdata{introde}{fi}(:,:,w),iprm); %subtract across each NTrode-choose-two comb of phase time series
                %                     B = reshape(B,[],size(phdata{introde}{fi}(:,:,w),2)); %reshape back into 2D
                %                     B = B(m(:),:); %get rid of duplicates
                %                     ispc(:,:,w) = B; %save result
                %                 end
                %                 clear phasedata
                %                 phdata{introde}{fi} = ispc;
                %             end
                %             ixpc.index{ianimal}{fi} = indices;
                disp(sprintf('nt%d freq %d of %d',introde,fi,num_frex));
                
            end
            
            
            ixpc.frex{ianimal}{introde} = frex;
            %% compute ITPC-Xtrial and Power-Xtrial with frequency-wise zscore
            % power is the squared magnitude from origin to the location in complex space
            % zscoring using the entire window..this is conservative as it includes the post-rip activity into the mean and std
            ixpc.basemeanpower{ianimal}{introde}(1,1,:) = mean(mean(abs(as{introde}(baselineind(1):baselineind(2),iEpTetIndsType{1},:)).^2,2),1);
            ixpc.basestdpower{ianimal}{introde}(1,1,:) = std(reshape(abs(as{introde}(baselineind(1):baselineind(2),iEpTetIndsType{1},:)).^2,[],1,num_frex),1);
            powrtmp = cellfun(@(x) mean(abs(as{introde}(:,x,:)).^2,2), iEpTetIndsType, 'un', 0);
            ixpc.poweroutput{ianimal}{introde}(:,:,:) = cat(2,powrtmp{:});
            ixpc.basenormpoweroutput{ianimal}{introde} = bsxfun(@rdivide, bsxfun(@minus, ixpc.poweroutput{ianimal}{introde},ixpc.basemeanpower{ianimal}{introde}), ixpc.basestdpower{ianimal}{introde});
            
            ixpc.basemeanphase{ianimal}{introde}(1,1,:) = mean(abs(mean(exp(1i*phdata{introde}(baselineind(1):baselineind(2),iEpTetIndsType{1},:)),2)),1);
            ixpc.basestdphase{ianimal}{introde}(1,1,:) = std(reshape(abs(mean(exp(1i*phdata{introde}(baselineind(1):baselineind(2),iEpTetIndsType{1},:)),2)),[],1,num_frex),1);
            itpctmp = cellfun(@(x) abs(mean(exp(1i*phdata{introde}(:,x,:)),2)), iEpTetIndsType, 'un', 0);
            ixpc.phaseoutput{ianimal}{introde}(:,:,:) = cat(2, itpctmp{:});
            ixpc.basenormphaseoutput{ianimal}{introde} = bsxfun(@rdivide, bsxfun(@minus, ixpc.phaseoutput{ianimal}{introde},ixpc.basemeanphase{ianimal}{introde}), ixpc.basestdphase{ianimal}{introde});
            %% compute differential tf maps
            ixpc.phaseoutput{ianimal}{introde}(:,6,:) = ixpc.phaseoutput{ianimal}{introde}(:,4,:) - ixpc.phaseoutput{ianimal}{introde}(:,5,:);
            ixpc.poweroutput{ianimal}{introde}(:,6,:) = ixpc.poweroutput{ianimal}{introde}(:,4,:) - ixpc.poweroutput{ianimal}{introde}(:,5,:);
            ixpc.phaseoutput{ianimal}{introde}(:,7,:) = ixpc.phaseoutput{ianimal}{introde}(:,2,:) - ixpc.phaseoutput{ianimal}{introde}(:,3,:);
            ixpc.poweroutput{ianimal}{introde}(:,7,:) = ixpc.poweroutput{ianimal}{introde}(:,2,:) - ixpc.poweroutput{ianimal}{introde}(:,3,:);
            ixpc.basenormphaseoutput{ianimal}{introde}(:,6,:) = ixpc.basenormphaseoutput{ianimal}{introde}(:,4,:) - ixpc.basenormphaseoutput{ianimal}{introde}(:,5,:);
            ixpc.basenormpoweroutput{ianimal}{introde}(:,6,:) = ixpc.basenormpoweroutput{ianimal}{introde}(:,4,:) - ixpc.basenormpoweroutput{ianimal}{introde}(:,5,:);
            ixpc.basenormphaseoutput{ianimal}{introde}(:,7,:) = ixpc.basenormphaseoutput{ianimal}{introde}(:,2,:) - ixpc.basenormphaseoutput{ianimal}{introde}(:,3,:);
            ixpc.basenormpoweroutput{ianimal}{introde}(:,7,:) = ixpc.basenormpoweroutput{ianimal}{introde}(:,2,:) - ixpc.basenormpoweroutput{ianimal}{introde}(:,3,:);
%         end
%         for introde = 1:ntets
            %% permutation testing... this will take TIME and a shit ton of RAM.. go get covfefe
            if runPermTest
                disp(sprintf('========== running perm tests ============= nt%d perm x %d',introde, n_permutes))
                % outbound vs inbound 
                ixpc.powerzmask{ianimal}{introde}{OvsIind} = []; ixpc.phasezmask{ianimal}{introde}{OvsIind} = [];
                ixpc.MCpower_minmax{ianimal}{introde}{OvsIind} = []; ixpc.MCphase_minmax{ianimal}{introde}{OvsIind} = [];
                [ixpc.powerzmask{ianimal}{introde}{OvsIind}, ixpc.phasezmask{ianimal}{introde}{OvsIind},...
                    ixpc.MCpower_minmax{ianimal}{introde}{OvsIind},ixpc.MCphase_minmax{ianimal}{introde}{OvsIind}]...
                = ITPCpermtest(ixpc.poweroutput{ianimal}{introde}(:,OvsIind,:), ixpc.phaseoutput{ianimal}{introde}(:,OvsIind,:), phdata{introde}, as{introde}, iEpTetIndsType{outboundInd}, iEpTetIndsType{inboundInd}, n_permutes);
                % correct-outbound vs mistake-outbound
                ixpc.powerzmask{ianimal}{introde}{CvsMind} = []; ixpc.phasezmask{ianimal}{introde}{CvsMind} = [];
                ixpc.MCpower_minmax{ianimal}{introde}{CvsMind} = []; ixpc.MCphase_minmax{ianimal}{introde}{CvsMind} = [];
                [ixpc.powerzmask{ianimal}{introde}{CvsMind}, ixpc.phasezmask{ianimal}{introde}{CvsMind},...
                    ixpc.MCpower_minmax{ianimal}{introde}{CvsMind},ixpc.MCphase_minmax{ianimal}{introde}{CvsMind}]...
                = ITPCpermtest(ixpc.poweroutput{ianimal}{introde}(:,CvsMind,:), ixpc.phaseoutput{ianimal}{introde}(:,CvsMind,:), phdata{introde}, as{introde}, iEpTetIndsType{correctOutInd}, iEpTetIndsType{mistakeOutInd}, n_permutes);
                % ixpc = ITPCpermtest(ixpc, phdata{introde}, as{introde}, ianimal, introde, 7, iEpTetIndsType{2}, iEpTetIndsType{3}, n_permutes);
                %% cluster-based multiple comp scratch
                %                 for iprm = 1:n_permutes;
                % %                     % take each permutation map, and transform to Z
                % %                     ipermmap = phasepermOutput(:,iprm,:);
                % %                     izpermmap = (ipermmap-ixpc.phasemean_h0{introde}{6}(:,1,:))./ixpc.powerstd_h0{introde}{6}(:,1,:);
                % %                     % threshold image at p-value
                % %                     izpermmap(abs(izpermmap)<zval) = 0;
                %                     % find clusters (need image processing toolbox for this!)
                %                     izpermmap = squeeze(zpermmaps(:,iprm,:));
                %                     izpermislands = bwconncomp(izpermmap);
                %                     if ~isempty(cell2mat(cellfun(@numel,izpermislands.PixelIdxList))); %numel(izpermislands.PixelIdxList)>0
                %                         % count sizes of clusters
                % %                         tempclustsizes = cellfun(@length,izpermislands.PixelIdxList);
                %                         % store size of biggest cluster
                %                         maxclustsize = max(cell2mat(cellfun(@numel,izpermislands.PixelIdxList)));
                % %                         [biggest,idx] = max(numPixels);
                % %                         BW(CC.PixelIdxList{idx})
                % %                         ixpc.MCmax_clust_size{ianimal}{introde} = [ixpc.MCmax_clust_size{ianimal}{introde}; max(tempclustsizes)];
                %                          ixpc.MCmax_clust_size{ianimal}{introde} = [ixpc.MCmax_clust_size{ianimal}{introde}; maxclustsize];
                %                     end
                %                 end
            end
        end
        ixpc.datatypes{ianimal} = [];
        ixpc.datatypes{ianimal} = iEpTetDataTypeFields;
        
% %         
%         demec = [2 3 7 4 6 14 1 8 9 10];
% %         sumec = [1 8 9 10];
%         por = [13 15 12];
%         hc = [16:30];
%         areascat = {demec, por, hc};
%         for iar = 1:length(areascat)
%             %phase
%             areaphaseoutputtmp = arrayfun(@(x) ixpc.basenormphaseoutput{ianimal}{x}(:,6,:), areascat{iar}, 'un', 0);
%             ixpc.areaphasemean{ianimal}{iar} = mean(cat(2,areaphaseoutputtmp{:}), 2);
%             
%             areaphasezmasktmp = arrayfun(@(x) ixpc.phasezmask{ianimal}{x}{6}, areascat{iar}, 'un', 0);
%             ixpc.areaphasezmaskmean{ianimal}{iar} = mean(cat(2,areaphasezmasktmp{:}), 2);
%             % power
%             areapoweroutputtmp = arrayfun(@(x) ixpc.basenormpoweroutput{ianimal}{x}(:,6,:), areascat{iar}, 'un', 0);
%             ixpc.areapowermean{ianimal}{iar} = mean(cat(2,areapoweroutputtmp{:}), 2);
%             
%             areapowerzmasktmp = arrayfun(@(x) ixpc.powerzmask{ianimal}{x}{6}, areascat{iar}, 'un', 0);
%             ixpc.areapowerzmaskmean{ianimal}{iar} = mean(cat(2,areapowerzmasktmp{:}), 2);
% 
%         end
%         ixpc.areas{ianimal} = {'mec', 'por', 'hc'};
    end
    toc
end
%% ---------------- Save RESULTS Output ---------------------------------------------------
if saveResultsOutput == 1;
    if ~isdir(sprintf('%s%s/', investInfo{3}, calcfunction));
        mkdir(sprintf('%s%s/', investInfo{3}, calcfunction));
    end
    save(sprintf('%s%s/%s_%dHz-%dHz_rCyc%d-%d.mat', investInfo{3}, calcfunction, resultfilename, min_freq, max_freq,range_cycles(1), range_cycles(2)), 'ixpc','-v7.3');
    disp(sprintf('%s saved', resultfilename))
end

%% ---------------- plot ITPC---------------------------------------------------------------------------------------------
if plotITPC
    %     clear F %save space on memory
    warning('off', 'MATLAB:contour:ConstantData')
    position = [.1 .1 .9 .8];
    SpacingHorizontal = 0.01;
    SpacingVertical = 0.02;
    Spacing = 0.00;
    Padding = 0.0;
    MarginLeft = 0.04;
    MarginRight = 0.04;
    MarginTop = 0.09;
    MarginBottom =  0.08;
    for ianimal = 1:length(ixpc.animals)
        animalinfo = animaldef(lower(animals{ianimal}));
        animalID = animalinfo{1,3}; %use anim prefix for name
        FFanimdir =  sprintf('%s',animalinfo{1,2});
        %% ---- loadtetinfostruct ----
        load([FFanimdir, animalID, 'tetinfo']);
        tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
        %         for iday = days
        %             ntrodesIndices = ixpc.index{ianimal}{iday};
        ntrodesIndices = ixpc.index{ianimal};
        ntets = size(unique(ntrodesIndices(:,3),'stable'),1);
        tets = unique(ntrodesIndices(:,3));
        if ntets ~= length(tets);
            error('data doesnt match indices')
        end
        %% ---- reorder the LFP traces by suparea|area|subarea tags (in that priority) ----
        [~, tagIndMap] = ismember(tets,tetinfoAll.index(:,3), 'rows');
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
        for iDT = 1:(length(ixpc.datatypes{ianimal}))%+2
            if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
                ifig1 = figure('Visible','off','units','normalized','position',position);
            else
                ifig1 = figure('units','normalized','position',position);
            end
            set(gcf,'color','white')
            sfrows = floor(sqrt(ntets));
            sfcols = ceil(sqrt(ntets));
            
            %% ---- loop across all tets for this day ----
            for introde = 1:ntets;
                introdeID = ntrodesIndices(numsumSortInds(introde),3);
                isupareatag = ntrodeTags{numsumSortInds(introde)}.suparea;
                iareatag = ntrodeTags{numsumSortInds(introde)}.area;
                isubareatag = ntrodeTags{numsumSortInds(introde)}.subarea;
                iNTcolor = icolors(introde,:);
                intfig = subaxis(sfrows,sfcols,introde, 'SpacingVertical', SpacingVertical, 'SpacingHorizontal', SpacingHorizontal, ...
                    'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight, 'MarginTop', MarginTop, ...
                    'MarginBottom', MarginBottom);
                %% ---- Main Plotting ----
                %                 contourf(intfig,timeWin,frex,squeeze(ixpc.output{ianimal}{iday}(numsumSortInds(introde),:,:))',num_frex,'linecolor','none');%, )
                %                 if iDT == 6
                %                     eval(sprintf('ixpc2plot = squeeze(ixpc.%soutput{ianimal}{4}(numsumSortInds(introde),:,:))'' - squeeze(ixpc.%soutput{ianimal}{5}(numsumSortInds(introde),:,:))'';', plotoutputtype, plotoutputtype));
                %
                % %                     ixpc2plot = (squeeze(ixpc.phaseoutput{ianimal}{4}(numsumSortInds(introde),:,:))' - squeeze(ixpc.phaseoutput{ianimal}{5}(numsumSortInds(introde),:,:))');
                %                 elseif iDT == 7
                %                     eval(sprintf('ixpc2plottmp = squeeze(ixpc.%soutput{ianimal}{2}(numsumSortInds(introde),:,:))'' - squeeze(ixpc.%soutput{ianimal}{3}(numsumSortInds(introde),:,:))'';', plotoutputtype, plotoutputtype));
                %                 else
                %                                 eval(sprintf('ixpc2plottmp = squeeze(ixpc.%soutput{ianimal}{iDT}(numsumSortInds(introde),:,:))'';', plotoutputtype));
                
                
                %                     eval(sprintf('basemeanFULL = squeeze(ixpc.basemean%s{ianimal}{numsumSortInds(introde)}(:,1,:))'';', plotoutputtype));
                %                     eval(sprintf('basestdFULL = squeeze(ixpc.basestd%s{ianimal}{numsumSortInds(introde)}(:,1,:))'';', plotoutputtype));
                %                     ixpc2plotFULL = (bsxfun(@rdivide, bsxfun(@minus, ixpc2plotFULL, basemeanFULL'), basemeanFULL'));
                %                 end
                if 0
%                     eval(sprintf('ixpc2plotFULL = squeeze(ixpc.area%smean{ianimal}{introde})'';', plotoutputtype));
                    eval(sprintf('ixpc2plotFULL = squeeze(ixpc.area%szmaskmean{ianimal}{introde})'';', plotoutputtype));
                    midpoint = (size(ixpc2plotFULL,2)-1)/2; %get middle index of window
                    ixpc2plot = ixpc2plotFULL(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
                    [~,bn] = contourf(intfig,plottimeWin,frex(1:end-1),ixpc2plot,100,'linecolor','none'); %
                    hold on;
%                     eval(sprintf('zmaskFULL = squeeze(ixpc.area%szmaskmean{ianimal}{introde})'';', plotoutputtype));
% %                     zmaskFULL = squeeze(ixpc.areaphasezmaskmean{ianimal}{introde})';
%                     midpoint = (size(ixpc2plotFULL,2)-1)/2; %get middle index of window
%                     zmask2plot = zmaskFULL(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
                    zmask = ixpc2plot;
                    zmask(abs(ixpc2plot)<zval) = 0;
                    [~,h] = contour(intfig,plottimeWin,frex(1:end-1),logical(zmask),1);
                    h.LineColor = [.85 .85 .85];
                else
                    
%                     eval(sprintf('ixpc2plotFULL = squeeze(ixpc.%soutput{ianimal}{numsumSortInds(introde)}(:,iDT,:))'';', plotoutputtype));
%                     midpoint = (size(ixpc2plotFULL,2)-1)/2; %get middle index of window
%                     ixpc2plot = ixpc2plotFULL(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
%                     [~,bn] = contourf(intfig,plottimeWin,frex(1:end-1),ixpc2plot,100,'linecolor','none'); %
                                    eval(sprintf('BaseNormixpc2plotFULL = squeeze(ixpc.basenorm%soutput{ianimal}{numsumSortInds(introde)}(:,iDT,:))'';', plotoutputtype));
                                    midpoint = (size(BaseNormixpc2plotFULL,2)-1)/2; %get middle index of window
                                    BaseNormixpc2plot = BaseNormixpc2plotFULL(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
                                    [~,bn] = contourf(intfig,plottimeWin,frex(1:end-1),BaseNormixpc2plot,100,'linecolor','none'); %
                end
                hold on;
                %
                
                if iDT == 6 || iDT == 7
                    eval(sprintf('ixpc2plotFULL = squeeze(ixpc.%soutput{ianimal}{numsumSortInds(introde)}(:,iDT,:))'';', plotoutputtype));
                    eval(sprintf('Zmaskixpc2plotFULL = squeeze(ixpc.%szmask{ianimal}{numsumSortInds(introde)}{iDT})'';', plotoutputtype));
                    eval(sprintf('MCminmax = abs(ixpc.MC%s_minmax{ianimal}{numsumSortInds(introde)}{iDT})'';',plotoutputtype));
                    ixpc2plot = ixpc2plotFULL(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
                    Zmaskixpc2plot = Zmaskixpc2plotFULL(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
                    
                    zmask = Zmaskixpc2plot;
                    zmask(abs(Zmaskixpc2plot)<zval) = 0;
                    
                    [~,h] = contour(intfig,plottimeWin,frex(1:end-1),logical(zmask),1);
                    h.LineColor = [.85 .85 .85];
                    %                     alphable = findobj(h, '-property', 'FaceAlpha');
                    %                     set(alphable, 'FaceAlpha', .5);
                    hold on;
                    MCthresh = MCminmax(ceil(length(MCminmax)*(1-pval)));
                    MCthreshmap = ixpc2plot;
                    MCthreshmap(abs(ixpc2plot)<MCthresh) = 0;
                    [~,mc] = contour(intfig,plottimeWin,frex(1:end-1),logical(MCthreshmap),1);
                    mc.LineColor = [.6 .6 .6];
                    caxis([-2 2]);
                else
                    caxis([-3 3]);
                end
                % figure
                %                 imagesc(ixpc2plot);
                %                 patch([-.8, .8, .8, -.8], [-10 -10 80 80], iNTcolor, 'edgecolor','none')
                hold on;
                %                 if strcmp(plotoutputtype, 'phase')
                set(gca,'ydir','normal','xlim',[plotwin(1) plotwin(2)], 'ylim',[frex(1) frex(end-1)])
                %                     caxis([-10 10])
                
                %                 caxis([-1 1]);
                %                 colorbar
                %                 caxis('auto')
                %                 else
                %                     set(gca,'ydir','normal','xlim',[-.5 .5])
                %                     calcfunction = 'power';
                %                 end
                %                 set(gca, 'YScale', 'log')
                colormap(usecolormap)
                
                
                if mod(introde, sfcols) ~= 1;
                                        set(gca, 'yscale','log','ytick',[])	
%                                         set(gca, 'yscale','log','ytick',ceil(logspace(log10(min(frex)),log10(max(frex)),6)))
                    set(gca, 'FontSize',8, 'FontName', 'Arial');%
                else
                                        set(gca, 'yscale','log','ytick',ceil(logspace(log10(min(frex)),log10(max(frex)),6)))
                    %                     set(gca, 'YTick',[round(linspace(min(frex),max(frex),6))])
                    set(gca, 'FontSize',8, 'FontName', 'Arial');%
                end
                if introde <= (sfrows-1)*sfcols;
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
            %% crap
            %                 zmaskinds = find(abs(Zmaskixpc2plot)<zval);
            %                 zmaskscaled = ones(size(Zmaskixpc2plot,1), size(Zmaskixpc2plot,2));
            %                 zmaskscaled(zmaskinds) = .1;
            %                 bn = pcolor(BaseNormixpc2plot);
            %                     bn.AlphaData = zmaskscaled; % make all zero bins transparent
            %                     bn.AlphaDataMapping = 'scaled';
            %                     bn.EdgeColor = 'none';
            %                     bn.FaceAlpha = 'interp';
            %                     bn.FaceColor = 'interp';
            %
            %                 hold on;
            %
            % %                 zmaskscaled;
            % %                 Bimage = get(bn, 'CData');
            %                 set(bn, 'AlphaData', zmaskscaled)%, 'AlphaDataMapping', 'none')
            %
            % %                 if strcmp(plotoutputtype, 'phase')
            %                     ixpc2plot = ixpc2plottmp2;
            %                 elseif strcmp(plotoutputtype, 'power') % ONLY DO THIS IF NOT A DIFFERENTIAL CONDITION! the diff 1/f will normalize itself
            %                     %                 ixpc2plottmp2 = zscore(ixpc2plottmp, [], 1);
            %                     %                 baseline = mean(ixpc2plottmp(:,1:midpoint+plotwin(1)*srate),2); %get baseline (start to middle index)
            %                     baselinemean = mean(ixpc2plotFULL,2);
            %                     baselinestd = std(ixpc2plotFULL,[],2);
            %                     %                 baselinemean = mean(ixpc2plottmp2,2);
            %                     %                 baselinestd = std(ixpc2plottmp2,[],2);
            %                     %                 ixpc2plot = 100 * bsxfun(@rdivide, bsxfun(@minus,ixpc2plottmp2,baselinemean), baselinemean); % normalize to percent change from baseline
            %                     ixpc2plot = bsxfun(@rdivide, bsxfun(@minus,ixpc2plot,baselinemean), baselinestd); % normalize to percent change from baseline
            %                 else
            %                     error('must be power or phase')
            %                 end
            %                intfig = figure
            
            %                 if iDT == 6
            % %                     % plot the distribution of values outside the MC null distribution
            % %                     [n1, x1] = hist(reshape(ixpc.MCpower_minmax{ianimal}{numsumSortInds(introde)},numel(ixpc.MCpower_minmax{ianimal}{numsumSortInds(introde)}), 1),100);
            % %                     [n2, x2] = hist(reshape(ixpc.poweroutput{ianimal}{numsumSortInds(introde)}(:,6,:),numel(ixpc.poweroutput{ianimal}{numsumSortInds(introde)}(:,6,:)), 1),100);
            % %                     figure; plot(x1,n1./max(n1),'r-')
            % %                     hold on; plot(x2,n2./max(n2),'b-'); hold off
            %                     a = sort(abs(reshape(ixpc.MCpower_minmax{ianimal}{numsumSortInds(introde)},numel(ixpc.MCpower_minmax{ianimal}{numsumSortInds(introde)}), 1)));
            %                     pvalthresh = a(ceil(length(a)*(1-pval)));
            % %                     imap = squeeze(ixpc.poweroutput{ianimal}{numsumSortInds(introde)}(:,6,:));
            %                     pthreshimap = ixpc2plot;
            %                     pthreshimap(abs(ixpc2plot)<pvalthresh) = 0;
            %                     contourf(intfig,plottimeWin,frex(1:end-1),pthreshimap,100,'linecolor','none'); %
            %                     contourf(intfig,plottimeWin,frex(1:end-1),ixpc2plot,100,'linecolor','none'); %
            %                     hold on;
            %                     contour(intfig,plottimeWin,frex(1:end-1),logical(pthreshimap),1,'k');
            %                     colorbar
            %                 end
            
            %% ---- super title and colorbar----
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig1);
            xlab = 'Time (s)';
            ylab = 'Frequency (Hz)';
            supxlabel = text(.5, .02, xlab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'Parent', sprtitleax,...
                'Units', 'normalized', 'horizontalAlignment', 'center');
            supylabel = text(.01, .5, ylab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'rotation', 90, ...
                'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
            if iDT == 6 || iDT == 7
                caxis([-2 2]);
            else
                caxis([-2 2]);
            end
            %                         caxis([-1 1])
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
            %             if iDT == 6
            %                 sprtit = sprintf('%s %s %s %s D%s %d-%dHz_rCyc%d-%d', plotoutputtype,'outB-inB', strjoin(epochEnvironment,'-'), animalID, strrep(num2str(days), '  ', '-'), min_freq,max_freq, range_cycles);
            %             elseif iDT == 7
            %                 sprtit = sprintf('%s %s %s %s D%s %d-%dHz_rCyc%d-%d', plotoutputtype,'outBCorr-outBMist', strjoin(epochEnvironment,'-'), animalID, strrep(num2str(days), '  ', '-'), min_freq,max_freq, range_cycles);
            %             else
            sprtit = sprintf('%s %s %s %s D%s %d-%dHz_rCyc%d-%d', plotoutputtype,ixpc.datatypes{ianimal}{iDT}, strjoin(epochEnvironment,'-'), animalID, strrep(num2str(days,'%-2d'), ' ', '-'), min_freq,max_freq, range_cycles);
            %             end
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
            clrbartit = text(posx(1)+posx(3)/2, posx(2)-MarginBottom/2, plotoutputtype, 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial','horizontalAlignment', 'center');
            
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
            end
            
            close all
        end
    end
end
%% ---------------- plot ISPC---------------------------------------------------------------------------------------------
%% ---------------- plot ISPC---------------------------------------------------------------------------------------------
%% ---------------- plot ISPC---------------------------------------------------------------------------------------------
if plotISPC
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
    for ianimal = 1:length(ixpc.animals)
        %% ---- loadtetinfostruct ----
        animalinfo = animaldef(lower(animals{ianimal}));
        animalID = animalinfo{1,3}; %use anim prefix for name
        FFanimdir =  sprintf('%s',animalinfo{1,2});
        load([FFanimdir, animalID, 'tetinfo']);
        tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
        for iday = days
            ntrodePairIndices = ixpc.index{ianimal}{iday};
            tetpairs = unique(ntrodePairIndices(:,[3 6]),'rows','stable');
            ntetPairs = size(tetpairs,1);
            %% ---- reorder the LFP traces by suparea|area|subarea tags (in that priority) ----
            ntrodePairIndicesLIN = [ntrodePairIndices(:,[1:3]); ntrodePairIndices(:,[4:6])];
            ntets = size(unique(ntrodePairIndicesLIN, 'stable'),1);
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
                contourf(timeWin,frex,squeeze(ixpc.phaseoutput{ianimal}{iday}(numsumSortInds(intrPair),:,:))',num_frex,'linecolor','none');
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
                sprtit = sprintf('%s %s %s D%d', calcfunction, epochEnvironment, animalID, iday);
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

