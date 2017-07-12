
% to do:
% ITPC - plot open field
% ISPC - preserve the indices

%%% First create the processed data structure with dfs_riptriglfp.m

close all
calcEventState = 0;
saveEventState = calcEventState;
calculateITPC = 1;
calculateISPC = 0;
calculateIXPC = calculateITPC || calculateISPC;
loadFilterOutput = calculateIXPC || calcEventState;
saveResultsOutput = calculateIXPC;
runPermTest = calculateIXPC;
loadResultsOutput = 1;
plotITPC = 0;
plotISPC = 0;
savefigs = 0;
pausefigs = 1;
calcfunction = 'ITPC'; %ISPC or ITPC
behavestruct = 'BehaveState';
plotoutputtype = 'power'; % 'power'
%% ---------------- plotting params --------------------------
colorSet = 'DR1';
if strcmp(plotoutputtype, 'phase')
    clims = [0 .7]; % mean ITPC goes from 0-1
    usecolormap = 'jet'; %cbrewer('seq', 'Greens', 255, 'linear');
    %     usecolormap = buildcmap('kryww');
else
    clims = [-1 6]; %zscore power
    usecolormap = 'jet';
    %     usecolormap = cbrewer('div', 'RdBu', 255, 'linear');
    %     usecolormap = flipud(usecolormap);
    % %     calcfunction = 'power';
end
% clims = [0 1]; %[0 .7]
% usecolormap = 'jet';

% usecolormap = buildcmap('rwb');
srate = 1500;
win = [-2 2]; %in seconds. The morlet wavelet needs to be centered at time zero.. so always make this window symmetric about zero
plotwin = [-.5 .5];
plottimeWin = plotwin(1):1/srate:plotwin(2);
% p-value
pval = 0.05;
% convert p-value to Z value
zval = abs(norminv(pval));
% number of permutations
n_permutes = 5; %1000

%% ---------------- Data Filters --------------------------
animals = {'JZ1'};
% animals = {'JZ1', 'D13'};
days = [1:6];
filtfunction = 'riptriglfp';
LFPtypes = {'eeg'};%, 'ripple'};%, 'theta', 'slowgamma', 'fastgamma'};
% LFPrangesHz = {'1-400', '140-250', '6-9', '20-50', '65 - 95'}; %need to make this a lookup from the filter mat
eventtype = 'rippleskons';
epochEnvironment = {'wtrack'};% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
% epochType = 'run';
eventSourceArea = 'ca1';
% ripAreas = {'ca1', 'mec', 'por'};
% ntAreas = {'ca1', 'sub', 'mec', 'por', 'v2l', 'ref'};
% consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 3;        % STD. how big your ripples are
% exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
% minvelocity = 0;
% maxvelocity = 4;
% outputDirectory = '/typhoon/droumis/analysis';

%% ---------------- Paths and Title strings ---------------------------------------------------
investInfo = animaldef(lower('Demetris'));
figdirectory = sprintf('%s%s/', investInfo{4}, plotoutputtype);
% filenamesave = sprintf('%s%s_%s_%s_%s', eventSourceArea, eventtype, epochEnvironment, cell2mat(animals), cell2mat(LFPtypes));
filenamesave = sprintf('%s%sSTD%d_%s_%s_%s_D%s', eventSourceArea, eventtype, minstdthresh, strjoin(epochEnvironment,'-'),...
    cell2mat(animals), cell2mat(LFPtypes),strjoin(arrayfun(@(x) num2str(x),days,'UniformOutput',false),'-'));
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
resultfilename = sprintf('%s_%s%s_%s_%s_D%s', calcfunction, eventSourceArea, eventtype, strjoin(epochEnvironment,'-'), cell2mat(animals), strrep(num2str(days), '  ', '-'));
filenameTitle = strrep(resultfilename,'_', ' ');
iEpTetDataTypeFields = {'all', 'corrOut', 'mistOut', 'outB', 'inB','outB-inB', 'corrOut-mistOut'};
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
            error('max event-time offset between pos and lfp times is more than 33ms (1 cam frame)')
        end
        if length(eventstate.state(:,1)) ~= length(iEpTetData)
            error('mismatch between number of state info and lfp sanples')
        end
        ixpc.eventstate = eventstate;
        ixcp.allNTDataCat{ianimal} = cat(3, iEpTetData{iEpTetIndsType{1}});
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

% wavelet parameters
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
if calculateIXPC
    tic
    fields = {'phaseoutput','poweroutput'};
    ixpc = rmfield(ixpc,fields);
    ixpc.type = calcfunction;
    ixpc.range_cycles = range_cycles;
    ixpc.max_freq = max_freq;
    ixpc.min_freq = min_freq;
    ixpc.num_frex = num_frex;
    ixpc.win = win;
    ixpc.srate = srate;
    % an even number of events will make nData even (since nsamps and nwave will always be...
    % odd due to the mirrored window. An even nData + an odd nwave minus 1 will be even... and therefore the ...
    % difference between this result and the next power of 2 will be even allowing an equal zero pad of odd number
    % and making things a lot easier... so always check if there's an
    % even number of events, and if not, remove the last one.
    evenEventnum = 0; %mod(size(ixpc.eventstate.state,1),2);
    % split the events out into groups based on correct out / mistake out
    corrOutEventInd = find(ixpc.eventstate.state(1:end-evenEventnum,6) == 1);
    mistOutEventInd = find(ixpc.eventstate.state(1:end-evenEventnum,7) == 1);
    outBEventInd = [corrOutEventInd; mistOutEventInd];
    inBEventInd = setdiff([1:length(ixpc.eventstate.state(1:end-evenEventnum,1))], outBEventInd)';
    iEpTetIndsType = {[1:length(ixpc.eventstate.state(1:end-evenEventnum,1))], corrOutEventInd, mistOutEventInd, outBEventInd, inBEventInd};
    %         iEpTetDataType{1} = iEpTetData;
    %         iEpTetDataType{2} = iEpTetData(corrOutEventInd);
    %         iEpTetDataType{3} = iEpTetData(mistOutEventInd);
    %         iEpTetDataType{4} = iEpTetData(outBEventInd);
    %         iEpTetDataType{5} = iEpTetData(inBEventInd);
    
    for ianimal = 1:length(F)
        % get event trials for this condition
        % cat stack the event snips into the 3rd dim
        
        nevents = size(ixcp.allNTDataCat{ianimal},3);
        nWave = length(timeWin);
        nNTrodes = 1; %size(intDataCat,1);
        nsamps = size(ixcp.allNTDataCat{ianimal},2);
        nData = nsamps*nevents;
        nConv = nWave+nData-1; % length of the result of convolution.
        nConv2pow = 2^nextpow2(nConv);
%         zpad2pow = (nConv2pow - nConv)/2;
        zpad2pow = nConv2pow - nConv;
%         if mod(zpad2pow,1) ~= 0
%             error('zpad2pow should be a whole number')
%         end
        for introde = 1:size(ixcp.allNTDataCat{ianimal},1)
            intDataCat = squeeze(ixcp.allNTDataCat{ianimal}(introde,:,:)); %squeeze will make nsampes x nevents.. 
            % FFT of data (doesn't change on frequency iteration)
            % dataX = fft(reshape(iEpTetDataCat(introde,:,:), 1,nData),nConv,2);
            % concat reshape all the events (peri-rip snips) for speed and
            % to reduce edge artifacts, then take the fft.
            % fft with the arg of nConv ensures that matlab will zero-pad
            % the half-wavelet duration on either side of the data series
            %             dataY = fft(reshape(iEpTetDataCat,nNTrodes,nData),nConv,2);
            dataY = fft(reshape(intDataCat,nNTrodes,nData),nConv2pow,2); %remember that reshape reshapes rowise so 
            % make sure that the each timeseries event is a column
            %             hz = linspace(0, 1500/2, floor(nConv/2)+1);
            %             plot(hz,2*abs(dataY(16,1:length(hz))/nData));
            
            %         %% fix ispc
            %         if calculateISPC
            %             % initialize output time-frequency data
            %             %                 ixpc.output{ianimal}{iday} = zeros(nchoosek(size(iEpTetDataCat,1),2), nsamps, num_frex);
            %             ixpc.phaseoutput{ianimal}{idatatype} = zeros(nchoosek(size(iEpTetDataCat,1),2), nsamps, num_frex);
            %         elseif calculateITPC
            %             % initialize output time-frequency data
            %             ixpc.phaseoutput{ianimal}{idatatype} = zeros(nNTrodes,nsamps, num_frex);
            %             %                 ixpc.output{ianimal}{iday} = zeros(nNTrodes,nsamps, num_frex);
            %             indices = [indices; F(ianimal).output{day}(eps(1)).index];
            %         end
            %%
            %         phdata = zeros(num_frex,nNTrodes,nsamps,nevents); %this crashes
            %         matlab because it exceeds maximum array size ... something like
            %         35GB...lol.. instead... i should probably loop over ntrodes and
            %         run the perm tests and everything then toss the trial data.
            %% loop over frequencies
            for fi=1:num_frex
                % create wavelet and get its FFT
                % creating the morlet wavelet by combingin the complex sine wave and the gaussian
                sine_wave = exp(2*1i*pi*frex(fi).*timeWin); %make a complex sine wave
                s = nWavecycles(fi)/(2*pi*frex(fi)); % std of the gaussian of the morlet wavelet. dependent on the freq and #cycles
                gaus_win = exp(-timeWin.^2./(2*s^2)); %the gaussian
                %                 wavelet  = exp(2*1i*pi*frex(fi).*timeWin) .* exp(-timeWin.^2./(2*s^2));
                wavelet = sine_wave .* gaus_win;
                waveletFFT = fft(wavelet,nConv2pow);
                %normalize wavelet to a maximum of 1. this will ensure that the units of convolution are the same as in the original data.
                waveletFFT = waveletFFT ./ max(waveletFFT);
                % run convolution (filtering) : Time-domain convolution in the frequency domain... because itz so much faster
                astmp = bsxfun(@times,dataY,waveletFFT); %multiply the power spectrum of the data and the wavelet
                astmp = ifft(astmp,nConv2pow,2); % take the inverse transform
                astmp = astmp(1,half_wave_size+1:end-half_wave_size-zpad2pow); %trim off the length of half the wavelet at the beginning and at the end
                as(:,:,fi) = reshape(astmp,nsamps,nevents); %reshape it back to ntrodes by nsamples by nevents
                phdata(:,:,fi) = angle(as(:,:,fi)); %get phase component of the analytic signal
                %             %% ISPC -- NEED TO FIX
                %             if calculateISPC % time series diff of all possible ntrode pairs.
                %                 %the idea here is to treat rows as pairs instead of individual ntrodes, as in ITPC
                %                 ispc = zeros(nchoosek(size(phdata{fi},1),2), nsamps, nevents);
%                                 m = logical(tril(ones(size(phdata{fi}(:,:,1),1)),-1)); %get indices of non-duplicates (below comb triangle)
%                                 [brow, bcol] = find(m); %get linear index of ntrode combination indices
%                                 indices = [F(ianimal).output{day}(eps(1)).index(bcol,:) F(ianimal).output{day}(eps(1)).index(brow,:)]; %convert to ntrodeID
                %                 for w = 1:size(phdata{fi},3); %loop through each event plane (ntrode x samples X event#)
                %                     iprm = permute(phdata{fi}(:,:,w),[3 2 1]); %transpose
                %                     B = bsxfun(@minus,phdata{fi}(:,:,w),iprm); %subtract across each NTrode-choose-two comb of phase time series
                %                     B = reshape(B,[],size(phdata{fi}(:,:,w),2)); %reshape back into 2D
                %                     B = B(m(:),:); %get rid of duplicates
                %                     ispc(:,:,w) = B; %save result
                %                 end
                %                 clear phasedata
                %                 phdata{fi} = ispc;
                %             end
                %             ixpc.index{ianimal}{fi} = indices;
                disp(sprintf('freq %d of %d',fi,num_frex));
                %             ixpc.frex{ianimal}{fi} = frex;
            
            %% AVERAGE ACROSS TRIALS (ITPC) OR SITES (ISPC)
%             for idatatype = 1:length(iEpTetIndsType)
                %             ixpc.phaseoutput{ianimal}{idatatype}(:,:,fi) = abs(mean(exp(1i*phdata{fi}(:,:,iEpTetIndsType{idatatype})),3));
                itpctmp = cellfun(@(x) abs(mean(exp(1i*phdata(:,x,fi)),2)), iEpTetIndsType, 'un', 0);
%                 itpc(:,:,fi) = cat(2, itpctmp{:});
                ixpc.phaseoutput{ianimal}{introde}(:,:,fi) = cat(2, itpctmp{:});
%                 phastmp{idatatype} = cellfun(@(x) abs(mean(exp(1i*x(:,:,iEpTetIndsType{idatatype})),3)), phdata, 'un', 0);
                %% POWER morlet wavelet power Zscored
                % power is the squared magnitude from origin to the location in complex space zscore using the entire window..
                % this is conservative as it includes the post-rip activity into the mean and std
                %             ixpc.poweroutput{ianimal}{idatatype}(:,:,fi) = mean(abs(as{fi}(:,:,iEpTetIndsType{idatatype})).^2,3);
                powrtmp = cellfun(@(x) mean(abs(as(:,x,fi)).^2,2), iEpTetIndsType, 'un', 0);
                ixpc.poweroutput{ianimal}{introde}(:,:,fi) = cat(2,powrtmp{:});
                %                 disp(sprintf('========== %s %s calculated=============',iEpTetDataTypeFields{idatatype}, calcfunction));
%             end
            %                         %% compute differential tf maps
            %                 ixpc.phaseoutput{ianimal}{6}(:,:,fi) = ixpc.phaseoutput{ianimal}{4}(:,:,fi) - ixpc.phaseoutput{ianimal}{5}(:,:,fi);
            %                 ixpc.poweroutput{ianimal}{6}(:,:,fi) = ixpc.poweroutput{ianimal}{4}(:,:,fi) - ixpc.poweroutput{ianimal}{5}(:,:,fi);
            %                 ixpc.phaseoutput{ianimal}{7}(:,:,fi) = ixpc.phaseoutput{ianimal}{2}(:,:,fi) - ixpc.phaseoutput{ianimal}{3}(:,:,fi);
            %                 ixpc.poweroutput{ianimal}{7}(:,:,fi) = ixpc.poweroutput{ianimal}{2}(:,:,fi) - ixpc.poweroutput{ianimal}{3}(:,:,fi);
            %             ixpc.waveletX{ianimal}{fi} = waveletFFT;
            
            
            end
            %% scratch
            %                 % power is the squared magnitude from origin to the location in complex space
            %                 temppow = mean(abs(as).^2,3);
            %                 %normalize the power to the pre-event period:
            % %                 baseidx   = dsearchn(EEG.times',[-400 -100]')
            %                 baseidx = [1 size(temppow,2)];
            %                 %currently using the mean of all trials at that time point for individual ntrodes.
            %                 %should i take the mean from each event instead? or the
            %                 %mean across every event and all times?
            %                 baseline = mean(temppow(:,baseidx(1):baseidx(2)),2);
            %                 basestd = std(temppow(:,baseidx(1):baseidx(2)),2);
            %                  % decibel normalized
            %                 ixpc.morletpowerDB{ianimal}{idatatype}(:,:,fi)  = 10*log10( bsxfun(@rdivide, temppow, baseline) );
            %                  % percent change normalized
            %                 ixpc.morletpowerPC{ianimal}{idatatype}(:,:,fi)  = 100 * bsxfun(@rdivide, bsxfun(@minus,temppow,baseline), baseline);
            
            %                 % alternative way to calculate power using hiLbert method:
            %                 %the hilbert only takes up to a 2D matrix so split up the
            %                 %trials
            %                 realas = real(as);
            %                 for itr = 1:size(realas,3)
            %                     itras = squeeze(realas(:,:,itr))'; % make the rows time as hilbert runs on columns for 2d input matrices
            %                     itrashilib(:,:,itr) = hilbert(itras)';
            %                 end
            %                 ixpc.hilbertpower{ianimal}{idatatype}(:,:,fi) = mean(abs(itrashilib).^2,3);
            %
            %                 % alternative way to calculate power using multitaper method:
            %
            
            %         end
            %         end
            %         return
            ixpc.datatypes{ianimal} = iEpTetDataTypeFields;
            %% compute differential tf maps
            ixpc.phaseoutput{ianimal}{introde}(:,6,:) = ixpc.phaseoutput{ianimal}{introde}(:,4,:) - ixpc.phaseoutput{ianimal}{introde}(:,5,:);
            ixpc.poweroutput{ianimal}{introde}(:,6,:) = ixpc.poweroutput{ianimal}{introde}(:,4,:) - ixpc.poweroutput{ianimal}{introde}(:,5,:);
            ixpc.phaseoutput{ianimal}{introde}(:,7,:) = ixpc.phaseoutput{ianimal}{introde}(:,2,:) - ixpc.phaseoutput{ianimal}{introde}(:,3,:);
            ixpc.poweroutput{ianimal}{introde}(:,7,:) = ixpc.poweroutput{ianimal}{introde}(:,2,:) - ixpc.poweroutput{ianimal}{introde}(:,3,:);
            %% permutation testing
            if runPermTest
%                 fields = {'phasemean_h0','phasestd_h0', 'phasezmap', 'powermean_h0', 'powerstd_h0', 'powerzmap','MCmax_clust_size', 'MCmin_max'};
                ixpc.phasemean_h0 = [];ixpc.phasestd_h0 = [];ixpc.phasezmap = [];
                ixpc.powermean_h0 = [];ixpc.powerstd_h0 = [];ixpc.powerzmap = [];
                
                clear permdata iprm permsetfull Apermset Bpermset phasepermOutput powerpermOutput
                disp('========== running perm test outbound-inbound=============')
                %% outbound vs inbound
                % create the shuffled distribution for each pixel
                permdataInds = [iEpTetIndsType{4}; iEpTetIndsType{5}];
                tic
                for iprm = 1:n_permutes;
                    disp(sprintf('nt%d outb-inb perm %d of %d',introde, iprm, n_permutes))
                    permsetfull = randperm(size(permdataInds,1));
                    Apermset = permsetfull(1:floor(length(permsetfull)/2));
                    Bpermset = permsetfull(end-floor((length(permsetfull)-1)/2):end);
                    phasepermOutput(:,iprm,:) = abs(mean(exp(1i*phdata(:,permdataInds(Apermset),:)),2)) - ...
                        abs(mean(exp(1i*phdata(:,permdataInds(Bpermset),:)),2));
                    powerpermOutput(:,iprm,:) = mean(abs(as(:,permdataInds(Apermset),:)).^2,2) - ...
                        mean(abs(as(:,permdataInds(Bpermset),:)).^2,2);
                   
                end
                toc
                % compute mean and standard deviation maps
                ixpc.phasemean_h0{introde}{6}(:,1,:) = squeeze(mean(phasepermOutput(:,:,:),2));
                ixpc.phasestd_h0{introde}{6}(:,1,:)  = squeeze(std(phasepermOutput(:,:,:),[],2));
                % Z-score the data
                ixpc.phasezmap{introde}{6}(:,1,:) = (ixpc.phaseoutput{ianimal}{introde}(:,6,:)-ixpc.phasemean_h0{introde}{6}) ./ ixpc.phasestd_h0{introde}{6};
                
                % compute mean and standard deviation maps
                ixpc.powermean_h0{introde}{6}(:,1,:) = squeeze(mean(powerpermOutput(:,:,:),2));
                ixpc.powerstd_h0{introde}{6}(:,1,:)  = squeeze(std(powerpermOutput(:,:,:),[],2));
                % Z-score the data
                ixpc.powerzmap{introde}{6}(:,1,:) = (ixpc.poweroutput{ianimal}{introde}(:,6,:)-ixpc.powermean_h0{introde}{6}) ./ ixpc.powerstd_h0{introde}{6};
%                 %% correct outbound vs mistake outbound
%                 disp('========== running perm test correctOut-mistakeOut=============')
%                 % create the shuffled distribution for each pixel
%                 clear permdata iprm permsetfull Apermset Bpermset phasepermOutput powerpermOutput
%                 permdataInds = [iEpTetIndsType{2}; iEpTetIndsType{3}];
%                 for iprm = 1:n_permutes;
%                     disp(sprintf('nt%d outb-inb perm %d of %d', introde, iprm, n_permutes))
%                     permsetfull = randperm(size(permdataInds,1));
%                     Apermset = permsetfull(1:floor(length(permsetfull)/2));
%                     Bpermset = permsetfull(end-floor((length(permsetfull)-1)/2):end);
%                     phasepermOutput(:,iprm,:) = abs(mean(exp(1i*phdata(:,permdataInds(Apermset),:)),2)) - ...
%                         abs(mean(exp(1i*phdata(:,permdataInds(Bpermset),:)),2));
%                     powerpermOutput(:,iprm,:) = mean(abs(as(:,permdataInds(Apermset),:)).^2,2) - ...
%                         mean(abs(as(:,permdataInds(Bpermset),:)).^2,2);
%                 end
%                 % compute mean and standard deviation maps
%                 ixpc.phasemean_h0{introde}{7}(:,1,:) = squeeze(mean(phasepermOutput(:,:,:),2));
%                 ixpc.phasestd_h0{introde}{7}(:,1,:)  = squeeze(std(phasepermOutput(:,:,:),[],2));
%                 % Z-score the data
%                 ixpc.phasezmap{introde}{7}(:,1,:) = (ixpc.phaseoutput{ianimal}{introde}(:,7,:)-ixpc.phasemean_h0{introde}{7}) ./ ixpc.phasestd_h0{introde}{7};
% 
%                 % compute mean and standard deviation maps
%                 ixpc.powermean_h0{introde}{7}(:,1,:) = squeeze(mean(powerpermOutput(:,:,:),2));
%                 ixpc.powerstd_h0{introde}{7}(:,1,:)  = squeeze(std(powerpermOutput(:,:,:),[],2));
%                 % Z-score the data
%                 ixpc.powerzmap{introde}{7}(:,1,:) = (ixpc.poweroutput{ianimal}{introde}(:,7,:)-ixpc.powermean_h0{introde}{7}) ./ ixpc.powerstd_h0{introde}{7};
%                 
                %% get info for multiple comparison calculation
                ixpc.MCmax_clust_size = [];ixpc.MCmin_max = [];
                for iprm = 1:n_permutes;
                    % take each permutation map, and transform to Z
                    threshimg = squeeze(phasepermOutput(:,iprm,:));
                    threshimg = (threshimg-squeeze(ixpc.phasemean_h0{introde}{6}(:,1,:)))./squeeze(ixpc.powerstd_h0{introde}{6}(:,1,:));
                    % threshold image at p-value
                    threshimg(abs(threshimg)<zval) = 0;
                    % find clusters (need image processing toolbox for this!)
                    islands = bwconncomp(threshimg);
                    if numel(islands.PixelIdxList)>0
                        % count sizes of clusters
                        tempclustsizes = cellfun(@length,islands.PixelIdxList);
                        % store size of biggest cluster
                        ixpc.MCmax_clust_size{ianimal}{introde} = [ixpc.MCmax_clust_size{ianimal}{introde}; max(tempclustsizes)];
                    end
                    % get extreme values (smallest and largest)
                    temp = sort( reshape(phasepermOutput(:,iprm,:),1,[]));
                    ixcp.MCmin_max{ianimal}{introde}(iprm,:) = [min(min(phasepermOutput(:,iprm,:))) max(max(phasepermOutput(:,iprm,:)))];
                end
            end
            
            
        end
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
%% ---------------- plot ITPC---------------------------------------------------------------------------------------------
%% ---------------- plot ITPC---------------------------------------------------------------------------------------------
if plotITPC
    %     clear F %save space on memory
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
%                 eval(sprintf('ixpc2plottmp = squeeze(ixpc.%soutput{ianimal}{iDT}(numsumSortInds(introde),:,:))'';', plotoutputtype));
                eval(sprintf('ixpc2plottmp = squeeze(ixpc.%soutput{ianimal}{}(numsumSortInds(introde):,iDT,:))'';', plotoutputtype));
%                 ixpc2plottmp = squeeze(ixpc.poweroutput{ianimal}{12}(:,6,:))';
                %                 end
                midpoint = (size(ixpc2plottmp,2)-1)/2; %get middle index of window
                ixpc2plottmp2 = ixpc2plottmp(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
                if strcmp(plotoutputtype, 'phase')
                    ixpc2plot = ixpc2plottmp2;
                elseif strcmp(plotoutputtype, 'power')
                    %                 ixpc2plottmp2 = zscore(ixpc2plottmp, [], 1);
                    %                 baseline = mean(ixpc2plottmp(:,1:midpoint+plotwin(1)*srate),2); %get baseline (start to middle index)
                    baselinemean = mean(ixpc2plottmp,2);
                    baselinestd = std(ixpc2plottmp,[],2);
                    %                 baselinemean = mean(ixpc2plottmp2,2);
                    %                 baselinestd = std(ixpc2plottmp2,[],2);
                    %                 ixpc2plot = 100 * bsxfun(@rdivide, bsxfun(@minus,ixpc2plottmp2,baselinemean), baselinemean); % normalize to percent change from baseline
                    ixpc2plot = bsxfun(@rdivide, bsxfun(@minus,ixpc2plottmp2,baselinemean), baselinestd); % normalize to percent change from baseline
                else
                    error('must be power or phase')
                end
%                intfig = figure
                contourf(intfig,plottimeWin,frex(1:end-1),ixpc2plot,num_frex,'linecolor','none'); %
                % threshold image at p-value, by setting subthreshold values to 0
                zmap(abs(zmap)<zval) = 0;
                % figure
                %                 imagesc(ixpc2plot);
                %                 patch([-.8, .8, .8, -.8], [-10 -10 80 80], iNTcolor, 'edgecolor','none')
                hold on;
                %                 if strcmp(plotoutputtype, 'phase')
                set(gca,'ydir','normal','xlim',[plotwin(1) plotwin(2)], 'ylim',[frex(1) frex(end)])
                %                     caxis([-10 10])
                caxis(clims);
                %                 else
                %                     set(gca,'ydir','normal','xlim',[-.5 .5])
                %                     calcfunction = 'power';
                %                 end
                %                 set(gca, 'YScale', 'log')
                colormap(usecolormap)
                if mod(introde, sfcols) ~= 1;
                    %                     set(gca, 'yscale','log','YTick', []);
                    %                     set(gca, 'yscale','log','ytick',ceil(logspace(log10(min(frex)),log10(max(frex)),6)))
                    set(gca, 'FontSize',8, 'FontName', 'Arial');%
                else
                    %                     set(gca, 'yscale','log','ytick',ceil(logspace(log10(min(frex)),log10(max(frex)),6)))
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
            %% ---- super title and colorbar----
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig1);
            xlab = 'Time (s)';
            ylab = 'Frequency (Hz)';
            supxlabel = text(.5, .02, xlab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'Parent', sprtitleax,...
                'Units', 'normalized', 'horizontalAlignment', 'center');
            supylabel = text(.01, .5, ylab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'rotation', 90, ...
                'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
            caxis(clims);
            %             caxis([-3 3])
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
            sprtit = sprintf('%s %s %s %s D%s %d-%dHz_rCyc%d-%d', plotoutputtype,ixpc.datatypes{ianimal}{iDT}, strjoin(epochEnvironment,'-'), animalID, strrep(num2str(days), '  ', '-'), min_freq,max_freq, range_cycles);
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

