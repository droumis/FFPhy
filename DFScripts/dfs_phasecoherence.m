
% to do:
% ITPC - plot open field
% ISPC - preserve the indices

%%% First create the processed data structure with dfs_riptriglfp.m

close all
calculateITPC = 0;
calculateISPC = 0;
calculateIXPC = calculateITPC || calculateISPC;
loadFilterOutput = calculateIXPC;
saveResultsOutput = calculateIXPC;
loadResultsOutput = 0;
plotITPC = 1;
plotISPC = 0;
savefigs = 1;
pausefigs = 0;
calcfunction = 'ITPC'; %ISPC or ITPC
behavestruct = 'BehaveState';
outputtype = 'phase'; % 'power'
%% ---------------- plotting params --------------------------
colorSet = 'DR1';
if strcmp(outputtype, 'phase')
    clims = [0 1]; %[0 .7]
else
    clims = [0 10];
    calcfunction = 'power';    
end
clims = [0 1]; %[0 .7]
usecolormap = 'jet';
win = [-.5 .5]; %in seconds. The morlet wavelet needs to be centered at time zero.. so always make this window symmetric about zero
srate = 1500;
%% ---------------- Data Filters --------------------------
animals = {'JZ1'};
% animals = {'JZ1', 'D13'};
days = [1:6];
filtfunction = 'riptriglfp';
LFPtypes = {'eeg', 'ripple'};%, 'theta', 'slowgamma', 'fastgamma'};
% LFPrangesHz = {'1-400', '140-250', '6-9', '20-50', '65 - 95'}; %need to make this a lookup from the filter mat
eventtype = 'rippleskons';
epochEnvironment = 'wtrack';% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
% epochType = 'run';
eventSourceArea = 'ca1';
% ripAreas = {'ca1', 'mec', 'por'};
% ntAreas = {'ca1', 'sub', 'mec', 'por', 'v2l', 'ref'};
% consensus_numtets = 1;   % minimum # of tets for consensus event detection
% minstdthresh = 5;        % STD. how big your ripples are
% exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
% minvelocity = 0;
% maxvelocity = 4;
% outputDirectory = '/typhoon/droumis/analysis';
%% ---------------- Paths and Title strings ---------------------------------------------------
investInfo = animaldef(lower('Demetris'));
figdirectory = sprintf('%s%s/', investInfo{4}, calcfunction);
filenamesave = sprintf('%s%s_%s_%s_%s', eventSourceArea, eventtype, epochEnvironment, cell2mat(animals), cell2mat(LFPtypes));
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
resultfilename = sprintf('%s_%s%s_%s_%s_D%s', calcfunction, eventSourceArea, eventtype, epochEnvironment, cell2mat(animals), strrep(num2str(days), '  ', '-'));
filenameTitle = strrep(resultfilename,'_', ' ');
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s%s/%s',investInfo{2},filtfunction, filename))
end
%% Calculate inter-trial phase clustering based on mike cohen's book

% wavelet parameters
min_freq =  2;
max_freq = 30;
num_frex = (max_freq - min_freq); %frequency resolution default 1
% set range for variable number of wavelet cycles
% range_cycles = [min_freq*win(2) max_freq*win(2)]; % is there a principled way to do this?
% more wavelets for the same frequency will reduce the temporal precision.
% more wavelets for the same frequency will increase the frequency precision
range_cycles = [2 10];
% win = [[1,-1*win(1)]; %seconds
timeWin = win(1):1/srate:win(2);
% frex = logspace(log10(min_freq),log10(max_freq),num_frex);
frex = linspace(min_freq,max_freq,num_frex+1);
nWavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex); %number of wavelet cycles per freq
% nWavecycles = repmat(range_cycles,1,num_frex);
% figure
% plot(nWavecycles,frex, 'o-')
% ylabel('frequency')
% xlabel('# cycles')
half_wave_size = (length(timeWin)-1)/2;
if calculateIXPC
    tic
    clear ixpc
    ixpc.type = calcfunction;
    ixpc.animals = [];
    ixpc.range_cycles = range_cycles;
    ixpc.max_freq = max_freq;
    ixpc.min_freq = min_freq;
    ixpc.num_frex = num_frex;
    ixpc.win = win;
    ixpc.srate = srate;
    
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
        %get the animal state for every event
        for iday = 1:length(days)
            day = days(iday);
            eps = find(~cellfun(@isempty,{F(ianimal).output{day}.index})); %get nonempty eps
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
        % split the events out into groups based on correct out / mistake out
        corrOutEventInd = find(eventstate.state(:,6) == 1);
        mistOutEventInd = find(eventstate.state(:,7) == 1);
        outBEventInd = [corrOutEventInd; mistOutEventInd];
        inBEventInd = setdiff([1:length(eventstate.state(:,1))], outBEventInd);
        iEpTetDataType{1} = iEpTetData;
        iEpTetDataType{2} = iEpTetData(corrOutEventInd);
        iEpTetDataType{3} = iEpTetData(mistOutEventInd);
        iEpTetDataType{4} = iEpTetData(outBEventInd);
        iEpTetDataType{5} = iEpTetData(inBEventInd);
        iEpTetDataTypeFields = {'all', 'corrOut', 'mistOut', 'outB', 'inB'};
        % to do: split the events out into groups based on trajectory (CR RC CL LC)
        for idatatype = 1:length(iEpTetDataType)
            iEpTetData = iEpTetDataType{idatatype};
            
            iEpTetDataCat = cat(3, iEpTetData{:}); % cat stack the rip snips into the 3rd dim
            % FFT parameters
            nWave = length(timeWin);
            nNTrodes = length(iEpTetDataCat(:,1,1));
            nsamps = length(iEpTetDataCat(1,:,1));
            nevents = length(iEpTetDataCat(1,1,:));
            nData = nsamps*nevents;
            nConv = nWave+nData-1; % length of the result of convolution. 
            
            % FFT of data (doesn't change on frequency iteration)
            % dataX = fft(reshape(iEpTetDataCat(introde,:,:), 1,nData),nConv,2);
            % concat reshape all the events (peri-rip snips) for speed and
            % to reduce edge artifacts, then take the fft.
            % fft with the arg of nConv ensures that matlab will zero-pad
            % the half-wavelet duration on either side of the data series
            dataY = fft(reshape(iEpTetDataCat,nNTrodes,nData),nConv,2);
%             hz = linspace(0, 1500/2, floor(nConv/2)+1);
%             plot(hz,2*abs(dataY(16,1:length(hz))/nData));
            indices = [];
            if calculateISPC
                % initialize output time-frequency data
                %                 ixpc.output{ianimal}{iday} = zeros(nchoosek(size(iEpTetDataCat,1),2), nsamps, num_frex);
                ixpc.phaseoutput{ianimal}{idatatype}= zeros(nchoosek(size(iEpTetDataCat,1),2), nsamps, num_frex);
            elseif calculateITPC
                % initialize output time-frequency data
                ixpc.phaseoutput{ianimal}{idatatype} = zeros(nNTrodes,nsamps, num_frex);
                %                 ixpc.output{ianimal}{iday} = zeros(nNTrodes,nsamps, num_frex);
                indices = [indices; F(ianimal).output{day}(eps(1)).index];
            end
            
            % loop over frequencies
            for fi=1:num_frex
                % create wavelet and get its FFT
                % creating the morlet wavelet by combingin the complex sine wave and the gaussian
                sine_wave = exp(2*1i*pi*frex(fi).*timeWin); %make a complex sine wave
                s = nWavecycles(fi)/(2*pi*frex(fi)); % std of the gaussian of the morlet wavelet. dependent on the freq and #cycles
                gaus_win = exp(-timeWin.^2./(2*s^2)); %the gaussian
%                 wavelet  = exp(2*1i*pi*frex(fi).*timeWin) .* exp(-timeWin.^2./(2*s^2));
                wavelet = sine_wave .* gaus_win;
                waveletX = fft(wavelet,nConv);
                %normalize wavelet to a maximum of 1. this will ensure that the units of convolution are the same as in the original data.
                waveletX = waveletX ./ max(waveletX);
                % run convolution (filtering) : Time-domain convolution in the frequency domain... because itz so much faster 
                as = bsxfun(@times,dataY,waveletX); %multiply the power spectrum of the data and the wavelet
                as = ifft(as,nConv,2); % take the inverse transform
                as = as(:,half_wave_size+1:end-half_wave_size); %trim off the length of half the wavelet at the beginning and at the end
                as = reshape(as,nNTrodes,nsamps,nevents); %reshape it back to ntrodes by nsamples by nevents
                phdata = angle(as); %get phase component of the analytic signal
                if calculateISPC % time series diff of all possible ntrode pairs.
                    %the idea here is to treat rows as pairs instead of individual ntrodes, as in ITPC
                    ispc = zeros(nchoosek(size(phdata,1),2), nsamps, nevents);
                    m = logical(tril(ones(size(phdata(:,:,1),1)),-1)); %get indices of non-duplicates (below comb triangle)
                    [brow, bcol] = find(m); %get linear index of ntrode combination indices
                    indices = [F(ianimal).output{day}(eps(1)).index(bcol,:) F(ianimal).output{day}(eps(1)).index(brow,:)]; %convert to ntrodeID
                    for w = 1:size(phdata,3); %loop through each event plane (ntrode x samples X event#)
                        iprm = permute(phdata(:,:,w),[3 2 1]); %transpose
                        B = bsxfun(@minus,phdata(:,:,w),iprm); %subtract across each NTrode-choose-two comb of phase time series
                        B = reshape(B,[],size(phdata(:,:,w),2)); %reshape back into 2D
                        B = B(m(:),:); %get rid of duplicates
                        ispc(:,:,w) = B; %save result
                    end
                    clear phasedata
                    phdata = ispc;
                end
                % compute IXPC
                ixpc.phaseoutput{ianimal}{idatatype}(:,:,fi) = abs(mean(exp(1i*phdata),3));
                % power is the squared magnitude from origin to the
                % dot-product location in complex space
                temppow = mean(abs(as).^2,3);
                %normalize the power to the pre-event period:
%                 baseidx   = dsearchn(EEG.times',[-400 -100]')
                baseidx = [1 size(temppow,2)]; 
                %currently using the mean of all trials at that time point for individual ntrodes.
                %should i take the mean from each event instead? or the
                %mean across every event and all times?
                ixpc.poweroutput{ianimal}{idatatype}(:,:,fi) = 10*log10( temppow./mean(temppow(:,baseidx(1):baseidx(2,:)),3));
                
                
                %alternative way to calculate power using himbert method: 
                ixpc.hilbertpower{ianimal}{idatatype}(:,:,fi) = mean(abs(hilbert(real(as))).^2,3);
                %                 ixpc.waveletX{ianimal}{iday} = waveletX;
                %                 ixpc.index{ianimal}{iday} = indices;
                ixpc.waveletX{ianimal}{fi} = waveletX;
                ixpc.index{ianimal}{fi} = indices;
                disp(sprintf('freq %d of %d',fi,num_frex));
            end
            ixpc.frex{ianimal} = frex;
            %             disp(sprintf('==========day %d calculated=============',iday));
            disp(sprintf('========== %s %s calculated=============',iEpTetDataTypeFields{idatatype}, calcfunction));
            %         end
        end
        ixpc.datatypes{ianimal} = iEpTetDataTypeFields;
    end
    toc
end
%% ---------------- Save RESULTS Output ---------------------------------------------------
if saveResultsOutput == 1;
    if ~isdir(sprintf('%s%s/', investInfo{3}, calcfunction));
        mkdir(sprintf('%s%s/', investInfo{3}, calcfunction));
    end
    save(sprintf('%s%s/%s_%dHz-%dHz_rCyc%d-%d', investInfo{3}, calcfunction, resultfilename, min_freq, max_freq,range_cycles), 'ixpc','-v7.3');
    disp(sprintf('%s saved', resultfilename))
end
%% ---------------- Load Results Output ---------------------------------------------------
if loadResultsOutput == 1;
    load(sprintf('%s%s/%s', investInfo{3}, calcfunction, resultfilename));
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
        ntrodesIndices = ixpc.index{ianimal}{1};
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
        for iDT = 7;%1:(length(ixpc.datatypes{ianimal}))+2
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
                if iDT == 6
                    eval(sprintf('ixpc2plot = squeeze(ixpc.%soutput{ianimal}{4}(numsumSortInds(introde),:,:))'' - squeeze(ixpc.%soutput{ianimal}{5}(numsumSortInds(introde),:,:))'';', outputtype, outputtype));
                    
%                     ixpc2plot = (squeeze(ixpc.phaseoutput{ianimal}{4}(numsumSortInds(introde),:,:))' - squeeze(ixpc.phaseoutput{ianimal}{5}(numsumSortInds(introde),:,:))');
                elseif iDT == 7
                    eval(sprintf('ixpc2plot = squeeze(ixpc.%soutput{ianimal}{2}(numsumSortInds(introde),:,:))'' - squeeze(ixpc.%soutput{ianimal}{3}(numsumSortInds(introde),:,:))'';', outputtype, outputtype));
                else
                    eval(sprintf('ixpc2plot = squeeze(ixpc.%soutput{ianimal}{iDT}(numsumSortInds(introde),:,:))'';', outputtype));
                end
                contourf(intfig,timeWin,frex(1:end-1),ixpc2plot,num_frex,'linecolor','none');%, )
% figure
%                 imagesc(ixpc2plot);
                %                 patch([-.8, .8, .8, -.8], [-10 -10 80 80], iNTcolor, 'edgecolor','none')
                hold on;
                if strcmp(outputtype, 'phase')
                    set(gca,'clim',clims,'ydir','normal','xlim',[-.5 .5])
                else
                    set(gca,'ydir','normal','xlim',[-.5 .5])
                    calcfunction = 'power';
                end
%                 set(gca, 'YScale', 'log')
                
                if mod(introde, sfcols) ~= 1;
                    set(gca, 'YTick', []);
                else
                    set(gca, 'YTick',[round(linspace(min(frex),max(frex),6))],'FontSize',8, 'FontName', 'Arial');%
                end
                if introde <= (sfrows-1)*sfcols;
                    set(gca, 'XTick', []);
                else
                    set(gca, 'XTick',[win(1):win(2)/2:win(2)],'FontSize',8, 'FontName', 'Arial');
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
            if iDT == 6
                sprtit = sprintf('%s %s %s %s D%s %d-%dHz_rCyc%d-%d', calcfunction,'outB-inB', epochEnvironment, animalID, strrep(num2str(days), '  ', '-'), min_freq,max_freq, range_cycles);
            elseif iDT == 7
                sprtit = sprintf('%s %s %s %s D%s %d-%dHz_rCyc%d-%d', calcfunction,'outBCorr-outBMist', epochEnvironment, animalID, strrep(num2str(days), '  ', '-'), min_freq,max_freq, range_cycles);
            else
            sprtit = sprintf('%s %s %s %s D%s %d-%dHz_rCyc%d-%d', calcfunction,ixpc.datatypes{ianimal}{iDT}, epochEnvironment, animalID, strrep(num2str(days), '  ', '-'), min_freq,max_freq, range_cycles);
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
            clrbartit = text(posx(1)+posx(3)/2, posx(2)-MarginBottom/2, calcfunction, 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial','horizontalAlignment', 'center');
            
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

