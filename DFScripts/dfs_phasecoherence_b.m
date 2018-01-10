


%% Dashboard
%%% First create the processed data structure similar to dfs_riptriglfp.m
close all
runFilterFramework = 0;
; saveFilterOutput = runFilterFramework;
; loadFilterOutput = 0;

%%% Then run phasecoherence on the LFP traces

calculateITPC = runFilterFramework;
calculateISPC = 0;
calculateIXPC = calculateITPC || calculateISPC;
calcEventState = calculateIXPC;
saveEventState = calcEventState;
saveResultsOutput = calculateIXPC;
runPermTest = 0; %calculateIXPC;
combineByArea = 0;
loadResultsOutput = 0;
plotITPC = 1;
plotPower = 0;
plotISPC = 0;
savefigs = 1;
pausefigs = 0;
calcfunction = 'ITPC'; %ITPC
behavestruct = 'BehaveState';
plotoutputtype = 'ITPC'; % 'power' or 'ITPC'
plotzmask = 0;
%% ----------------  params --------------------------
colorSet = 'DR1';
figspecs = 'itpc';
ripSet = 'DR4';
waveSet = 'DR1';
indsSet = 'DR1';
% if strcmp(plotoutputtype, 'phase')
%     clims = [0 .7]; % mean ITPC goes from 0-1
%     usecolormap = 'jet'; %cbrewer('seq', 'Greens', 255, 'linear');
%     usecolormap = buildcmap('kryww');
% else
% clims = [-8 8]; %zscore power
coloraxis = [-10 10];
%     usecolormap = 'jet';
usecolormap = flipud(cbrewer('div', 'RdBu', 255, 'linear')); %red high white neutral blue low
%         usecolormap = flipud(usecolormap);
% %     calcfunction = 'power';
% end
% clims = [0 1]; %[0 .7]
% usecolormap = 'jet';

% usecolormap = buildcmap('rwb');
srate = 1500;
win = [-1.5 1.5]; %in seconds. The morlet wavelet needs to be centered at time zero.. so always make this window symmetric about zero
plotwin = [-.5 .5];
basewin = plotwin(2); %baseline window to use for 1/f normalization
plottimeWin = plotwin(1):1/srate:plotwin(2);
%baseline indices. length of 1/2 plot window, centered on first plot ind
% i.e. 0.5 sec basewindow starting at .25 sec before the start of the
% plotwindow
baseind = ([(win(2)-plotwin(2))-basewin/2 (win(2)-plotwin(2))+basewin/2]).*srate;
% baseind = [1 abs(win(1))*srate-(abs(plotwin(1))*srate)]; %the entire duration before the plotwin
% basedur = abs(win(1))+abs(win(2));
% baselineind = ([1 (srate*basedur)+1]);
% p-value
pval = 0.05;
% convert p-value to Z value
zval = abs(norminv(pval));
% number of permutations
n_permutes = 10; %1000 takes 24 hours for JZ1 wtrack days 1:6
% OvsIind = 6;
% outboundInd = 4;
% inboundInd = 5;
% CvsMind = 7;
% correctOutInd = 2;
% mistakeOutInd = 3;
%% ----------------- wavelet parameters ------------------------------------
min_freq =  4;
max_freq = 30;
frexres = 1; %frequency resolution default 4
numfrex = floor((max_freq - min_freq)/frexres);
% set range for variable number of wavelet cycles
% range_cycles = [min_freq*win(2) max_freq*win(2)]; % is there a principled way to do this?
% more wavelets for the same frequency will reduce the temporal precision.
% more wavelets for the same frequency will increase the frequency precision
% range_cycles = [4 floor(max_freq/10)];
range_cycles = (([min_freq max_freq])/min_freq)*2; % sliding window
% range_cycles = [2];
% win = [[1,-1*win(1)]; %seconds
timeWin = win(1):1/srate:win(2);
frex = logspace(log10(min_freq),log10(max_freq),numfrex);
% frex = linspace(min_freq,max_freq,numfrex+1);
% nWavecycles = linspace(range_cycles(1),range_cycles(end),numfrex); %number of wavelet cycles per freq
nWavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),numfrex); %number of wavelet cycles per freq
% wp.nWavecycles = repmat(range_cycles,1,num_frex);
% figure
% plot(wp.nWavecycles,frex, 'o-')
% ylabel('frequency')
% xlabel('# cycles')
half_wave_size = (length(timeWin)-1)/2;
%% ---------------- Data Filters --------------------------
% animals = {'D12', 'JZ1', 'JZ2'};
animals = {'JZ2'};
days = [];
filtfunction = 'riptriglfp';
LFPtypes = {'eeg'};%, 'ripple', 'theta', 'lowgamma', 'fastgamma'}; %
eventtype = 'rippleskons';
epochType = {'run'};
epochEnvironment = {'wtrack'};% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
eventSourceArea = 'ca1';
ntAreas = {'ca1','sub','mec', 'por', 'v2l', 'ref'}; %which areas to include in datafilrer
rippms = getRipParams(ripSet); % set parameters for ripple detection like speed, numtetrodes, etc.
%% ---------------- Paths and Title strings ---------------------------------------------------
investInfo = animaldef(lower('Demetris'));
filtOutputDirectory = sprintf('%s%s/', investInfo{2}, filtfunction);
figdirectory = sprintf('%s%s/%s/', investInfo{4}, calcfunction, plotoutputtype);
% filenamesave = sprintf('%s%s_%s_%s_%s', eventSourceArea, eventtype, epochEnvironment, cell2mat(animals), cell2mat(LFPtypes));
filenamesave = sprintf('%s%sSTD%d_%s_%s_%s', eventSourceArea, eventtype, rippms.minstdthresh, strjoin(epochEnvironment,'-'),...
    strjoin(animals,'-'), strjoin(LFPtypes, '-'));%strjoin(arrayfun(@(x) num2str(x,'%-2d'),days,'UniformOutput',false),'-'));
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
resultfilename = filenamesave; %sprintf('%s_%s%s_%s_%s_D%s', calcfunction, eventSourceArea, eventtype, strjoin(epochEnvironment,'-'), cell2mat(animals), strrep(num2str(days, '%-2d'),' ', '-'));
filenameTitle = strrep(resultfilename,'_', ' ');
iEpTetDataTypeFields = {'all', 'corrOut', 'mistOut', 'outB', 'inB','outB-inB', 'corrOut-mistOut'};

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
    eval(['F = setfilterfunction(F, [''dfa_'' filtfunction], {[eventSourceArea eventtype],' strjoin(arrayfun(@(x) sprintf('''%s'',', cell2mat(x)), LFPtypes,'UniformOutput',false)) '},''eventtype'', [eventSourceArea eventtype], ''LFPtypes'', LFPtypes, ''win'', win);']);
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
%% ---------------- Load Results Output ---------------------------------------------------
if loadResultsOutput == 1;
    load(sprintf('%s%s/%s_%dHz-%dHz_rCyc%d-%d.mat', investInfo{3}, calcfunction, resultfilename, min_freq, max_freq,range_cycles(1), range_cycles(2)));
end
%% ---------------- Calc Event State ---------------------------------------------------
if calcEventState
    [ixpc] = calcEventPerformanceState(F, animals, behavestruct);
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
    ixpc.phOut = [];    ixpc.powerOut = [];
    ixpc.type = calcfunction; ixpc.range_cycles = range_cycles; ixpc.max_freq = max_freq; ixpc.min_freq = min_freq;
    ixpc.num_frex = numfrex; ixpc.win = win; ixpc.srate = srate;
    for ian = 1:length(F)
        ixpc.index{ian} = [];
        andays = find(~cellfun(@isempty,F(ian).output)); %get nonempty eps
        iEpTetData = [];
        for iday = 1:length(andays)
            day = andays(iday);
            eps = find(~cellfun(@isempty,{F(ian).output{day}.index})); %get nonempty eps
            ixpc.index{ian} = [ixpc.index{ian}; F(ian).output{day}(eps(1)).index];
            for iep = eps;
                iEpTetData = cat(2,iEpTetData,F(ian).output{day}(iep).data{1,1}); %cat the event snips across eps, days
            end
        end
        % trim events too close to epoch boundaries
        if size(ixpc.eventstate{ian}.state(:,1),1) ~= size(iEpTetData,2)
            iEpTetData = iEpTetData(logical(ixpc.removevec{ian}));
            if size(ixpc.eventstate{ian}.state(:,1),1) ~= size(iEpTetData,2)
                error('mismatch between number of state info and lfp sanples')
            end
        end
        ixpc.allNTDataCat{ian} = cat(3, iEpTetData{:});
%         % split the events out into groups based on correct out / mistake out
%         corrOutBInd = find(cell2mat(cellfun(@(x) strcmp(x,'corrOut'), strsplit(iEpTetDataTypeFields, ' '), 'UniformOutput', false)));
%         mistOutBInd = find(cell2mat(cellfun(@(x) strcmp(x,'mistOut'), strsplit(iEpTetDataTypeFields, ' '), 'UniformOutput', false)));
%         corrOutEventInd = find(ixpc.eventstate{ian}.state(1:end,corrOutBInd) == 1);
%         mistOutEventInd = find(ixpc.eventstate{ian}.state(1:end,mistOutBInd) == 1);
%         outBEventInd = [corrOutEventInd; mistOutEventInd];
%         inBEventInd = setdiff([1:length(ixpc.eventstate{ian}.state(1:end,1))], outBEventInd)';
%         iEpTetIndsType = {[1:length(ixpc.eventstate{ian}.state(1:end,1))], corrOutEventInd, mistOutEventInd, outBEventInd, inBEventInd};
        iEpTetIndsType = getIndsByType(indsSet,ixpc, ian, iEpTetDataTypeFields);
        wp = getWaveParams(waveSet, ixpc, ian, timeWin);
%         wp.nevents = size(ixpc.allNTDataCat{ian},3);
%         wp.nWave = length(timeWin);
%         wp.nTimeSeries = 1; %keep as 1 right now even if theres a bunch of tetrodes.. used to be called 'wp.nNTrodes'
%         wp.nsamps = size(ixpc.allNTDataCat{ian},2);
%         wp.nData = wp.nsamps*wp.nevents;
%         wp.nConv = wp.nWave+wp.nData-1; % length of the result of convolution.
%         wp.nConv2pow = 2^nextpow2(wp.nConv); %find the next power of 2 (for FFT optimization purposes)
%         wp.zpad2pow = wp.nConv2pow - wp.nConv;
%         wp.nNTrodes = size(ixpc.allNTDataCat{ian},1);
        
        for int = 1:wp.nNTrodes
            ixpc.frex{ian}{int} = frex;
%             intDataCat = squeeze(ixpc.allNTDataCat{ian}(int,:,:)); %squeeze will make nsampes x wp.nevents..
            intDataCat = squeeze(ixpc.allNTDataCat{ian}(int,:,:)); %squeeze will make nsampes x wp.nevents..
            % concat reshape all the events (peri-rip snips) for speed and to reduce edge artifacts, then take the fft.
            % fft with the arg of wp.nConv ensures that matlab will zero-pad the half-wavelet duration on either side of the data series
%                         dataY = fft(reshape(iEpTetDataCat,wp.nNTrodes,wp.nData),wp.nConv,2);
            dataY = fft(reshape(intDataCat,wp.nTimeSeries,wp.nData),wp.nConv2pow,2); % reshape reshapes row-wise into a 1 x d vector
%             dataY = fft(reshape(ixpc.allNTDataCat{ian},wp.nNTrodes,wp.nData),wp.nConv2pow,2); % reshape reshapes row-wise into a 1 x d vector
            %             dataY{introde} = fft(reshape(squeeze(ixpc.allNTDataCat{ianimal}(introde,:,:)),wp.nNTrodes,wp.nData),wp.nConv2pow,2); % reshape reshapes row-wise into a 1 x d vector
            %         %% fix ispc
            %         if calculateISPC
            %             % initialize output time-frequency data
            %             %                 ixpc.output{ianimal}{iday} = zeros(nchoosek(size(iEpTetDataCat,1),2), wp.nsamps, num_frex);
            %             ixpc.phOut{ianimal}{idatatype} = zeros(nchoosek(size(iEpTetDataCat,1),2), wp.nsamps, num_frex);
            [as{int}, ph{int}] = computeAnalyticSignal(dataY, wp.nConv2pow, timeWin, frex, numfrex, nWavecycles, wp.zpad2pow, half_wave_size, wp.nsamps, wp.nevents);
            ixpc = computeITPC(ixpc,ph,ian,int,iEpTetIndsType,baseind,numfrex);
            ixpc = computePower(ixpc,as,ian,int,iEpTetIndsType,baseind,numfrex);
            %% permutation testing... this will take TIME and a shit ton of RAM.. go get covfefe
            if runPermTest
                permTestType = 'outboundVSinbound';
                ixpc = permutationTest(ixpc, as, ph, ian, int, iEpTetIndsType, iEpTetDataTypeFields, n_permutes, permTestType);
                permTestType = 'correctOutvsmistakeOut';
                ixpc = permutationTest(ixpc, as, ph, ian, int, iEpTetIndsType, iEpTetDataTypeFields, n_permutes, permTestType);
            end
        end
        ixpc.datatypes{ian} = [];
        ixpc.datatypes{ian} = iEpTetDataTypeFields;
        %% combine by area
        if combineByArea
            [nt] = getNTinfo(ixpc.index, animals, ian);
            areas2proc = unique(nt.strAreas);
                for iar = 1:length(areas2proc)
                    areaID = areas2proc{iar};
                    areainds = find(strcmp(areaID, nt.strAreas));
                    %itpc
                    areaITPCoutputtmp = arrayfun(@(x) ixpc.logbasedITPCOut{ian}{x}(:,6,:), areainds, 'un', 0);
                    ixpc.areaITPCmean{ian}{iar} = mean(cat(2,areaITPoutputtmp{:}), 2);
                    areaITPCzmasktmp = arrayfun(@(x) ixpc.ITPCzmask{ian}{x}{6}, areainds, 'un', 0);
                    ixpc.areaITPCzmaskmean{ian}{iar} = mean(cat(2,areaITPCzmasktmp{:}), 2);
                    % power
                    areapoweroutputtmp = arrayfun(@(x) ixpc.logbasedpowerOut{ian}{x}(:,6,:), areainds, 'un', 0);
                    ixpc.areapowermean{ian}{iar} = mean(cat(2,areapoweroutputtmp{:}), 2);
                    areapowerzmasktmp = arrayfun(@(x) ixpc.powerzmask{ian}{x}{6}, areainds, 'un', 0);
                    ixpc.areapowerzmaskmean{ian}{iar} = mean(cat(2,areapowerzmasktmp{:}), 2);
                end
                ixpc.areas{ian} = areas2proc;
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
    disp(sprintf('SAVED RESULTS ++++++++++ %s%s/%s_%dHz-%dHz_rCyc%d-%d.mat', investInfo{3}, calcfunction, resultfilename, min_freq, max_freq,range_cycles(1), range_cycles(2)))
end

%% ---------------- plot ITPC---------------------------------------------------------------------------------------------
if plotITPC
    warning('off', 'MATLAB:contour:ConstantData')
    fig = getFigspecs(figspecs);
    for ian = 1:length(ixpc.animals)
        [nt] = getNTinfo(ixpc.index, animals, ian, colorSet);
        animalID = ixpc.animals{ian};
        for iDT = 1:(length(ixpc.datatypes{ian}))%+2
            [ifig, sfrows, sfcols] = prepFig(savefigs, pausefigs, nt, fig);
            %% ---- loop across all tets for this day ----
            for int = 1:nt.nNTrodes;
                introdeID = nt.ntrodesIndices(nt.numsumSortInds(int),3);
                isupareatag = nt.ntrodeTags{nt.numsumSortInds(int)}.suparea;
                iareatag = nt.ntrodeTags{nt.numsumSortInds(int)}.area;
                isubareatag = nt.ntrodeTags{nt.numsumSortInds(int)}.subarea;
                iNTcolor = nt.icolors(int,:);
                intfig = subaxis(sfrows,sfcols,int, 'SpacingVertical', fig.SpacingVertical, 'SpacingHorizontal', fig.SpacingHorizontal, ...
                    'Padding', fig.Padding, 'MarginLeft', fig.MarginLeft, 'MarginRight', fig.MarginRight, 'MarginTop', fig.MarginTop, ...
                    'MarginBottom', fig.MarginBottom);
                %% ---- Main Plotting ----
                eval(sprintf('idata2plot = squeeze(ixpc.logbased%sout{ian}{nt.numsumSortInds(int)}(:,iDT,:))'';', plotoutputtype));        
                idata2plot = trim2win(idata2plot, srate, plotwin);
                try
                    [~,bn] = contourf(intfig, plottimeWin,frex(1:end-1),idata2plot,fig.contourRes,'linecolor','none'); %
                catch
                    [~,bn] = contourf(intfig, plottimeWin,frex,idata2plot,fig.contourRes,'linecolor','none'); %
                end
                caxis(coloraxis)
%                 caxis(coloraxis);
                hold on;
                if (iDT == 6 || iDT == 7) && plotzmask % plot the significance countours
                    eval(sprintf('zmask2plot = squeeze(ixpc.%szmask{ian}{nt.numsumSortInds(int)}{iDT})'';', plotoutputtype));
                    eval(sprintf('MCminmax = abs(ixpc.MC_%s_minmax{ian}{nt.numsumSortInds(int)}{iDT})'';',plotoutputtype));
                    eval(sprintf('irawdata = squeeze(ixpc.%sout{ian}{nt.numsumSortInds(int)}(:,iDT,:))'';', plotoutputtype));
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
                    [~,mc] = contour(intfig,plottimeWin,frex(1:end-1),logical(MCthreshmap),1);
%                     mc.LineColor = [.6 .6 .6];
                    mc.LineColor = 'm';%[.85 .85 .85];
                else
                    
                end
                
                
%                 if 0 % non baseline normalized
%                     %                     eval(sprintf('ixpc2plotFULL = squeeze(ixpc.area%smean{ianimal}{introde})'';', plotoutputtype));
%                     eval(sprintf('ixpc2plotFULL = squeeze(ixpc.area%szmaskmean{ianimal}{introde})'';', plotoutputtype));
%                     midpoint = (size(ixpc2plotFULL,2)-1)/2; %get middle index of window
%                     ixpc2plot = ixpc2plotFULL(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
%                     [~,bn] = contourf(intfig,plottimeWin,frex(1:end-1),ixpc2plot,100,'linecolor','none'); %
%                     hold on;
%                     %                     eval(sprintf('zmaskFULL = squeeze(ixpc.area%szmaskmean{ianimal}{introde})'';', plotoutputtype));
%                     % %                     zmaskFULL = squeeze(ixpc.areaphasezmaskmean{ianimal}{introde})';
%                     %                     midpoint = (size(ixpc2plotFULL,2)-1)/2; %get middle index of window
%                     %                     zmask2plot = zmaskFULL(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
%                     zmask = ixpc2plot;
%                     zmask(abs(ixpc2plot)<zval) = 0;
%                     [~,h] = contour(intfig,plottimeWin,frex(1:end-1),logical(zmask),1);
%                     h.LineColor = [.85 .85 .85];
%                 else %baseline normalized
                    %                     eval(sprintf('ixpc2plotFULL = squeeze(ixpc.%soutput{ianimal}{numsumSortInds(introde)}(:,iDT,:))'';', plotoutputtype));
                    %                     midpoint = (size(ixpc2plotFULL,2)-1)/2; %get middle index of window
                    %                     ixpc2plot = ixpc2plotFULL(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
                    %                     [~,bn] = contourf(intfig,plottimeWin,frex(1:end-1),ixpc2plot,100,'linecolor','none'); %
% %                     eval(sprintf('BaseNormixpc2plotFULL = squeeze(ixpc.base%sOut{ianimal}{numsumSortInds(introde)}(:,iDT,:))'';', plotoutputtype));
%                     eval(sprintf('basedixpc = squeeze(ixpc.based%sOut{ian}{nt.numsumSortInds(int)}(:,iDT,:))'';', 'ph'));
%                     midpoint = (size(basedixpc,2)-1)/2; %get middle index of window
%                     basedixpc = basedixpc(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
% %                     ixpc2plot = (basedixpc - min(min(basedixpc))) ./ max(max(basedixpc));
%                     [~,bn] = contourf(intfig,plottimeWin,frex(1:end-1),ixpc2plot,100,'linecolor','none'); %
% %                     set(gca, 'yscale','log')
%                 end
%                 hold on;
%                 
%                 if iDT == 6 || iDT == 7 % plot the significance countours
%                     plotoutputtype = 'phase';
%                     eval(sprintf('ixpc2plotFULL = squeeze(ixpc.based%sOut{ian}{nt.numsumSortInds(int)}(:,iDT,:))'';', 'ph'));
%                     eval(sprintf('Zmaskixpc2plotFULL = squeeze(ixpc.%szmask{ian}{nt.numsumSortInds(int)}{iDT})'';', plotoutputtype));
%                     eval(sprintf('MCminmax = abs(ixpc.MC%s_minmax{ian}{nt.numsumSortInds(int)}{iDT})'';',plotoutputtype));
%                     ixpc2plot = ixpc2plotFULL(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
%                     Zmaskixpc2plot = Zmaskixpc2plotFULL(:,midpoint + plotwin(1)*srate : midpoint + plotwin(2)*srate); %get plot window
%                     
%                     zmask = Zmaskixpc2plot;
%                     zmask(abs(Zmaskixpc2plot)<zval) = 0;
%                     
%                     [~,h] = contour(intfig,plottimeWin,frex(1:end-1),logical(zmask),1);
%                     h.LineColor = [.85 .85 .85];
%                     %                     alphable = findobj(h, '-property', 'FaceAlpha');
%                     %                     set(alphable, 'FaceAlpha', .5);
%                     hold on;
%                     MCthresh = MCminmax(ceil(length(MCminmax)*(1-pval)));
%                     MCthreshmap = ixpc2plot;
%                     MCthreshmap(abs(ixpc2plot)<MCthresh) = 0;
%                     [~,mc] = contour(intfig,plottimeWin,frex(1:end-1),logical(MCthreshmap),1);
%                     mc.LineColor = [.6 .6 .6];
%                     caxis([-1 1]); % i need to normalize the coloring
%                 else
%                     caxis([-1 1]); % i need to normalize the coloring
%                 end
%                 % figure
%                 %                 imagesc(ixpc2plot);
%                 %                 patch([-.8, .8, .8, -.8], [-10 -10 80 80], iNTcolor, 'edgecolor','none')
%                 hold on;
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
%                 colormap(usecolormap)
                
                
                if mod(int, sfcols) ~= 1;
                    set(gca, 'yscale','log','ytick',[])
                    %                                         set(gca, 'yscale','log','ytick',ceil(logspace(log10(min(frex)),log10(max(frex)),6)))
                    set(gca, 'FontSize',8, 'FontName', 'Arial');%
                else
                    set(gca, 'yscale','log','ytick',ceil(logspace(log10(min(frex)),log10(max(frex)),6)))
                    %                     set(gca, 'YTick',[round(linspace(min(frex),max(frex),6))])
                    set(gca, 'FontSize',8, 'FontName', 'Arial');%
                end
                if int <= (sfrows-1)*sfcols;
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
                                    
%             colormap(usecolormap)
            clrbar = colorbar('location','eastoutside', 'FontSize',6,'FontName', 'Arial');%, 'FontWeight','bold');
            caxis(coloraxis)
            posx1=get(gca,'position');
            posx=get(clrbar,'Position');
            posx(1)= 1-fig.MarginRight+.01;
            posx(2)= fig.MarginBottom;
            posx(3)= .01;
            posx(4)= 1 - fig.MarginTop - fig.MarginBottom;
            set(clrbar,'Position',posx)
            set(gca,'position',posx1)
            %             sprtit = sprintf('%s D%s E%s T%d', ianimalinfo{1}, strjoin(arrayfun(@(x) num2str(x),days','UniformOutput',false),'-'), strjoin(arrayfun(@(x) num2str(x),idayEpTet(:,2)','UniformOutput',false),'-'), idayEpTet(1,3));
            %             if iDT == 6
            %                 sprtit = sprintf('%s %s %s %s D%s %d-%dHz_rCyc%d-%d', plotoutputtype,'outB-inB', strjoin(epochEnvironment,'-'), animalID, strrep(num2str(days), '  ', '-'), min_freq,max_freq, range_cycles);
            %             elseif iDT == 7
            %                 sprtit = sprintf('%s %s %s %s D%s %d-%dHz_rCyc%d-%d', plotoutputtype,'outBCorr-outBMist', strjoin(epochEnvironment,'-'), animalID, strrep(num2str(days), '  ', '-'), min_freq,max_freq, range_cycles);
            %             else
            sprtit = sprintf('%s %s %s %s D%s %d-%dHz_rCyc%d-%d', plotoutputtype,ixpc.datatypes{ian}{iDT}, strjoin(epochEnvironment,'-'), animalID, strrep(num2str(days,'%-2d'), ' ', '-'), min_freq,max_freq, range_cycles);
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
            clrbartit = text(posx(1)+posx(3)/2, posx(2)-fig.MarginBottom/2, plotoutputtype, 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial','horizontalAlignment', 'center');
            
            %% ---- pause, save figs ----
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
    for ian = 1:length(ixpc.animals)
        %% ---- loadtetinfostruct ----
        animalinfo = animaldef(lower(animals{ian}));
        animalID = animalinfo{1,3}; %use anim prefix for name
        FFanimdir =  sprintf('%s',animalinfo{1,2});
        load([FFanimdir, animalID, 'tetinfo']);
        tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
        for iday = days
            ntrodePairIndices = ixpc.index{ian}{iday};
            tetpairs = unique(ntrodePairIndices(:,[3 6]),'rows','stable');
            ntetPairs = size(tetpairs,1);
            %% ---- reorder the LFP traces by suparea|area|subarea tags (in that priority) ----
            ntrodePairIndicesLIN = [ntrodePairIndices(:,[1:3]); ntrodePairIndices(:,[4:6])];
            wp.nNTrodes = size(unique(ntrodePairIndicesLIN, 'stable'),1);
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
                contourf(timeWin,frex,squeeze(ixpc.phOut{ian}{iday}(numsumSortInds(intrPair),:,:))',numfrex,'linecolor','none');
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

