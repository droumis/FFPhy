
% to do:
% ITPC - plot open field
% ISPC - preserve the indices

%%% First create the processed data structure with dfs_riptriglfp.m

close all
saveFilterOutput = 0;
loadFilterOutput = 0;
calculateIXPC = 0;
calculateITPC = 0;
calculateISPC = 0;
% useCombinedEpochs = calculateIXPC;
saveResultsOutput = 0;
loadResultsOutput = 0;
% plotNTrodesAcrossDays = 0;
plotITPC = 1;
plotISPC = 0;
savefigs = 0;
pausefigs = 1;


%% ---------------- plotting params --------------------------
colorSet = 'DR1';
clims = [0 .7];
position = [.1 .1 .9 .8]; %[.1 .1 .9 .8]
SpacingHorizontal = 0.01;
SpacingVertical = 0.02;
% Spacing = 0.00;
Padding = 0.0;
MarginLeft = 0.04;
MarginRight = 0.04;
MarginTop = 0.09;
MarginBottom =  0.08;
usecolormap = 'jet';
win = [.5 .5]; %in seconds
% indwin = win*1500;
%% ---------------- Data Filters --------------------------
% wavelet parameters
num_frex = 60;
min_freq =  4;
max_freq = 60;
% set range for variable number of wavelet cycles
range_cycles = [ 4 10 ];
win = [-.5 .5];
srate = 1500;
timeWin = win(1):1/srate:win(2);
calcfunction = 'ITPC'; %ISPC or ITPC
animals = {'JZ1'};
% animals = {'JZ1', 'D13'};
days = [1];
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
outputDirectory = '/typhoon/droumis/analysis';
%% ---------------- Paths and Title strings ---------------------------------------------------
currfigdirectory = sprintf('%s/figures/%s/', outputDirectory, calcfunction);
filenamesave = sprintf('%s%s_%s_%s_%s', eventSourceArea, eventtype, epochEnvironment, cell2mat(animals), cell2mat(LFPtypes));
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
filenameTitle = strrep(filename,'_', '\_');
resultfilename = sprintf('%s_%s%s_%s_%s_D%s_%d-%dHz', calcfunction, eventSourceArea, eventtype, epochEnvironment, cell2mat(animals), strrep(num2str(days), '  ', '-'), min_freq,max_freq);
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s/filter_output/%s/%s',outputDirectory,filtfunction, filename))
end
%% Calculate inter-trial phase clustering based on mike cohen's book

% other wavelet parameters
frex = logspace(log10(min_freq),log10(max_freq),num_frex);
wavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex);
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
        for iday = days
%             if useCombinedEpochs
                eps = find(~cellfun(@isempty,{F(ianimal).output{iday}.index})); %get nonempty eps
                iEpTetData = [];
                for iep = eps;
                    iEpTetData = cat(2,iEpTetData,F(ianimal).output{iday}(iep).data{1,1}); %cat the event snips across eps
                end
%             else
%                 error('stahpit')
                %     iEpTetData = F.output{1, 1}(ep2use).data{1, 1}
%             end
            %             iEpTetDataCat = cat(3, iEpTetData{:}); % cat stack the rip snips into the 3rd dim
            
            iEpTetDataCat = cat(3, iEpTetData{:}); % cat stack the rip snips into the 3rd dim
            % FFT parameters
            nWave = length(timeWin);
            nNTrodes = length(iEpTetDataCat(:,1,1));
            nsamps = length(iEpTetDataCat(1,:,1));
            nevents = length(iEpTetDataCat(1,1,:));
            nData = nsamps*nevents;
            nConv = nWave+nData-1;
            
            % FFT of data (doesn't change on frequency iteration)
            % dataX = fft(reshape(iEpTetDataCat(introde,:,:), 1,nData),nConv,2);
            dataY = fft(reshape(iEpTetDataCat,nNTrodes,nData),nConv,2); %concat peri-rip snips
            indices = [];
            if calculateISPC
                % initialize output time-frequency data
                ixpc.output{ianimal}{iday} = zeros(nchoosek(size(iEpTetDataCat,1),2), nsamps, num_frex);
            elseif calculateITPC
                % initialize output time-frequency data
                ixpc.output{ianimal}{iday} = zeros(nNTrodes,nsamps, num_frex);
                indices = [indices; F(ianimal).output{iday}(eps(1)).index];
            end
            
            % loop over frequencies
            for fi=1:num_frex
                % create wavelet and get its FFT
                s = wavecycles(fi)/(2*pi*frex(fi));
                wavelet  = exp(2*1i*pi*frex(fi).*timeWin) .* exp(-timeWin.^2./(2*s^2));
                waveletX = fft(wavelet,nConv);
                
                % run convolution
                as = ifft(bsxfun(@times,dataY,waveletX),nConv,2);
                as = as(:,half_wave_size+1:end-half_wave_size);
                as = reshape(as,nNTrodes,nsamps,nevents);
                phasedata = angle(as); %get phase component of the analytic signal
                if calculateISPC % time series diff of all possible ntrode pairs. 
                    %the idea here is to treat rows as pairs instead of individual ntrodes, as in ITPC
                    ispc = zeros(nchoosek(size(phasedata,1),2), nsamps, nevents);
                    m = logical(tril(ones(size(phasedata(:,:,1),1)),-1)); %get indices of non-duplicates (below comb triangle)
                    [brow, bcol] = find(m); %get linear index of ntrode combination indices
                    indices = [F(ianimal).output{iday}(eps(1)).index(bcol,:) F(ianimal).output{iday}(eps(1)).index(brow,:)]; %convert to ntrodeID
                    for w = 1:size(phasedata,3); %loop through each event plane (ntrode x samples X event#)
                        iprm = permute(phasedata(:,:,w),[3 2 1]); %transpose
                        B = bsxfun(@minus,phasedata(:,:,w),iprm); %subtract across each NTrode-choose-two comb of phase time series
                        B = reshape(B,[],size(phasedata(:,:,w),2)); %reshape back into 2D
                        B = B(m(:),:); %get rid of duplicates
                        ispc(:,:,w) = B; %save result
                    end
                    clear phasedata
                    phasedata = ispc;
                end
                % compute IXPC
                ixpc.output{ianimal}{iday}(:,:,fi) = abs(mean(exp(1i*phasedata),3));
                ixpc.waveletX{ianimal}{iday} = waveletX;
                ixpc.index{ianimal}{iday} = indices;
                disp(sprintf('day %d freq %d of %d',iday,fi,num_frex));
            end
            disp(sprintf('==========day %d calculated=============',iday));
        end
    end
    toc
end
%% ---------------- Save RESULTS Output ---------------------------------------------------
if saveResultsOutput == 1;
    if ~isdir(sprintf('%s/results_output/%s/', outputDirectory, calcfunction));
        mkdir(sprintf('%s/results_output/%s/', outputDirectory, calcfunction));
    end
    save(sprintf('%s/results_output/%s/%s',outputDirectory, calcfunction, resultfilename), 'ixpc','-v7.3');
    disp(sprintf('%s saved', resultfilename))
end
%% ---------------- Load Results Output ---------------------------------------------------
if loadResultsOutput == 1;
    load(sprintf('%s/results_output/%s/%s',outputDirectory, calcfunction, resultfilename));
end
%% ---------------- plot ITPC---------------------------------------------------------------------------------------------
%% ---------------- plot ITPC---------------------------------------------------------------------------------------------
%% ---------------- plot ITPC---------------------------------------------------------------------------------------------
if plotITPC
    clear F %save space on memory
    for ianimal = 1:length(ixpc.animals)
        %% ---- loadtetinfostruct ----
        animalinfo = animaldef(lower(animals{ianimal}));
        animalID = animalinfo{1,3}; %use anim prefix for name
        FFanimdir =  sprintf('%s',animalinfo{1,2});
        load([FFanimdir, animalID, 'tetinfo']);
        tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
        for iday = days
            ntrodesIndices = ixpc.index{ianimal}{iday};
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
            if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
                ifig = figure('Visible','off','units','normalized','position',[.1 .1 .9 .8]);
            else
                ifig = figure('units','normalized','position',position);
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
                contourf(intfig,timeWin,frex,squeeze(ixpc.output{ianimal}{iday}(numsumSortInds(introde),:,:))',num_frex,'linecolor','none');%, )
%                 patch([-.8, .8, .8, -.8], [-10 -10 80 80], iNTcolor, 'edgecolor','none')
                hold on;
                set(gca,'clim',clims,'ydir','normal','xlim',[-.5 .5])
                set(gca, 'YScale', 'log')

                if mod(introde, sfcols) ~= 1;
                    set(gca, 'YTick', []);
                else
                    set(gca, 'YTick',[round(linspace(min(frex),max(frex),6))],'FontSize',8, 'FontName', 'Arial');%
                end
                if introde <= (sfrows-1)*sfcols;
                    set(gca, 'XTick', []);
                else
                    set(gca, 'XTick',[win(1):win(2)/2:win(2)],'FontSize',8, 'FontName', 'Arial');%
                end
                %% ---- Source ripple line, subplot title ----
                Xline = [0 0];
                Yline = [min_freq max_freq];
                line(Xline, Yline, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                ititle = title(sprintf('nT%d %s %s',introdeID, iareatag, num2str(isubareatag)));
                set(ititle,'FontSize',10,'Color', iNTcolor, 'FontName', 'Arial','FontWeight','bold'); %
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
            sprtit = sprintf('%s %s %s D%d %d-%dHz', calcfunction, epochEnvironment, animalID, iday, min_freq,max_freq);
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
            clrbartit = text(posx(1)+posx(3)/2, posx(2)-MarginBottom/2, 'ITPC', 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial','horizontalAlignment', 'center');
            %% ---- pause, save figs ----
%             return
            if pausefigs
                pause
            end
            if savefigs
                if ~isdir(currfigdirectory);
                    mkdir(currfigdirectory);
                end
%                 if ~isdir([currfigdirectory filenamesave]);
%                     mkdir([currfigdirectory filenamesave]);
%                 end
                sprtitsave = strrep(sprtit,' ', '_'); %put '_' character in place of whitespace for filename
                set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
                currfigfile = sprintf('%s/%s',currfigdirectory, sprtitsave);
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


