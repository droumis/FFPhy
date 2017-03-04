warning('off', 'MATLAB:Figure:RecursionOnClose') %suppress this specific warning
warning('off', 'MATLAB:MKDIR:DirectoryExists') %suppress this specific warning



%DR load and plot ripple triggered LFP across regions

%>>>>>>>>>EXCLUDE RIPPLES OUTSIDE OF EPOCH TIME RANGE -- TO DO
%make window a bit smaller
%make population images
%look through all of them
%add other epochs
%add deep vs superficial tag and colorize categories

%%

%% plots selected tetrodes' EEG side-by-side around reward times, looking at ripples

% function out = DR_rippleTrig_LFP(directory,animal,day,epoch,regions,rippleregion, rippletype,windowsize)
clear all; close all;
animal = 'D13';
day = 1;
epoch = 2;
regions = [{'ca1'},{'mec'},{'por'}]; %regions to markup ripple times
% regions = [{'ca1'},{'mec'},{'por'},{'v2l'}];
rippleregion = 'ca1'; %region of ripples you want to align to
rippletype = 'kons'; %cons or kons or '' for old method
windowsize = 0.4; %seconds
LFPtypes = [{'eeg'} {'ripple'}]; %[{'eeg'} {'ripple'}] only option right now
plotscatters = 0;
plotLFPtraces = 1; % display figures
pausefigs = 1; % pause figs to manually inspect before or instead of saving
savefigs = 0; % saves pngs of individual LFPtype traces.. then createmergedplots will create new combined images; This will suppress visualization of figures to prevent fig tweaking
createmergedplots = 0; %system call to ImageMagick Convert append tool... usage may differ across linux distros
nstdres = 1; %increment/RESolution of nstd range.. default 1
minstdthresh = 5;

% ------------- param-def --------------------
    % get date string for saving plots
animalinfo = animaldef(lower(animal)); %get the animal info/data paths
directory =  animalinfo{1,2}; 
[~, dirpath] = fileattrib(directory);
[dirstrind] = regexp(dirpath.Name,'/');
figdirectory = strcat(dirpath.Name(1:(dirstrind(end))), 'figs'); %create path name to animal/figs from animal/filterframework path in animaldef
formatOut = 'yyyymmdd-HHMMSS';
dateprintstr = datestr(now,formatOut); %grab current date-time to add to filename string
    % plot aesthetics
regionclr = ([170,11,10; 80,215,167; 0,40,185; 255,170,55])./255;
patchclr = (1.-regionclr)/1.3 + regionclr;   % scaled lightening of regionclr for patch ([.9 .7 .7 ; 0.2 .6 .9; 0 0 1; 0 1 0]);
    %other stuff
srcRegionind = find(strcmp(rippleregion, regions),1); % index of the rippleregion in the regions cellstring
RippleEEGind = find(strcmp('ripple', LFPtypes),1); % index of the rippleregion in the regions cellstring
EEGind = find(strcmp('eeg', LFPtypes),1); % index of the rippleregion in the regions cellstring
tetinfo = loaddatastruct(directory,animal,'tetinfo',day);


%% STEP 1: get the consensus ripple times from the designated regions
if ~isempty(rippletype);
    ripplefile=sprintf('%s%s%sripples%s%02d.mat',directory,animal,rippleregion,rippletype, day); %cons or kons or dons
else %not implemented yet
    disp('old style rips not yet implemented beyond this point')
    pause
end
load(ripplefile);

clear ripout
for iArea = 1:length(regions)
    ripout{iArea} = DR_kk_getconstimes(directory,animal,[day epoch],sprintf('%srippleskons',regions{iArea}),1,'consensus_numtets', 1, 'maxvelocity',4, 'minstdthresh', minstdthresh); %%kk usually did concensus numtets 3 but i only had 2 tets for some regions
    [eventtimes, eventindices] = dr_vec2list(ripout{iArea}{day}{epoch}.cons,ripout{iArea}{day}{epoch}.time); %idk why we do this?? .. start end for every ripple is already extracted into the eventkons struct and it matches this
    %update:: i think this is done bc the start and end in the eventkons struct is for rips derived with the params used to run the extraction.. so if you extracted all rips > 2std but now want just >3std.. you have to recreate the start/endtimes
    ripout{iArea}{day}{epoch}.rippleTimes = eventtimes; %list of start and end times for each ripple for each area
    ripout{iArea}{day}{epoch}.rippleIndices = eventindices; %indices into the eeg traces for start and end of each ripple for each area
    ripout{iArea}{day}{epoch}.area = regions{iArea}; 
end

%% MOVE THIS TO NEW SCRIPT
%% return indices of ripples for each nstd range
% for i = 1:length(regions)
%     for p = 1:length(ripout{srcRegionind}{day}{epoch}.rippleIndices(:,1)) %for each ripple in src region
%         powertraces{i}{day}{epoch}{p} = ripout{i}{day}{epoch}.powertrace(ripout{srcRegionind}{day}{epoch}.rippleIndices(p,1):ripout{srcRegionind}{day}{epoch}.rippleIndices(p,2))';   %take the powertrace of each region within the range of each src rip and calculate the nstds of each (peak - baseline)/std
%         ripnstd{i}{day}{epoch}(p,1) = ((max(powertraces{i}{day}{epoch}{p}) - ripout{i}{day}{epoch}.baseline)/ripout{i}{day}{epoch}.std); %get the peak nstds above baseline with curr ripple window for src region
%     end
% end

% %% get nstds correlation for src vs all regions... ONLY needed for plotscatters
% for i = 1:length(regions)
%     ripnstdcorr(i,1) = corr(ripnstd{srcRegionind}{day}{epoch},ripnstd{i}{day}{epoch});
% end
%% plot scatters of each nstd range
% if 0
%     for i = 1:length(regions)
%         figure(i);
%         s = floor(min(nstd{srcRegionind}{day}{epoch}));
%         while s <= ceil(max(nstd{srcRegionind}{day}{epoch}));
%             %         validstdripinds = ones(nstd{srcRegionind}{day}{epoch} > s & nstd{srcRegionind}{day}{epoch} < (s + nstdres))*s; %  return indices of src ripples with a nstd value within the defined range
%             %         validstdripinds{s}.nstdrange =  [s (s + nstdres)];
%             ripsincurrstd = (ripnstd{srcRegionind}{day}{epoch} > s & nstd{srcRegionind}{day}{epoch} < (s+ nstdres));
%             %             scatter(nstd{srcRegionind}{day}{epoch}(ripsincurrstd), nstd{i}{day}{epoch}(ripsincurrstd), 'MarkerEdgeColor','none','MarkerFaceColor',rand(1,3),'LineWidth',1.5); hold on; %ca1 vs mec, ca1 vs por nstd peak
%             scatter(nstd{2}{day}{epoch}(ripsincurrstd), nstd{3}{day}{epoch}(ripsincurrstd), 'MarkerEdgeColor','none','MarkerFaceColor',rand(1,3),'LineWidth',1.5); hold on; %mec vs por nstd peak
%             s = s + nstdres;
%         end
%         if pausefigs; pause; end
%         if savefigs; disp('need to write this saving code'); end
%     end
% end

%% STEP 2: For each EEG + Ripple Type... select nTrodes and Lookup Rip time windows

for iLFPtype = 1:length(LFPtypes); % For each LFP type (wideband EEG, ripple band, etc), load all of the regions LFP files into eegstruct
    for iArea = 1:length(regions) % Parse 'regions' into tetfilters
        nTrodeFilter{iArea} = sprintf('(isequal($area,''%s''))', regions{iArea});
        tmpNTrode = evaluatefilter(tetinfo,nTrodeFilter{iArea}); % Get tets corresponding to current area
        selected_nTrodes{iArea} = unique(tmpNTrode((tmpNTrode(:,1) == day),3))';
        if ~selected_nTrodes{iArea};
            disp(sprintf('no valid ntrodes for %s detected... continue?', regions{iArea})); 
            pause;
        end
        eegstruct{iLFPtype}{iArea} = loadeegstruct(directory,animal,LFPtypes{iLFPtype},day,epoch,selected_nTrodes{iArea}); %load each region's lfp data {region}{day}{epoch}{ntrode}.FIELDS
    end
    tet4epochinfo = find(~cellfun(@isempty,eegstruct{iLFPtype}{1}{day}{epoch}),1); %find a non-empty cell from the first region
    num_samp=length(eegstruct{iLFPtype}{1}{day}{epoch}{tet4epochinfo}.data); % get number of samples
    samprate=eegstruct{iLFPtype}{1}{day}{epoch}{tet4epochinfo}.samprate; %get sampling rate
    epochLFP_starttime=double(eegstruct{iLFPtype}{1}{day}{epoch}{tet4epochinfo}.starttime); %get epoch start time
    
    %update... nan filled any gaps in the lfp data so i can now actually use this method of getting timestamp
    LFPtimes=(epochLFP_starttime:1/samprate:(epochLFP_starttime+(1/samprate)*num_samp))'; % prob should be using the adjusted timestamps instead, no?

    %% Lookup all of the ripple time windows
    % DR FIX.. USE THE LOOKUP MEX INSTEAD OF KNNSEARCH... is this faster?
    ripsStartTimesind = cellfun(@(x) lookup(x{day}{epoch}.rippleTimes(:,1),LFPtimes), ripout, 'UniformOutput', false); %map all ripple start times from all the regions to nearest LFP timestamp indicies
    ripsEndTimesind = cellfun(@(x) lookup(x{day}{epoch}.rippleTimes(:,2),LFPtimes), ripout, 'UniformOutput', false); %map all ripple start times from all the regions to nearest LFP timestamp indicies
%     ripsStartTimesind = cellfun(@(x) knnsearch(LFPtimes, x{day}{epoch}.rippleTimes(:,1)), ripout, 'UniformOutput', false); %map all ripple start times from all the regions to nearest LFP timestamp indicies
%     ripsEndTimesind = cellfun(@(x) knnsearch(LFPtimes, x{day}{epoch}.rippleTimes(:,2)), ripout, 'UniformOutput', false); %map all ripple start times from all the regions to nearest LFP timestamp indicies
    tmpripStartEndtimesInd = cellfun(@(x,y) [x y], ripsStartTimesind, ripsEndTimesind, 'UniformOutput', false);
    ripsStartEndTimesind = cellfun(@(x) x(x(:,1) > (windowsize*samprate) & x(:,2) < (length(LFPtimes)-(windowsize*samprate)),:), tmpripStartEndtimesInd, 'UniformOutput', false); %keep all the ripples that occur after windowsize and before (end of epoch minus windowsize)
    ripsamplengths = cellfun(@(x) x(:,2)-x(:,1),ripsStartEndTimesind,'UniformOutput', false); %get # samples length of each ripple
    ripsStartTimes = cellfun(@(x) LFPtimes(x(:,1)),ripsStartEndTimesind,'UniformOutput', false); %get rip start timestamp values from LFPtimes indices
    ripsEndTimes = cellfun(@(x) LFPtimes(x(:,2)),ripsStartEndTimesind,'UniformOutput', false); %get rip start timestamp values from LFPtimes indices
    WindowStartEndTimes = [ripsStartTimes{srcRegionind}(:)-windowsize ripsStartTimes{srcRegionind}(:)+windowsize]; %get window range for all the rips from the src regions

    %% STEP 3: Gather the ripple window data from all the regions 
    for currrip=1:length(ripsStartTimes{srcRegionind}) %for each ripple from source region
        clear YripLFPdata %precautionary 
        for iArea = 1:length(regions);
            nTrode_count = 0; 
            for iNTrode = selected_nTrodes{iArea}
                nTrode_count = nTrode_count + 1;
                if strcmp(LFPtypes{iLFPtype},'eeg') || strcmp(LFPtypes{iLFPtype},'ripple')
                    % Gather LFPtype data for each area within the current rip window
                    YripLFPdata{iArea}(:,nTrode_count) = double(eegstruct{iLFPtype}{iArea}{day}{epoch}{iNTrode}.data(ripsStartEndTimesind{srcRegionind}(currrip,1)-(windowsize*samprate):ripsStartEndTimesind{srcRegionind}(currrip,1)+(windowsize*samprate)));
                else
                    disp(sprintf('incompatible LFPtype: %s .. currently only tested with eeg and ripple', LFPtypes{iLFPtype}))
                    return
                end
                lfptraceregion{iLFPtype}{iArea} = ones(length(YripLFPdata{iArea}(1,:)),1)*iArea;
            end
            
            %MOVE THIS TO NEW SCRIPT
%             %calculate all the nstds for each region for all tetrodes in the time windows of this src-region-ripple..
%             %instead of using the std and baseline from the powertrace used for extracting ripples.. recalculate those here within each pre-rip window so it can be used for full band or other band LFP
%             baseline{iLFPtype}{iArea}{day}{epoch} = mean(YripLFPdata{iArea}(1:(windowsize*samprate),:),1);
%             currripstartind = (windowsize*samprate)+1;
%             currripendind = ((windowsize*samprate)+1+ripsamplengths{srcRegionind}(currrip));
%             if currripendind > length(YripLFPdata{iArea}(:,1)) %if the computed rip end is longer than the rip window, truncate the end to the window
%                 currripendind = (windowsize*samprate)*2+1;
%             end
%             try
%                 allstd{iLFPtype}{iArea}{day}{epoch}(currrip,:) = std(YripLFPdata{iArea}(currripstartind:currripendind,:),1);
%                 nstd{iLFPtype}{iArea}{day}{epoch}(currrip,:) = (max(YripLFPdata{iArea}(currripstartind:currripendind,:),[],1) - baseline{iLFPtype}{iArea}{day}{epoch})./allstd{iLFPtype}{iArea}{day}{epoch}(currrip,:);
%                 %computer FFT and spectrogram for each snippet
%                 
%             catch
%                 return
%                 allstd{LFP}{i}{day}{epoch}(currrip,:) = [];
%                 nstd{LFP}{i}{day}{epoch}(currrip,:) = [];
%                 
%             end
        end
        YripLFPdataMAT{iLFPtype}{currrip} = cell2mat(YripLFPdata); %stack all the traces from all the regions next to each other
        Xwindowtimes{iLFPtype}{currrip} = WindowStartEndTimes(currrip,1):1/samprate:WindowStartEndTimes(currrip,2);
        lfptraceLUTregion = cell2mat(lfptraceregion{iLFPtype}'); % [1 x length(nTrodes)] vec mapping each lfp trace to a region in 'regions' by index for visualization
    end
    
%% MOVE THIS TO NEW SCRIPT
     %% plot scatters of each nstd range  NEED TO CLEAN this section up and save figs and format
%     if plotscatters && iLFPtype == length(LFPtypes); %if it's on the last lfp type loop
%         close all
%         fig1 = figure(1);
%         fig2 = figure(2);
%         fig3 = figure(3);
%         fig4 = figure(4);
%         s = floor(min(ripnstd{srcRegionind}{day}{epoch})); %find the min nstd of the average nstd for each rip for the src rip region
%         while s <= ceil(max(ripnstd{srcRegionind}{day}{epoch})); %iterate until reaching the max average rip nstd
%             ripsincurrstd = (ripnstd{srcRegionind}{day}{epoch} >= s & ripnstd{srcRegionind}{day}{epoch} < (s+ nstdres)); % all rips from the kons powertrace within curr nstd range
%             %             scatter(nstd{srcRegionind}{day}{epoch}(ripsincurrstd), nstd{i}{day}{epoch}(ripsincurrstd), 'MarkerEdgeColor','none','MarkerFaceColor',rand(1,3),'LineWidth',1.5); hold on; %ca1 vs mec, ca1 vs por nstd peak
%             figure(fig1)
%             scatter(max(nstd{EEGind}{2}{day}{epoch}(ripsincurrstd),[],2), max(nstd{EEGind}{3}{day}{epoch}(ripsincurrstd),[],2), 'MarkerEdgeColor','none','MarkerFaceColor',rand(1,3),'LineWidth',1.5); hold on; %mec vs por nstd peak
%             figure(fig2)
%             scatter(max(nstd{RippleEEGind}{2}{day}{epoch}(ripsincurrstd),[],2), max(nstd{RippleEEGind}{3}{day}{epoch}(ripsincurrstd),[],2), 'MarkerEdgeColor','none','MarkerFaceColor',rand(1,3),'LineWidth',1.5); hold on; %mec vs por nstd peak
%             figure(fig3)
%             scatter(ripnstd{srcRegionind}{day}{epoch}(ripsincurrstd), max(nstd{EEGind}{2}{day}{epoch}(ripsincurrstd),[],2), 'MarkerEdgeColor','none','MarkerFaceColor',rand(1,3),'LineWidth',1.5); hold on; %mec vs por nstd peak
%             figure(fig4)
%             scatter(ripnstd{srcRegionind}{day}{epoch}(ripsincurrstd), max(nstd{EEGind}{3}{day}{epoch}(ripsincurrstd),[],2), 'MarkerEdgeColor','none','MarkerFaceColor',rand(1,3),'LineWidth',1.5); hold on; %mec vs por nstd peakak
%             s = s + nstdres;
%         end
%         if pausefigs; return; end
%         if savefigs; disp('need to write this saving code'); end
%     end

%% plot LFP traces for all areas for each ripple time window
    if plotLFPtraces
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
                    lfpoffset = 0;
                else %compute offset from the absMax+absMin to place the current trace under the last trace
                    lfpoffset(p) = lfpoffset(p-1) + abs(min(YripLFPdataMAT{iLFPtype}{currrip}(:,p-1))) + abs(max(YripLFPdataMAT{iLFPtype}{currrip}(:,p))); %use e.g. abs(max(YripLFPdataMAT{LFP}{currrip}(:,p)))/3 to overlay traces more
                    plot(Xwindowtimes{iLFPtype}{currrip},YripLFPdataMAT{iLFPtype}{currrip}(:,p) - lfpoffset(p), 'Color',regionclr(lfptraceLUTregion(p),:), 'LineWidth',1); hold on;
                end
            end
            %% Plot Ripple times as patchs and other accessory ploting stuff
            for iArea = 1:length(regions); %for each region
                ripsinwinInds = (ripsStartTimes{iArea}>WindowStartEndTimes(currrip,1) & ripsStartTimes{iArea}<WindowStartEndTimes(currrip,2)); %if the start time is within the current window
                ripsinwin{currrip}{iArea} = [ripsStartTimes{iArea}(ripsinwinInds) ripsEndTimes{iArea}(ripsinwinInds)];
                if ~isempty(ripsinwin{currrip}{iArea}) %if there are any ripples from this region in this window
                    for m = 1:length(ripsinwin{currrip}{iArea}(:,1)) %for each ripple within the current window
                        Ylfpranges4region{iArea} = [-lfpoffset(find(lfptraceLUTregion == iArea,1,'last'))-(abs(min(YripLFPdataMAT{iLFPtype}{currrip}(:,find(lfptraceLUTregion == iArea,1,'last'))))) -lfpoffset(find(lfptraceLUTregion == iArea,1,'first'))+max(YripLFPdataMAT{iLFPtype}{currrip}(:,find(lfptraceLUTregion == iArea,1,'first')))];
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
            Yextent = [-lfpoffset(end)-(abs(min(YripLFPdataMAT{iLFPtype}{currrip}(:,end)))) max(max(YripLFPdataMAT{iLFPtype}{currrip}))];
            line([centerripStartTime centerripStartTime], Yextent ,'Color',regionclr(srcRegionind,:), 'LineStyle', '--','LineWidth',1.5) %line for the center trigger-ripple
            set(gca,'children',flipud(get(gca,'children'))) %send the patch behind the LFP traces
            ylim([Yextent(1) Yextent(2)])
            xlim([WindowStartEndTimes(currrip,1) WindowStartEndTimes(currrip,2)])
            xl = xlim;
            for k=1:length(regions)
                text(xl(1)-diff(WindowStartEndTimes(currrip,:))/10, -lfpoffset(find(lfptraceLUTregion == k,1,'first')), regions{k}, 'Color', regionclr(k,:),'FontSize',13,'FontWeight','bold')
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
        end
    end
end

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