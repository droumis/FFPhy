
% Similar to DFSsj_HPexpt_getripalignspiking and DFSsj_HPexpt_getripalignspiking_sleep
% Except I am not combining across epochs. This is to save a structure with the ripple aligned responses,
% so that it can be loaded through the filter framework. Can also do statistics separately for epochs.
% For all cells.


% Ripple modulation of cells. Time filter version of sj_HPexpt_ripalign_singlecell_getrip4.
% Will call DFAsj_getripalign.m
% Also see DFSsj_plotthetamod.m and DFSsj_HPexpt_xcorrmeasures2. Will gather data like these

clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - Individual cells
%
savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
savefile = [savedir 'HPc_ripplemod_all_epoch']; % All cells

% savefile = [savedir 'HP_ripplemod_all_epoch']; % All cells
% %savefile = [savedir 'HP_ripplemod_PFC_epoch']; area = 'PFC'; clr = 'b'; % PFC
% %savefile = [savedir 'HP_ripplemod_CA1_epoch']; area = 'CA1';  clr = 'r';% CA1
% %savefile = [savedir 'HP_ripplemod_PFC_speed_epoch']; area = 'PFC'; clr = 'b'; % PFC - low speed criterion
%
% savefig1=0;
%
%
% % Plot options
% plotanimidx =  []; % To pick animals for plotting
% plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days
%
%
% % If runscript, run Datafilter and save data
% if runscript == 1
%
%     %Animal selection
%     %-----------------------------------------------------
%     %animals = {'HPa','HPb'};
%     animals = {'HPc'};
%
%     %Filter creation
%     %-----------------------------------------------------
%
%     % Epoch filter
%     % -------------
%     dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
%     % Either Only do 1st w-track. 2 or 1 epochs per day
%     % Or do Wtr1 and Wtr2, 2 epochs per day
%     allepochfilter = 'isequal($type, ''run'') || isequal($type, ''sleep'')'; % All epochs
%     %runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
%     %sleepepochfilter = 'isequal($type, ''sleep'')'; % Only pre and post sleep marked as sleep
%
%     % Cell filter
%     % -----------
%     cellfilter = 'strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') || strcmp($area, ''PFC'')';
%     %cellfilter = '~isempty($spikewidth)'; % All defined cells
%
%     % Time filter
%     % -----------
%     riptetfilter = '(isequal($descrip, ''riptet''))';
%
%     % The following are ripple time filter options. Instead of using time filters, within the function,
%     % call getripples using the riptetfilter option to get the riptimes. Since you want a vector of riptimes.
%     % You can also use inter-ripple-interval of 1 sec within the function, etc.
%
%     % a) Using getriptimes: generates 1 vector with 1's for ripple times for all tetrodes.
%     % Should be similar to getripples/ getripplees_direct
%
%     %     timefilterrun_rip = {{'DFTFsj_getvelpos', '(($absvel <= 5))'},...
%     %         {'DFTFsj_getriptimes','($nripples > 0)','tetfilter',riptetfilter,'minthresh',3}};
%
%     % b) getripples. This is similar to getripltimes, but gives start and end of each ripple time
%     % instead of a time filter. Also has an additional condition of at least 50 ms ripple.
%
%     % Iterator
%     % --------
%     iterator = 'singlecellanal';
%
%     % Filter creation
%     % ----------------
%     modf = createfilter('animal',animals,'days',dayfilter,'epochs',allepochfilter, 'cells',...
%         cellfilter, 'iterator', iterator);
%
%     disp('Done Filter Creation');
%
%     % Set analysis function
%     % ----------------------
%
%     modf = setfilterfunction(modf,'DFAsj_getripalignspiking',{'spikes', 'ripples', 'tetinfo', 'pos'}); %
%     %modf = setfilterfunction(modf,'DFAsj_getripalignspiking',{'spikes', 'ripples', 'tetinfo', 'pos'},'dospeed',1); %
%     % Going to call getripples_tetinfo within function
%
%     % Run analysis
%     % ------------
%     modf = runfilter(modf);
%     disp('Finished running filter script');
%     %--------------------- Finished Filter Function Run -------------------
%
%     if savedata == 1
%         clear figopt1 runscript plotdays plotanimidx savedata
%         save(savefile,'-v7.3');
%         %save(savefile);
%     end
%
% else
%
%     load(savefile);
%
% end % end runscript
%
% if ~exist('savedata')
%     return
% end


% -------------------------  Filter Format Done -------------------------



% %----------------------------------
% %Whether to gather data or to load previously gathered data
% %--------------------------------------------------------------------
gatherdata = 0; savegatherdata = 0;
gatherdatafile = [savedir 'HPc_ripplemod_all_epoch_gather'];

%gatherdatafile = [savedir 'HP_ripplemod_all_epoch_gather'];
%gatherdatafile = [savedir 'HP_ripplemod_PFC_epoch_gather']; % PFC cells to Hipp ripples
%gatherdatafile = [savedir 'HP_ripplemod_CA1_epoch_gather']; % CA1 cells to Hipp ripples
%gatherdatafile = [savedir 'HP_ripplemodspeed__epochPFC_gather']; % PFC cells to Hipp ripples - low speed criterion


if gatherdata
    
    % Parameters if any
    % -----------------
    
    % -------------------------------------------------------------
    
    cnt=0; % Count how many cells will be kept based on nspikes in output: >0
    allanimindex=[]; alldataraster=[]; alldatahist = []; all_Nspk=[];
    
    for an = 1:length(modf)
        for i=1:length(modf(an).output{1})
            % Check for empty output - If Cell defined in rpoch and Nspks in ripple response wndow > 0
            % if (modf(an).output{1}(i).Nspikes > 0)
            cnt=cnt+1;
            anim_index{an}(cnt,:) = modf(an).output{1}(i).index;
            % Only indexes
            animindex=[an modf(an).output{1}(i).index]; % Put animal index in front
            allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
            % Data
            alldataraster{cnt} = modf(an).output{1}(i).rip_spks_cell; % Only get raster and histogram response
            alldatahist{cnt} = modf(an).output{1}(i).rip_spkshist_cell;
            all_Nspk(cnt) = modf(an).output{1}(i).Nspikes;
            alldataraster_rdm{cnt} = modf(an).output{1}(i).rdm_spks_cell; % Only get raster and histogram response
            alldatahist_rdm{cnt} = modf(an).output{1}(i).rdm_spkshist_cell;
            % trialResps: Summed Nspks/trial in respective window
            alldatatrialResps{cnt} = modf(an).output{1}(i).trialResps;
            alldatatrialResps_bck{cnt} = modf(an).output{1}(i).trialResps_bck;
            alldatatrialResps_rdm{cnt} = modf(an).output{1}(i).trialResps_rdm;
            %end
            if cnt==1
                pret =  modf(an).output{1}(i).pret;
                postt = modf(an).output{1}(i).postt;
                binsize = modf(an).output{1}(i).binsize;
                rwin = modf(an).output{1}(i).rwin;
                bckwin = modf(an).output{1}(i).bckwin;
                bins_resp  = modf(an).output{1}(i).bins_resp;
                bins_bck = modf(an).output{1}(i).bins_bck;
                timeaxis = modf(an).output{1}(i).timeaxis;
            end
            % end
        end
    end
    
    
    
    %Save data separately for epochs.
    % ----------------------------------------------------------------------------
    
    allripplemod = struct;
    for i=1:length(alldatahist)
        animdaytetepcell=allanimindex(i,[1 2 3 4 5]);
        allripplemod_idx(i,:) = animdaytetepcell;
        allripplemod(i).index = animdaytetepcell;
        allripplemod(i).dayeptetcell = animdaytetepcell(2:end);
        try
            allripplemod(i).hist = alldatahist{i}*(1000/binsize); % Convert to firing rate in Hz
        catch
            keyboard;
        end
        allripplemod(i).raster = alldataraster{i};
        allripplemod(i).Nspk = all_Nspk(i);
        %allripplemod(i).hist_rdm = alldatahist_rdm{i}*(1000/binsize); % Convert to firing rate in Hz
        %allripplemod(i).raster_rdm = alldataraster_rdm{i};
        % Trial Resps
        allripplemod(i).trialResps = alldatatrialResps{i}';
        allripplemod(i).trialResps_rdm = alldatatrialResps_rdm{i}';
        allripplemod(i).trialResps_bck = alldatatrialResps_bck{i}';
        
    end
    
    % Calculations/ Stats. Stats between response and bck
    % Similar to ...getrip4
    
    for i=1:length(allripplemod)
        
        currNspk = allripplemod(i).Nspk;
        currhist = allripplemod(i).hist; %currraster = allripplemod(i).raster;
        %currhist_rdm = allripplemod(i).hist_rdm; %currraster_rdm = allripplemod(i).raster_rdm;
        
        % Get the bckresp, rresp and rdmresp again - using firing rate matrices
        rresp = currhist(:,bins_resp);
        bckresp = currhist(:,bins_bck);
        %rresp_rdm = currhist_rdm(:,bins_resp);
        
        % Bck
        avgbckresp_trial = mean(bckresp,2); avgrespbck_trial = avgbckresp_trial; % Mean in bck for each ripple
        %trialRespsh_bck = sum(bckresp,2);  % Nspikes/ trial in resp window
        avgbckhist = mean(bckresp); % Avg bck histogram
        mean_bckresp = mean(mean(bckresp)); % Single value
        distr_bckresp = bckresp(:); %All values taken by bins in background
        
        % Response
        avgresp_trial = mean(rresp,2); % Mean for each ripple
        %trialRespsh = sum(rresp,2); % Nspikes/ trial in resp window
        avgresphist = mean(rresp,1); % Avg resp histogram
        mean_rresp = mean(mean(rresp)); % Single value
        
        % RandomResponse
        %avgresprdm_trial = mean(rresp_rdm,2); % Mean for each random ripple
        %trialRespsh_rdm = sum(rresp_rdm,2); % Nspikes/ trial in resp window
        
        sig_shuf = 0; sig_ttest = 0;
        
        % 0) Simple t-test
        % -----------------
        [sig_ttest, p] = ttest2(avgbckresp_trial, avgresp_trial);
        
        
        % % 1) Significance test - USE SHUFFLE BETWEEN MEAN RESP AND MEAN BCK FOR EACH TRIAL
        % % ---------------------------------------------------------------------------
        % Get the actual mean difference between bck and resp
        Dm = mean(avgresp_trial) - mean(avgbckresp_trial); % DONT want absolute value. want to shuffle.
        
        % Shuffle mean bck and mean resp 1000 times
        comb = [avgresp_trial; avgbckresp_trial];
        ntr = size(comb,1);
        nshuffles = 1000;
        for shufidx = 1:nshuffles
            order = randperm(ntr);
            shuffle = comb(order,:);
            
            shufresp = shuffle(1:ntr/2,:); shufavgresp = mean(shufresp);
            shufbck = shuffle((ntr/2)+1:ntr,:); shufavgbckresp = mean(shufbck);
            %if shufidx==1, figure; hold on; end
            %plot(shufavgresp,'b'); plot(shufavgrdmresp,'g')
            Dshuf(shufidx) = shufavgresp - shufavgbckresp;
        end
        % Can Plot the distribution of shuffled values
        histD = histc(Dshuf,min(Dshuf):0.1:max(Dshuf));
        %figure; hold on; plot([min(Dshuf):0.1:max(Dshuf)],histD)
        
        % Get significance by comparing Dm to Dshuf
        % ------------------------------------------
        if Dm>=0
            pshuf = length(find(Dshuf>Dm))/nshuffles;
            if Dm > prctile(Dshuf,95)
                sig_shuf = 1;
            end
            type = 'exc'; peakresp = max(avgresphist); % Peak in response histogram
        else
            pshuf = length(find(Dshuf<Dm))/nshuffles;
            if Dm < prctile(Dshuf,5)
                sig_shuf = 1;
            end
            type = 'inh'; peakresp = min(avgresphist); % Trough in response histogram
        end
        % Get the p-value of shuffle. The modulation index for shuffle will be the prctile value
        % ----------------------------
        modln_shuf = 100 - (pshuf*100); %eg p=0.05 => prctile is 95
        % Get %tage change over baseline and peak/trough - Save with sign- +ve or -ve modln. Not abs
        % ----------------------------------------------
        modln = 100*(Dm)/mean(avgbckresp_trial); % Mean %tage change above/below baseline
        peakchange = peakresp - mean(avgbckhist); % mean(avgbckhist) is same as mean(avgbckresp_trial)
        modln_peak = 100*(peakchange)/mean(avgbckresp_trial); % Peak/trough %tage change above/below baseline
        
        % Save
        % -----
        allripplemod(i).Dm = Dm;
        allripplemod(i).pshuf = pshuf; allripplemod(i).p = pshuf;
        allripplemod(i).sig_shuf = sig_shuf; allripplemod(i).h = sig_shuf; % Sig or not
        allripplemod(i).sig_ttest = sig_ttest; % Sig or not
        allripplemod(i).modln_shuf = modln_shuf ; % %tile value of shuffle
        allripplemod(i).modln_peak = modln_peak; % %tage peak change over baseline
        allripplemod(i).modln = modln; % %tage mean change over baseline
        allripplemod(i).type = type; % exc or inh
        allripplemod(i).anim = allripplemod(i).index(1); allanim(i) = allripplemod(i).index(1);
        allripplemod(i).days = allripplemod(i).index(2); alldays(i) = allripplemod(i).index(2);
        
        % Mean resp from histogram in respective window
        %allripplemod(i).avghisttrialResps = avgresp_trial;
        %allripplemod(i).avghisttrialResps_rdm = avgresprdm_trial;
        %allripplemod(i).avghisttrialResps_bck = avgrespbck_trial;
        
        allsig_shuf(i) = sig_shuf;
        allsig_ttest(i) = sig_ttest;
    end
    
    
    % Save
    % -----
    if savegatherdata == 1
        clear modf
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
    
    if 0
        % ALSO, SAVE A DATA STRUCTRE, eg. "HParipplemod01". Will have to save each day separately
        % -----------------------------------------------------------------------------------------
        animvector = allripplemod_idx(:,1); dayvector = allripplemod_idx(:,2);
        animals = unique(allripplemod_idx(:,1));
        
        for an = 1:length(animals)
            curranim = animals(an);
            %prefix='HPc'; currsavedir = '/data25/sjadhav/HPExpt/HPc_direct/';
            switch curranim
                case 1
                    prefix='HPa'; currsavedir = '/data25/sjadhav/HPExpt/HPa_direct/';
                    
                case 2
                    prefix='HPb'; currsavedir = '/data25/sjadhav/HPExpt/HPb_direct/';
                    
                case 3
                    prefix='HPc'; currsavedir = '/data25/sjadhav/HPExpt/HPc_direct/';
            end
            
            %prefix = modf(1).animal{curranim};
            
            curranimidxs = find(animvector == curranim);
            curralldays = dayvector(curranimidxs);
            days = unique(curralldays);
            for d=1:length(days)
                ripplemod = [];
                currday = days(d);
                animdayidxs = curranimidxs(find(curralldays==currday));
                
                for i=1:length(animdayidxs)
                    c = animdayidxs(i); % The current index for Anim curranim and Day currday
                    [dayeptetcell] = allripplemod_idx(c,2:end);
                    day=dayeptetcell(1); ep=dayeptetcell(2); tet=dayeptetcell(3); cell=dayeptetcell(4);
                    if currday~=day
                        disp('Error - something wrong in day parsing - Mismatch');
                        keyboard;
                    else
                        cmd=sprintf('ripplemod{%d}{%d}{%d}{%d} = allripplemod(%d);',day, ep, tet, cell, c);
                        %cmd=sprintf('%sripplemod%02d{%d}{%d}{%d}{%d} = allripplemod(%d);', prefix, day, day, ep, tet, cell, c);
                        %eg. HParipplemod01{day}{ep}{tet}{cell} = allripplemod(c);
                        eval(cmd);
                    end
                end
                
                % Save current day for current animal
                currsavefile = sprintf('%sripplemod%02d',prefix,currday);
                currsavefile = [currsavedir currsavefile]
                save(currsavefile,'ripplemod');
            end % curr day
        end % curr anim
        
    end % end if
   
    
end % end gather data

length(find(allsig_shuf==1))
length(find(allsig_ttest==1))

%add to "super" variable
% super{reg}=modf;


% ------------------------------
% Plotting for individual cells
% ------------------------------

figdir = '/data25/sjadhav/HPExpt/Figures/RippleMod/Egs/'; saveg1=0;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
forppr=0;
if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end

figopt1=0;
if (figopt1)
    for i=1:length(allripplemod)
        curridx = allripplemod(i).index;
        currhist = allripplemod(i).hist;
        currraster = allripplemod(i).raster;
        sig_shuf = allripplemod(i).sig_shuf; % Sig or not
        sig_ttest = allripplemod(i).sig_ttest; % Sig or not
        modln_shuf = allripplemod(i).modln_shuf; % %tile value of shuffle
        modln_peak = allripplemod(i).modln_peak; % %tage peak change over baseline
        modln = allripplemod(i).modln; % %tage mean change over baseline
        
        rip_spkshist_cellsort_PFC = currhist; rip_spks_cellsort_PFC = currraster;
        day = curridx(2); tet = curridx(3); cell = curridx(4);
        if curridx(1)==1, prefix = 'HPa'; else, prefix = 'HPb'; end
        
        % 1) Raster
        % ----------
        figure; hold on; redimscreen_figforppt1;
        set(gcf, 'Position',[205 136 723 892]);
        %set(gcf, 'Position',[205 658 723 446]);
        subplot(2,1,1); hold on;
        spkcount = [];
        for c=1:length(rip_spks_cellsort_PFC)
            tmps = rip_spks_cellsort_PFC{c};
            plotraster(tmps,(length(rip_spks_cellsort_PFC)-c+1)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','b');
            % Get count of spikes in response window
            if ~isempty(tmps)
                subset_tmps = find(tmps>=rwin(1) & tmps<=rwin(2));
                spkcount = [spkcount; length(subset_tmps)];
            end
        end
        set(gca,'XLim',[-pret postt]);
        set(gca,'XTick',[-pret:200:postt],'XTickLabel',num2str([-pret:200:postt]'));
        xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
        ylabel('Ripple number','FontSize',yfont,'Fontweight','normal');
        set(gca,'YLim',[0 size(rip_spkshist_cellsort_PFC,1)]);
        % Plot Line at 0 ms and rwin
        ypts = 0:1:size(rip_spkshist_cellsort_PFC,1);
        xpts = 0*ones(size(ypts));
        plot(xpts , ypts, 'k--','Linewidth',2);
        % Plot lines at rwin
        xpts = rwin(1)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
        xpts = rwin(2)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
        xpts = bckwin(1)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
        xpts = bckwin(2)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
        
        title(sprintf('%s Day %d Tet %d Cell %d Nspkwin %d', prefix, day, tet, cell, sum(spkcount)),...
            'FontSize',tfont,'Fontweight','normal');
        if saveg1==1,
            figfile = [figdir,'RippleAlignRaster_',prefix,'_Day',num2str(day),'_Tet',num2str(PFCtet),'_Cell',num2str(PFCcell)];
            print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
        end
        
        
        % Hist
        % ----
        %figure; hold on; redimscreen_figforppt1;
        %set(gcf, 'Position',[205 136 723 446]);
        subplot(2,1,2); hold on;
        xaxis = -pret:binsize:postt;
        plot(xaxis,mean(rip_spkshist_cellsort_PFC),'Linewidth',3);
        plot(xaxis,mean(rip_spkshist_cellsort_PFC)+sem(rip_spkshist_cellsort_PFC),'b--','Linewidth',1);
        plot(xaxis,mean(rip_spkshist_cellsort_PFC)-sem(rip_spkshist_cellsort_PFC),'b--','Linewidth',1);
        
        set(gca,'XLim',[-pret postt]);
        xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
        ylabel('Fir rate (Hz)','FontSize',yfont,'Fontweight','normal');
        set(gca,'XTick',[-pret:200:postt],'XTickLabel',num2str([-pret:200:postt]'));
        ylow = min(mean(rip_spkshist_cellsort_PFC)-sem(rip_spkshist_cellsort_PFC));
        yhigh = max(mean(rip_spkshist_cellsort_PFC)+sem(rip_spkshist_cellsort_PFC));
        set(gca,'YLim',[ylow-0.1 yhigh+0.1]);
        ypts = ylow-0.1:0.1:yhigh+0.1;
        xpts = 0*ones(size(ypts));
        % Plot Line at 0 ms - Onset of stimulation
        plot(xpts , ypts, 'k--','Linewidth',2);
        % Plot lines at rwin and bckwi
        xpts = rwin(1)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
        xpts = rwin(2)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
        xpts = bckwin(1)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
        xpts = bckwin(2)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
        
        if sig_ttest ==1, str = '*'; else, str = ''; end
        if sig_shuf ==1, str_shuf = '*'; else, str_shuf = ''; end
        title(sprintf('%s Day%d Tet%d Cell%d: M %g%s Prc %g%s', prefix, day, tet, cell, roundn(modln_peak,-1),...
            str, roundn(modln_shuf,-2), str_shuf),'FontSize',tfont,'Fontweight','normal');
        if saveg1==1,
            figfile = [figdir,'RippleAlignHist_',prefix,'_Day',num2str(day),'_Tet',num2str(PFCtet),'_Cell',num2str(PFCcell)];
            print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
        end
        
        
        keyboard;
        
    end % end i
end % end if figopt









% ------------------
% Population Figures
% ------------------

% Population Data
% ----------------------------
cntsig = 0; cntnosig = 0;
allsigmodln_peak = []; allsigmodln_shuf = [];  allsighist = [];
allnosigmodln_peak = []; allnosigmodln_shuf = [];  allnosighist = [];

days = unique(alldays);
anim = unique(allanim);
ncells_days = zeros(length(days),1);
ncells_days_sig = zeros(length(days),1);

for i = 1:length(allripplemod)
    currday = allripplemod(i).days;
    ncells_days(currday) = ncells_days(currday)+1;
    
    if ( allsig_shuf(i) == 1)
        cntsig = cntsig+1;
        ncells_days_sig(currday) = ncells_days_sig(currday)+1;
        allsigmodln_peak(cntsig) = allripplemod(i).modln_peak;
        allsigmodln_shuf(cntsig) = allripplemod(i).modln_shuf;
        allsighist(cntsig,:) = mean(allripplemod(i).hist,1);
    else
        cntnosig = cntnosig+1;
        allnosigmodln_peak(cntnosig) = allripplemod(i).modln_peak;
        allnosigmodln_shuf(cntnosig) = allripplemod(i).modln_shuf;
        allnosighist(cntnosig,:) = mean(allripplemod(i).hist,1);
    end
end

% Normalize histogram by mean fr per row
allsignormhist = allsighist./repmat(max(allsighist,[],2),1,size(allsighist,2));
allnosignormhist = allnosighist./repmat(max(allnosighist,[],2),1,size(allnosighist,2));

forppr = 0;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Figures/ThetaMod/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end


if 1
    % 1) Norm histogram of all Sig and Non-Sig
    % ----------------------------------------
    
    figure; hold on; redimscreen_figforppt1;
    set(gcf, 'Position',[205 136 723 446]);
    xaxis = -pret:binsize:postt;
    plot(xaxis,mean(allsignormhist),'r','Linewidth',3);
    %plot(xaxis,mean(allsignormhist)+sem(allsignormhist),'r--','Linewidth',1);
    %plot(xaxis,mean(allsignormhist)-sem(allsignormhist),'r--','Linewidth',1);
    plot(xaxis,mean(allnosignormhist),'b','Linewidth',3);
    %plot(xaxis,mean(allnosignormhist)+sem(allnosignormhist),'b--','Linewidth',1);
    %plot(xaxis,mean(allnosignormhist)-sem(allnosignormhist),'b--','Linewidth',1);
    
    set(gca,'XLim',[-pret postt]);
    title(sprintf('Norm histogram for sig and nonsig ripple mod %s units',area));
    xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
    ylabel('Norm. Fir rate','FontSize',yfont,'Fontweight','normal');
    set(gca,'XTick',[-pret:200:postt],'XTickLabel',num2str([-pret:200:postt]'));
    ylow = min([mean(allsignormhist), mean(allnosignormhist)]);
    yhigh = max([mean(allsignormhist),mean(allnosignormhist)]);
    set(gca,'YLim',[ylow-0.05 yhigh+0.05]);
    ypts = ylow-0.2:0.1:yhigh+0.2;
    xpts = 0*ones(size(ypts));
    % Plot Line at 0 ms - Onset of stimulation
    plot(xpts , ypts, 'k--','Linewidth',2);
    % Plot lines at rwin and bckwi
    xpts = rwin(1)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
    xpts = rwin(2)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
    xpts = bckwin(1)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
    xpts = bckwin(2)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
    
    if savefig1==1,
        figfile = [figdir,area,'_Ripplemod_MeanHist'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
end


if 1
    % 2) No of sig modulated cells over days: %tage and number
    % ------------------------------------
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
    
    persig_days = 100*ncells_days_sig./ncells_days;
    plot(persig_days,[clr 'o'],'MarkerSize',18,'LineWidth',2);
    title(sprintf('%tage of sig modulated units'));
    xlabel(['Day'],'FontSize',xfont,'Fontweight','normal');
    ylabel(['Percentage of cells'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'YLim',[0 max(persig_days)+5]);
    
    if savefig1==1,
        figfile = [figdir,area,'_Ripplemod_PerSigDays'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
    
    Nsig_days = ncells_days_sig;
    plot(Nsig_days,[clr 'o'],'MarkerSize',18,'LineWidth',2);
    title(sprintf('No. of sig modulated units'));
    xlabel(['Day'],'FontSize',xfont,'Fontweight','normal');
    ylabel(['Number of cells'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'YLim',[0 max(Nsig_days)+2]);
    
    if savefig1==1,
        figfile = [figdir,area,'_Ripplemodmod_NSigDays'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
end























%
% if 0
%     % plot individual phase histogram of all units
%
%     norm = 1;
%
%     figure
%     titlestring=sprintf('%s %s phase hist of individual units // %s',animals{1},region,referencestring);
%     title(titlestring,'FontSize',14,'FontWeight','bold')
%     counter=1;
%     for k=1:length(caf.celloutput)
%         if counter==81
%             counter=1;
%             figure
%             titlestring=sprintf('%s %s phase hist of individual units // %s',animals{1},region,referencestring);
%             title(titlestring,'FontSize',14,'FontWeight','bold')
%         end
%         subplot(8,10,counter)
%         bins_plot = caf.celloutput(k).bins(1:(end-1));
%         bins_plot = bins_plot + (bins(2)-bins(1))/2;
%         phasehist = caf.celloutput(k).phasehist(1:(end-1));
%         phasehist_norm = phasehist/sum(phasehist);
%         if norm == 1
%             if size(phasehist_norm,1) < size(phasehist_norm,2)
%                 phasehist_norm = phasehist_norm';
%             end
%             %plot
%             h = bar([bins_plot bins_plot+2*pi],[phasehist_norm ; phasehist_norm],'histc');
%             title(num2str(caf.celloutput(k).index))
%             axis tight
%             ylim([0 .2])
%         else
%             if size(phasehist,1) < size(phasehist,2)
%                 phasehist = phasehist';
%             end
%             %plot
%             h = bar([bins_plot bins_plot+2*pi],[phasehist ; phasehist],'histc');
%             title(num2str(caf.celloutput(k).index),'FontSize',12,'FontWeight','bold')
%             axis tight
%             ylim([0 250])
%         end
%
%         set(h(1),'facecolor',clr)
%         set(h(1),'edgecolor',clr)
%
%         % plot guide lines
%         hold on
%         plot([pi,pi],[0 9999],'k--','LineWidth',1.5)
%         plot([-pi,-pi],[0 9999],'k--','LineWidth',1.5)
%         plot([3*pi,3*pi],[0 9999],'k--','LineWidth',1.5)
%         plot([0,0],[0 9999],'k:','LineWidth',1.5)
%         plot([2*pi,2*pi],[0 9999],'k:','LineWidth',1.5)
%
%         counter=counter+1;
%     end
% end
%



% % bar
% count = histc(allspikephases, bins);
% out = bar(bins, count, 'hist');
% set(out,'facecolor','k')
% title('aggregate theta modulation');
%
% % lineplot
% dischargeprob=count./sum(count);
% plot(bins(1:(end-1)),dischargeprob(1:(end-1)),'k','LineWidth',2);
%
% [m ph] = modulation(allspikephases);





