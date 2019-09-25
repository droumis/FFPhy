% Ver4 : Starting 10Feb2014 - Sync codes with everyone

% Ver2, Dec 2013 - Implement Min. NSpike condition for PFC cells. See Line 47 and Line 240

% Ripple modulation of cells, especilally PFC cells. Time filter version of sj_HPexpt_ripalign_singlecell_getrip4.
% Will call DFAsj_getripalign.m
% Also see DFSsj_plotthetamod.m and DFSsj_HPexpt_xcorrmeasures2. Will gather data like these

clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 1; % Figure Options - Individual cells

%savedir = '/data15/gideon/ProcessedData/';
savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/';
savefile = [savedir,'HPall_ripplemod_std3_speed4_ntet2_PFC_DR_testing']; area = 'PFC'; clr = 'b'; % PFC - low speed criterion
savefig1=0;
saveDR = 1;
figdir = '/mnt/data25/sjadhav/HPExpt/Figures_DR/';
gatherdata = 1; savegatherdata = 1;
figname = 'HPall_inhib';
gatherdatafile = [savedir 'HPall_ripplemod_std3_speed4_ntet2_PFC_postsleep_gather'];

% [y, m, d] = datevec(date);
% s
% if runscript==1
% %     val=6;savefile = [savedir,'HPa_ripplemod_PFC_DRTEST_std3_speed4_ntet2_',num2str(m),'-',num2str(d),'-',num2str(y)]; area = 'PFC'; clr = 'b'; % PFC - low speed criterion
%     
%     %val=7;savefile = [savedir,'HP_ripplemod_CA1_alldata_std3_speed4_ntet2_',num2str(m),'-',num2str(d),'-',num2str(y)]; area = 'CA1'; clr = 'r'; % CA1
% % else
% %     val=6;savefile = [savedir,'HP_ripplemod_PFC_alldata_std3_speed4_ntet2_2-12-2014']; area = 'PFC'; clr = 'b'; % PFC. 74 Mod, 137 UnMod, 211 total
% % %     %val=7;savefile = [savedir,'HP_ripplemod_CA1_alldata_std3_speed4_ntet2_2-12-2014']; area = 'CA1'; clr = 'r'; % CA1.  255 Mod, 75 UnMod: 330 total 
% % % 
% end


% Pre version4
% -------------
%val=1; savefile = [savedir,'HP_ripplemod_PFC_alldata-',num2str(m),'-',num2str(d),'-',num2str(y)]; area = 'PFC'; clr = 'b'; % PFC
%val=2; savefile = [savedir 'HP_ripplemod_CA1_alldata']; area = 'CA1';  clr = 'r';% CA1
%val=3;savefile = [savedir 'HP_ripplemod_PFC_alldata_speed']; area = 'PFC'; clr = 'b'; % PFC - low speed criterion
%val=4;savefile = [savedir 'HP_ripplemod_PFC_alldata_stdev5']; area = 'PFC'; clr = 'b'; % PFC - low speed criterion
%val=5;savefile = [savedir 'HP_ripplemod_PFC_alldata_singletrack']; area = 'PFC'; clr = 'b'; % PFC - single track



% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


%If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
%     animals = {'HPa','HPb','HPc'};
    animals = {'HPa'};
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    % Either Only do 1st w-track. 2 or 1 epochs per day
    % Or do Wtr1 and Wtr2, 2 epochs per day
%     if val==5
%         runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''ytr'')';
%     else
%         runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
%         runepochfilter = 'isequal($environment, ''presleep'')';
        runepochfilter = 'isequal($environment, ''postsleep'')';

%     end
%     %sleepepochfilter = 'isequal($type, ''sleep'')'; % Only pre and post sleep marked as sleep
    
    % Cell filter
    % -----------
%     switch val
%         
%         case 6
            cellfilter = 'strcmp($area, ''PFC'') && ($numspikes > 100)'; % PFC cells with spiking criterion
       
%         case 7
%             %cellfilter = ' (strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && ~strcmp($tag2, ''CA1Int'') && ~strcmp($tag2, ''iCA1Int'')';
%             cellfilter = '(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && ($meanrate < 7)';
%             
%        % Pre version4
%        % -------------
%         case 1
%             cellfilter = 'strcmp($area, ''PFC'') && ($numspikes > 100)'; % PFC cells with spiking criterion
%             %cellfilter = 'strcmp($area, ''PFC'')'; % This includes all, including silent cells
%         case 2
%             %cellfilter = ' (strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && ~strcmp($tag2, ''CA1Int'') && ~strcmp($tag2, ''iCA1Int'')';
%             cellfilter = '(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && ($meanrate < 7)';
%         case 3
%             cellfilter = 'strcmp($area, ''PFC'') && ($numspikes > 100)'; % PFC cells with spiking criterion
%             %cellfilter = 'strcmp($area, ''PFC'')'; % This includes all, including silent cell
%         case 4
%             cellfilter = 'strcmp($area, ''PFC'') && ($numspikes > 100)'; % PFC cells with spiking criterion
%             %cellfilter = 'strcmp($area, ''PFC'')'; % This includes all, including silent cells
%         case 5
%             cellfilter = 'strcmp($area, ''PFC'') && ($numspikes > 100)'; % PFC cells with spiking criterion
%             %cellfilter = 'strcmp($area, ''PFC'')'; % This includes all, including silent cells
            
        
%     end
    
    % cellfilter = 'strcmp($tag2, ''PFC'')';
    %cellfilter = 'strcmp($tag2, ''CA1Pyr'') && ($numspikes > 100)'; % This includes all.
    % For more exclusive choosing, use $tag. Take care of number of spikes while gathering data
    
    % Time filter
    % -----------
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % The following are ripple time filter options. Instead of using time filters, within the function,
    % call getripples using the riptetfilter option to get the riptimes. Since you want a vector of riptimes.
    % You can also use inter-ripple-interval of 1 sec within the function, etc.
    
    % a) Using getriptimes: generates 1 vector with 1's for ripple times for all tetrodes.
    % Should be similar to getripples/ getripplees_direct
    
    %     timefilterrun_rip = {{'DFTFsj_getvelpos', '(($absvel <= 5))'},...
    %         {'DFTFsj_getriptimes','($nripples > 0)','tetfilter',riptetfilter,'minthresh',3}};
    
    % b) getripples. This is similar to getripltimes, but gives start and end of each ripple time
    % instead of a time filter. Also has an additional condition of at least 50 ms ripple.
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
        modf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter, 'cells', cellfilter, 'iterator', iterator);
%     modf = createfilter('animal',animals,'epochs',runepochfilter, 'cells',...
%         cellfilter, 'iterator', iterator);
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
%     
%    switch val      
%        
%        % Version 4
%        % ---------
%        case {6,7}
        modf = setfilterfunction(modf,'DFAsj_getripalignspiking4',{'spikes', 'ripples', 'tetinfo', 'pos'},'dospeed',1,'lowsp_thrs',4,'minrip',2,'minstd',3); % Default stdev is 3
%         
%        % Pre version4
%        % -------------
%        case 1
%         modf = setfilterfunction(modf,'DFAsj_getripalignspiking',{'spikes', 'ripples', 'tetinfo', 'pos'}); % Default stdev is 3
%        case 2
%         modf = setfilterfunction(modf,'DFAsj_getripalignspiking',{'spikes', 'ripples', 'tetinfo', 'pos'}); % 
%        case 3
%         modf = setfilterfunction(modf,'DFAsj_getripalignspiking',{'spikes', 'ripples', 'tetinfo', 'pos'},'dospeed',1,'lowsp_thrs',4); %
%        case 4
%         modf = setfilterfunction(modf,'DFAsj_getripalignspiking',{'spikes', 'ripples', 'tetinfo', 'pos'},'minstd',5); %
%        case 5
%         modf = setfilterfunction(modf,'DFAsj_getripalignspiking',{'spikes', 'ripples', 'tetinfo', 'pos'}); % Default stdev is 3
%        
%     end
    % Going to call getripples_tetinfo within function
    
    % Run analysis
    % ------------
    modf = runfilter(modf);
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata saveDR gatherdata savegatherdata figname savedir gatherdatafile
        save(savefile);
    end
    
else
    
    load(savefile);
    
end % end runscript

if ~exist('savedata')
    return
end


% -------------------------  Filter Format Done -------------------------



% ----------------------------------
% Whether to gather data or to load previously gathered data
% --------------------------------------------------------------------

% [y, m, d] = datevec(date);


% % Version 4 onwards
% % Filename of gatherdata. If generating, put current date. If not, then load a file from previous date.
% % -----------------------------------------------------------------------------------------------------
% switch val
%     case 6
%         if gatherdata
            
%         else
%             gatherdatafile = [savedir 'HP_ripplemod_PFC_alldata_std3_speed4_ntet2_Nspk50_gather_2-12-2014']; % PFC cells to Hipp ripples - low speed criterion
%         end   
%     case 7
%         if gatherdata
%             gatherdatafile = [savedir 'HP_ripplemod_CA1_alldata_std3_speed4_ntet2_Nspk50_gather_',num2str(m),'-',num2str(d),'-',num2str(y)], % PFC cells to Hipp ripples - low speed criterion
%         else
%             gatherdatafile = [savedir 'HP_ripplemod_CA1_alldata_std3_speed4_ntet2_Nspk50_gather_2-12-2014'], % PFC cells to Hipp ripples - low speed criterion
%         end 
%         
%         % Pre version4
%         % -------------
%     case 1
%         gatherdatafile = [savedir 'HP_ripplemod_PFC_alldata_gather_var2']; % PFC cells to Hipp ripples
%     case 2
%         gatherdatafile = [savedir 'HP_ripplemod_CA1_alldata_gather_var']; % CA1 cells to Hipp ripples
%     case 3
%         gatherdatafile = [savedir 'HP_ripplemod_PFC_alldata_speed_gather_var']; % PFC cells to Hipp ripples - low speed criterion
%     case 4
%         gatherdatafile = [savedir 'HP_ripplemod_PFC_alldata_stdev5_gather_var_Nspk50']; % PFC cells to Hipp ripples - low speed criterion
%     case 5
%         gatherdatafile = [savedir 'HP_ripplemod_PFC_alldata_singletrack_gather_var']; % PFC cells to Hipp ripples - low speed criterion
%         
% end

if gatherdata
    
    % Parameters if any
    % -----------------
    
    % -------------------------------------------------------------
    
    cnt=0; % Count how many cells will be kept based on nspikes in output: >0
    allanimindex=[]; alldataraster=[]; alldatahist = []; all_Nspk=[];
    
    for an = 1:length(modf)
        for i=1:length(modf(an).output{1})
            % Check for empty output - If Cell defined in rpoch and Nspks in ripple response wndow > 0
            if (modf(an).output{1}(i).Nspikes > 0)
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
                % Nspikes summed across trials in response and bckgnd window
                all_Nspk_resp(cnt) = sum(modf(an).output{1}(i).trialResps);
                all_Nspk_bck(cnt) = sum(modf(an).output{1}(i).trialResps_bck);
                % Properties
                allcellfr(cnt) = modf(an).output{1}(i).cellfr;
                
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
            end
        end
        
    end
    
    % Consolidate single cells across epochs. Multiple methods: see also DFSsj_getcellinfo and DFSsj_xcorrmeasures2
    % ----------------------------------------------------------------------------
    
    allripplemod = struct;
    
    % Method 1
    % ---------------------------------------------
    %     cntcells=0;
    %     animdaytetcell = unique(allanimindex(:,[1 2 4 5]),'rows'); % Collapse across epochs
    %     for ind = 1:size(animdaytetcell,1)
    %         currhist=[]; currraster = [];
    %         for a = 1:size(allanimindex,1)
    %             if animdaytetcell(ind,:)==allanimindex(a,[1 2 4 5]) % Epoch match
    %                 currhist = [currhist; alldatahist{a}];
    %                 currraster = [currraster; alldataraster{a}];
    %             end
    %         end
    %         cntcells = cntcells + 1;
    %         allripplemod_idx(cntcells,:)=animdaytetcell(ind,:);
    %         allripplemod(cntcells).index=animdaytetcell(ind,:);
    %         allripplemod(cntcells).hist=currhist;
    %         allripplemod(cntcells).raster=currraster;
    %         allripplemod(cntcells).Nspk=Nspk;
    %     end
    
    
    % Method 2
    % ---------
    dummyindex=allanimindex;  % all anim-day-epoch-tet-cell indices
    cntcells=0;
    for i=1:length(alldatahist)
        animdaytetcell=allanimindex(i,[1 2 4 5]);
        ind=[];
        while rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))~=0          % collect all rows (epochs)
            ind = [ind rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))];        % finds the first matching row
            dummyindex(rowfind(animdaytetcell,dummyindex(:,[1 2 4 5])),:)=[0 0 0 0 0]; % after adding index, remove the corresponding row
            % so you could find the next one
        end
        
        % Gather everything for the current cell across epochs
        currhist=[]; currraster=[]; currNspk=0; currNspk_resp=0; currNspk_bck=0; curr_cellfr=[];
        currhist_rdm=[]; currraster_rdm=[];
        currtrialResps=[]; currtrialResps_rdm=[]; currtrialResps_bck=[];
        for r=ind
            currNspk = currNspk + all_Nspk(r);
            currNspk_resp = currNspk_resp + all_Nspk_resp(r);
            currNspk_bck = currNspk_bck + all_Nspk_bck(r);
            currhist = [currhist; alldatahist{r}];
            currraster = [currraster, alldataraster{r}];
            currhist_rdm = [currhist_rdm; alldatahist_rdm{r}];
            currraster_rdm = [currraster_rdm, alldataraster_rdm{r}];
            currtrialResps = [currtrialResps, alldatatrialResps{r}];
            currtrialResps_rdm = [currtrialResps_rdm, alldatatrialResps_rdm{r}];
            currtrialResps_bck = [currtrialResps_bck, alldatatrialResps_bck{r}];
            curr_cellfr = [curr_cellfr; allcellfr(r)];
        end
        
        % Condition for Nspk. Version 1 had a min of 50 for entire window. Increase it to 100,
        % and can also add a condition for spikes in (resp+bck) window. Need to try a few values
        if (currNspk >= 50)
            %if ((currNspk_resp+currNspk_bck) >= 40)
            %if (currNspk >= 100) && ((currNspk_resp+currNspk_bck) >= 40)
            cntcells = cntcells + 1;
            
            % DONT SHIFT ANY DAYS
            % For Ndl-GIdeon;s animal, shift days
%             if animdaytetcell(1)==4
%                 animdaytetcell(2)=animdaytetcell(2)-7; % Day starts from no. 8
%             end
            allripplemod_idx(cntcells,:)=animdaytetcell;
            allripplemod(cntcells).index=animdaytetcell;
            allripplemod(cntcells).hist=currhist*(1000/binsize); % Convert to firing rate in Hz
            allripplemod(cntcells).raster=currraster;
            allripplemod(cntcells).Nspk=currNspk;
            allripplemod(cntcells).hist_rdm=currhist_rdm*(1000/binsize); % Convert to firing rate in Hz
            allripplemod(cntcells).raster_rdm=currraster_rdm;
            % Trial Resps
            allripplemod(cntcells).trialResps = currtrialResps';
            allripplemod(cntcells).trialResps_rdm = currtrialResps_rdm';
            allripplemod(cntcells).trialResps_bck = currtrialResps_bck';
            % Properties
            allripplemod(cntcells).cellfr = nanmean(curr_cellfr);
        end
    end
    
    
    
    % Calculations/ Stats. Stats between response and bck
    % Similar to ...getrip4
    % -----------------------------------------------------------
    for i=1:cntcells
        
        curr_cellfr = allripplemod(i).cellfr;
        currhist = allripplemod(i).hist; %currraster = allripplemod(i).raster;
        currhist_rdm = allripplemod(i).hist_rdm; %currraster_rdm = allripplemod(i).raster_rdm;
        
        % Get the bckresp, rresp and rdmresp again - using firing rate matrices
        rresp = currhist(:,bins_resp);
        bckresp = currhist(:,bins_bck);
        rresp_rdm = currhist_rdm(:,bins_resp);
        
        
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
        avgresprdm_trial = mean(rresp_rdm,2); % Mean for each random ripple
        %trialRespsh_rdm = sum(rresp_rdm,2); % Nspikes/ trial in resp window
        avgrdmhist = mean(rresp_rdm,1); % Avg resp histogram
        mean_rdmresp = mean(mean(rresp_rdm)); % Single value
        
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
        
        % Get significance by comparing Dm to Dshuf. One-tailed test
        % --------------------------------------------------------------
        clear maxpeak minpeak 
        peakindex = [];
        if Dm>=0
            pshuf = length(find(Dshuf>Dm))/nshuffles;
            if Dm > prctile(Dshuf,95)
                sig_shuf = 1;
            end
            type = 'exc'; peakresp = max(avgresphist); % Peak in response histogram
            peakresp_rdm = max(avgrdmhist);
            %DR. save the bin of the peak of the average across trials.. if
            %this is an excited cell
            [maxpeak peakindex] = max(avgresphist);
            
        else
            pshuf = length(find(Dshuf<Dm))/nshuffles;
            if Dm < prctile(Dshuf,5)
                sig_shuf = 1;
            end
            type = 'inh'; peakresp = min(avgresphist); % Trough in response histogram
            peakresp_rdm = min(avgrdmhist);
            [minpeak peakindex] = min(avgresphist);
        end
        % Get the p-value of shuffle. The modulation index for shuffle will be the prctile value
        % ----------------------------
        modln_shuf = 100 - (pshuf*100); %eg p=0.05 => prctile is 95
        % Get %tage change over baseline and peak/trough: Save with sign- +ve or -ve modln. Not abs
        % ----------------------------------------------
        modln_raw = Dm;
        modln = 100*(Dm)/mean(avgbckresp_trial); % Mean %tage change above/below baseline
        peakchange = peakresp - mean(avgbckhist); % mean(avgbckhist) is same as mean(avgbckresp_trial)
        modln_peak = 100*(peakchange)/mean(avgbckresp_trial); % Peak/trough %tage change above/below baseline
        modln_div = 100*peakresp/mean(avgbckresp_trial);
        
        
        % Get same values for random times as well
        Dm_rdm = mean(avgresprdm_trial) - mean(avgbckresp_trial);
        modln_rdm = 100*(Dm_rdm)/mean(avgbckresp_trial); % Mean %tage change above/below baseline
        peakchange_rdm = peakresp_rdm - mean(avgbckhist); % mean(avgbckhist) is same as mean(avgbckresp_trial)
        modln_peak_rdm = 100*(peakchange_rdm)/mean(avgbckresp_trial); % Peak/trough %tage change above/below baseline
       
        % Save
        % -----
        allripplemod(i).peakindex = peakindex; %DR
        allripplemod(i).Dm = Dm;
        allripplemod(i).pshuf = pshuf; allripplemod(i).p = pshuf;
        allripplemod(i).sig_shuf = sig_shuf; allripplemod(i).h = sig_shuf; % Sig or not
        allripplemod(i).sig_ttest = sig_ttest; % Sig or not
        allripplemod(i).modln_shuf = modln_shuf ; % %tile value of shuffle
        allripplemod(i).modln_peak = modln_peak; % %tage peak change over baseline
        allripplemod(i).modln_div = modln_div;
        allripplemod(i).modln = modln; % %tage mean change over baseline
        allripplemod(i).modln_raw = modln_raw; % Raw change in firing rate. The statistic
        allripplemod(i).type = type; % exc or inh
        allripplemod(i).anim = allripplemod(i).index(1); allanim(i) = allripplemod(i).index(1);
        allripplemod(i).days = allripplemod(i).index(2); alldays(i) = allripplemod(i).index(2);
        
        % Mean resp from histogram in respective window
        allripplemod(i).avghisttrialResps = avgresp_trial;
        allripplemod(i).avghisttrialResps_rdm = avgresprdm_trial;
        allripplemod(i).avghisttrialResps_bck = avgrespbck_trial;
        
        allsig_shuf(i) = sig_shuf;
        allsig_ttest(i) = sig_ttest;
        
        % Properties
        allripplemod(i).cellfr = curr_cellfr;
        allripplemod(i).Nrip = length(avgresp_trial);
        
        % Rdm resp modulation
        allripplemod(i).modln_peak_rdm = modln_peak_rdm; % %tage peak change over baseline
        allripplemod(i).modln_rdm = modln_rdm; % %tage mean change over baseline
        
    end
    

    
    
%     
%     
%     
%     % GIDEON
%     %varRange=[400:700];
%     varRange=[500:700];
%     bckRange=[1:200];
%     for ii=1:cntcells
%         ii;
%         curRast=allripplemod(ii).raster;
%         numRips=length(curRast);
%         % creating a raster matrix of 0/1 with 1ms resolution 
%         % it goes from -550 t0 +550ms relative to ripples
%         curRastMat=zeros(numRips,1100);
%         for i=1:numRips,curRastMat(i,round(curRast{i}+551))=1;end    
%         
%         % --- Old, no longer nec ---
%         meanRespt=smooth(mean(curRastMat(:,50:end-50)),50);
%         respVar=var(meanRespt(varRange)); bckVar=var(meanRespt(bckRange));
%         var_respbck = respVar/bckVar; % Var in win divided by var in bck
%         bckVarShufs=[]; respbckratioShufs=[];
%         % --- Old, no longer nec ---
%         respVarShufs=[];
%         numRuns=1000;
%         allShufMeans=[];
%         % shuffling each trial in the raster separately in time,
%         % cyclically. Doing that numRuns times (currently 1000).
%         for runs=1:numRuns
%             
%             curRastMatShuf=zeros(numRips,1100);            
%             for i=1:numRips
%                 shiftVal=round(rand(1)*1000);
%                 curRastMatShuf(i,1+mod(round(curRast{i}+551+shiftVal),1100))=1;
%             end
%             
%             % shuffled "psth"
%             meanRespShuf=smooth(mean(curRastMatShuf(:,50:end-50)),50);
%             allShufMeans=[allShufMeans; meanRespShuf'];
%             %             plot(smooth(mean(curRastMatShuf(:,50:end-50)),50),'r')
%             %             hold on;
%             
%             % for each shuffled psth, calculate the variance (mean squared
%             % distance from mean)
%             respVarShufs=[respVarShufs var(meanRespShuf(varRange))];
%             % --- Old, no longer nec ---
%             bckVarShufs=[bckVarShufs var(meanRespShuf(bckRange))];
%             respbckratioShufs=[respbckratioShufs var(meanRespShuf(varRange))./var(meanRespShuf(bckRange))];
%             % --- Old, no longer nec ---
%         end
%         
%         meanResp=smooth(mean(curRastMat(:,50:end-50)),50);
%         meanRespRange=meanResp(varRange);
%         meanRespShuf=(mean(allShufMeans));
%         meanRespShufRange= meanRespShuf(varRange);
%         
%         %old measure for response: the variance of the "real" psth
%         respVar=var(meanResp(varRange));
%         %new measure: instead of mean squared distance from its own mean,
%         %mean squared distance from the mean shuffled psth's
%         respVar2=mean((meanRespRange-meanRespShufRange').^2);
%         
%         rasterShufP=1-sum(respVar>respVarShufs)/numRuns;
%         rasterShufP2=1-sum(respVar2>respVarShufs)/numRuns;      
%         
%         % --- Old, no longer nec ---
%         rasterShufP_respbck=1-sum(var_respbck>respbckratioShufs)/numRuns;
%         %          plot(smooth(mean(curRastMat(:,50:end-50)),50),'k','linewidth',2)
%         % --- Old, no longer nec ---
%         
%         allripplemod(ii).rasterShufP=rasterShufP;
%         allripplemod(ii).rasterShufP2=rasterShufP2;
%         
%         varRespAmp=respVar/mean(respVarShufs);
%         varRespAmp2=respVar2/mean(respVarShufs);
%         
%         allripplemod(ii).varRespAmp=varRespAmp;
%         allripplemod(ii).varRespAmp2=varRespAmp2;
%         
%         % --- Old, no longer nec ---
%         % Ratio of var in Resp / ratio of var in bck: how much did it change in resp relative to bck
%         varBckAmp=bckVar/mean(bckVarShufs);
%         var_changerespbck = varRespAmp/varBckAmp;
%         %Save
%         allripplemod(ii).var_respbck=var_respbck; allripplemod(ii).rasterShufP_respbck=rasterShufP_respbck;
%         allripplemod(ii).var_changerespbck=var_changerespbck;
%         % --- Old, no longer nec ---
%     end
%     
%     
    
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
end % end gather data

%length(find(allsig_shuf==1))
%length(find(allsig_ttest==1))

%______________________________________________________________________________________________________________________
%DR find latency to peak
day1peakbins = []; day2peakbins = []; day3peakbins = []; day4peakbins = []; day5peakbins = []; day6peakbins = []; day7peakbins = []; day8peakbins = []; peakcellsall = [];
for i = 1:length(allripplemod);
    %     if (~isempty(allripplemod(i).peakindex)) && (strcmp(allripplemod(i).type,'exc') || strcmp(allripplemod(i).type,'inh')) && (allripplemod(i).sig_shuf == 1);
    %         if (~isempty(allripplemod(i).peakindex)) && (strcmp(allripplemod(i).type,'exc')) && (allripplemod(i).sig_shuf == 1);
    if (~isempty(allripplemod(i).peakindex)) && (strcmp(allripplemod(i).type,'exc'));
        if (allripplemod(i).days == 1);
            day1peakbins = [day1peakbins; (allripplemod(i).peakindex)];
            peakcellsall = [peakcellsall; 1 i];
        elseif (allripplemod(i).days == 2);
            day2peakbins = [day2peakbins; (allripplemod(i).peakindex)];
            peakcellsall = [peakcellsall; 2 i];
        elseif (allripplemod(i).days == 3);
            day3peakbins = [day3peakbins; (allripplemod(i).peakindex)];
            peakcellsall = [peakcellsall; 3 i];
        elseif (allripplemod(i).days == 4);
            day4peakbins = [day4peakbins; (allripplemod(i).peakindex)];
            peakcellsall = [peakcellsall; 4 i];
        elseif (allripplemod(i).days == 5);
            day5peakbins = [day5peakbins; (allripplemod(i).peakindex)];
            peakcellsall = [peakcellsall; 5 i];
        elseif (allripplemod(i).days == 6);
            day6peakbins = [day6peakbins; (allripplemod(i).peakindex)];
            peakcellsall = [peakcellsall; 6 i];
        elseif (allripplemod(i).days == 7);
            day7peakbins = [day7peakbins; (allripplemod(i).peakindex)];
            peakcellsall = [peakcellsall; 7 i];
        elseif (allripplemod(i).days == 8);
            day8peakbins = [day8peakbins; (allripplemod(i).peakindex)];
            peakcellsall = [peakcellsall; 8 i];
        end
    end
end

figure;
bar([1:8], [mean(day1peakbins) mean(day2peakbins) mean(day3peakbins) mean(day4peakbins) mean(day5peakbins) mean(day6peakbins) mean(day7peakbins) mean(day8peakbins)])
hold on;
scatter(ones(length(day1peakbins),1), day1peakbins,'r')
hold on;
scatter(ones(length(day2peakbins),1)*2, day2peakbins,'r')
hold on;
scatter(ones(length(day3peakbins),1)*3, day3peakbins,'r')
hold on;
scatter(ones(length(day4peakbins),1)*4, day4peakbins,'r')
hold on;
scatter(ones(length(day5peakbins),1)*5, day5peakbins,'r')
hold on;
scatter(ones(length(day6peakbins),1)*6, day6peakbins,'r')
hold on;
scatter(ones(length(day7peakbins),1)*7, day7peakbins,'r')
hold on;
scatter(ones(length(day8peakbins),1)*8, day8peakbins,'r')
%  bar(1, mean(posRAllAn), 'FaceColor', [0 .7 .93], 'EdgeColor', 'none')

 if saveDR==1;
     mkdir('/mnt/data25/sjadhav/HPExpt/Figures_DR//RippleEfficacy/');
                figfile = [figdir,'/RippleEfficacy/',figname];
                print('-dpng', figfile);
 end
    
%______________________________________________________________________________________________________________________
% return
% ------------------------------
% Plotting for individual cells
% ------------------------------

figdir = '/data25/sjadhav/HPExpt/Figures/AllRippleMod/'; saveg1=0;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
forppr=0;
if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',40);
    tfont = 40;
    xfont = 40;
    yfont = 40;
end

figopt1=1; tmpcnt=0;
if (figopt1)
    for i=1:cntcells
        curridx = allripplemod(i).index;
        %if curridx(1)==4
        currhist = allripplemod(i).hist;
        currraster = allripplemod(i).raster;
        sig_shuf = allripplemod(i).sig_shuf; % Sig or not
        sig_ttest = allripplemod(i).sig_ttest; % Sig or not
        modln_shuf = allripplemod(i).modln_shuf; % %tile value of shuffle
        modln_peak = allripplemod(i).modln_peak; % %tage peak change over baseline
        modln = allripplemod(i).modln; % %tage mean change over baseline
        
        currNspk = allripplemod(i).Nspk;
        
        % New
%         modln_var = allripplemod(i).varRespAmp;
%         p_var = allripplemod(i).rasterShufP;
%         modln_var2 = allripplemod(i).varRespAmp2;
%         p_var2 = allripplemod(i).rasterShufP2;
        sigvar=1; 
%         if p_var<0.05, sigvar=1; end
        sigvar2=1; 
%         if p_var2<0.05, sigvar2=1; end
%         cellfr = allripplemod(i).cellfr;
%         
        rip_spkshist_cellsort_PFC = currhist; rip_spks_cellsort_PFC = currraster;
        day = curridx(2); tet = curridx(3); cell = curridx(4);
        switch curridx(1)
            case 1
                prefix = 'HPa';
            case 2
                prefix = 'HPb';
            case 3
                prefix = 'HPc';
            case 4
                prefix = 'Ndl';
        end
        
        % To control plotting
        %if sigvar==0 && sig_shuf==1
        %if curridx(1)==4 && curridx(2)==10 && curridx(3)==17 && curridx(4)==4
        if sigvar2==1
            
            tmpcnt = tmpcnt+1,
            % 1) Raster
            % ----------
            figure; hold on; redimscreen_figforppt1;
            set(gcf,'Position',[100 130 1000 950]);
            %set(gcf, 'Position',[205 658 723 446]);
            subplot(2,1,1); hold on;
            spkcount = [];
            for c=1:length(rip_spks_cellsort_PFC)
                tmps = rip_spks_cellsort_PFC{c};
                
                %plotraster(tmps,(length(rip_spks_cellsort_PFC)-c+1)*ones(size(tmps)),0.8,[],'LineWidth',2,'Color','k');
                plot(tmps,(length(rip_spks_cellsort_PFC)-c+1)*ones(size(tmps)),'k.','MarkerSize',16);
                %sj_plotraster(tmps,(length(rip_spks_cellsort_PFC)-c+1)*ones(size(tmps)),0.8,'k');
                
                % Get count of spikes in response window
                if ~isempty(tmps)
                    subset_tmps = find(tmps>=rwin(1) & tmps<=rwin(2));
                    spkcount = [spkcount; length(subset_tmps)];
                end
            end
            set(gca,'XLim',[-pret postt]);
            set(gca,'XTick',[])
            %set(gca,'XTick',[-pret:200:postt],'XTickLabel',num2str([-pret:200:postt]'));
            xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
            ylabel('SWR number','FontSize',yfont,'Fontweight','normal');
            set(gca,'YLim',[0 size(rip_spkshist_cellsort_PFC,1)]);
            % Plot Line at 0 ms and rwin
            ypts = 0:1:size(rip_spkshist_cellsort_PFC,1);
            xpts = 0*ones(size(ypts));
            plot(xpts , ypts, 'k--','Linewidth',2);
            % Plot lines at rwin
            %             xpts = rwin(1)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
            %             xpts = rwin(2)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
            %             xpts = bckwin(1)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
            %             xpts = bckwin(2)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
            %
%             title(sprintf('%s Day %d Tet %d Cell %d Nspkwin %d FR %g', prefix, day, tet, cell, sum(spkcount), roundn(cellfr,-1)),...
%                 'FontSize',tfont,'Fontweight','normal');
            if saveg1==1,
                figfile = [figdir,'RippleAlignRaster_',prefix,'_Day',num2str(day),'_Tet',num2str(tet),'_Cell',num2str(cell)];
                print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
            end
            %close(1);
            
            % Hist
            % ----
            %figure; hold on; redimscreen_figforppt1;
            %set(gcf, 'Position',[205 136 723 446]);
            subplot(2,1,2); hold on;
            xaxis = -pret:binsize:postt;
            plot(xaxis,mean(rip_spkshist_cellsort_PFC),'k-','Linewidth',4);
            %plot(xaxis,mean(rip_spkshist_cellsort_PFC)+sem(rip_spkshist_cellsort_PFC),'b--','Linewidth',1);
            %plot(xaxis,mean(rip_spkshist_cellsort_PFC)-sem(rip_spkshist_cellsort_PFC),'b--','Linewidth',1);
            
            set(gca,'XLim',[-pret postt]);
            xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal');
            ylabel('Firing rate (Hz)','FontSize',yfont,'Fontweight','normal');
            %set(gca,'XTick',[-pret:200:postt],'XTickLabel',num2str([-pret:200:postt]'));
            set(gca,'XTick',[-500,-250,0,250,500],'XTickLabel',num2str([-500,-250,0,250,500]'));
            ylow = min(mean(rip_spkshist_cellsort_PFC)-sem(rip_spkshist_cellsort_PFC));
            yhigh = max(mean(rip_spkshist_cellsort_PFC)+sem(rip_spkshist_cellsort_PFC));
            set(gca,'YLim',[ylow-0.1 yhigh+0.1]);
            ypts = ylow-0.1:0.1:yhigh+0.1;
            xpts = 0*ones(size(ypts));
            % Plot Line at 0 ms - Onset of stimulation
            plot(xpts , ypts, 'k--','Linewidth',2);
            % Plot lines at rwin and bckwi
            %             xpts = rwin(1)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
            %             xpts = rwin(2)*ones(size(ypts)); plot(xpts , ypts, 'k--','Linewidth',1);
            %             xpts = bckwin(1)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
            %             xpts = bckwin(2)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
            
            if sig_ttest ==1, str = '*'; else, str = ''; end
            if sig_shuf ==1, str_shuf = '*'; else, str_shuf = ''; end
            if sigvar ==1, str_var = '*'; else, str_var = ''; end
            if sigvar2 ==1, str_var2 = '*'; else, str_var2 = ''; end
            %         title(sprintf('%s Day%d Tet%d Cell%d: M %g%s Prc %g%s', prefix, day, tet, cell, roundn(modln_peak,-1),...
            %             str, roundn(modln_shuf,-2), str_shuf),'FontSize',tfont,'Fontweight','normal');
            title(sprintf('M %g%s Mvar %g%s Mvar2 %g%s', roundn(modln_peak,-1),...
                str_shuf, roundn(modln_var,-1), str_var, roundn(modln_var2,-1), str_var2),'FontSize',tfont,'Fontweight','normal');
            if saveg1==1,
                figfile = [figdir,'RippleAlignHist_',prefix,'_Day',num2str(day),'_Tet',num2str(tet),'_Cell',num2str(cell)];
                print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
            end
            
            keyboard;
            
            close all
            
        end % if curridx and sig_shuf
        
        
        
    end % end cntcells
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
% Force days to got from 1-10. For Ndl, you will push days by 7
days = 1:10;
ncells_days = zeros(length(days),1);
ncells_days_sig = zeros(length(days),1);

for i = 1:length(allripplemod)
    curranim = allripplemod(i).anim;
    currday = allripplemod(i).days;
    if curranim==4
        currday = currday-7; % For Ndl, days start from 8
    end
    ncells_days(currday) = ncells_days(currday)+1;
    
    % Change to var measure
    p_var2 = allripplemod(i).rasterShufP2;
    sigvar=0; if p_var2<0.05, sigvar=1; end
    varRespAmp2=[];    varRespAmp2=allripplemod(i).varRespAmp2;
    
    
    %if ( allsig_shuf(i) == 1)
    if ( sigvar == 1)
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

figdir = '/data25/sjadhav/HPExpt/Figures/31Oct/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',40);
    tfont = 40;
    xfont = 40;
    yfont = 40;
end


if 1
    % 1) Norm histogram of all Sig and Non-Sig
    % ----------------------------------------
    
    figure; hold on; redimscreen_figforppt1;
    set(gcf, 'Position',[205 136 723 446]);
    xaxis = -pret:binsize:postt;
    plot(xaxis,mean(allsignormhist),'r','Linewidth',5);
    plot(xaxis,mean(allnosignormhist),'b','Linewidth',5);
    legend('Ripple Mod Units','Non-Mod Units');
    plot(xaxis,mean(allsignormhist)+sem(allsignormhist),'r--','Linewidth',1);
    plot(xaxis,mean(allsignormhist)-sem(allsignormhist),'r--','Linewidth',1);
    plot(xaxis,mean(allnosignormhist)+sem(allnosignormhist),'b--','Linewidth',1);
    plot(xaxis,mean(allnosignormhist)-sem(allnosignormhist),'b--','Linewidth',1);
    
    set(gca,'XLim',[-pret postt]);
    title(sprintf('%s Run: Norm histogram for ripple aligned response',area),'FontSize',tfont);
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
    xpts = bckwin(1)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);%
    xpts = bckwin(2)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
    
    %%Nsig = length(find(allsig_shuf==1));
    %%Nnosig = length(find(allsig_shuf==0));
    text(-pret+50, 0.85, sprintf('Mod: %d',cntsig),'FontSize',xfont,'Fontweight','normal');
    text(-pret+50, 0.8, sprintf('UnMod: %d',cntnosig),'FontSize',xfont,'Fontweight','normal');
    
    figfile = [figdir,area,'_Run_RippleAlign_PoplnHist']
    if savefig1==1,
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
end


if 1
    % 2) No of sig modulated cells over days: %tage and number
    % ------------------------------------
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
    
    persig_days = 100*ncells_days_sig./ncells_days;
    plot(persig_days,['ro'],'MarkerSize',30,'LineWidth',2);
    title(sprintf('%s Run: Fraction Ripple mod cells',area),'FontSize',tfont);
    xlabel(['Day'],'FontSize',xfont,'Fontweight','normal');
    ylabel(['Percentage of cells'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'XLim',[0.5 length(persig_days)+0.5])
    set(gca,'YLim',[0 max(persig_days)+5]);
    
    figfile = [figdir,area,'_Run_RippleAlign_PertageSigDays']
    if savefig1==1,
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
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
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
    
    N_days = ncells_days;
    plot(N_days,[clr 'o'],'MarkerSize',18,'LineWidth',2);
    title(sprintf('No. of total units'));
    xlabel(['Day'],'FontSize',xfont,'Fontweight','normal');
    ylabel(['Number of total cells'],'FontSize',yfont,'Fontweight','normal');
    %set(gca,'YLim',[0 max(Nsig_days)+2]);
    
    if savefig1==1,
        figfile = [figdir,area,'_Ripplemodmod_NDays'];
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
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






