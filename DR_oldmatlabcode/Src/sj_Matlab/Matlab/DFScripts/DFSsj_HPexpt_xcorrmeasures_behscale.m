
% For getting behavioral time scale correlations and also sleep
% correlations

%xcorrmeasures2 - Only do cell-cell correlations. No need for place field
%overlap. Get rid of that part of the code. Also incorporate taking all
%cells regardless of how many epochs they are defined in. Combine later.


% Old Code Comments
% --------------
% Also see DFSsj_placefields1 and mkarlsso/xcorrscript
% First use to get and plot overlap, trajdata and mapdata for run epochs:
% Make a DFA_calcoverlap for this purpose
% Overlap: Dont Plot all at once for all days - too many plots!
% Then add calls to DFA_calcxcorrmeasures_linfields. Get xcorr measures for
% run and sleep epochs. Sleep filter can be speed


clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 1; % Figure Options -

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';

%savefile = [savedir 'HP_xcorrmeasures_behscale_CA1PFC']; area = ''; % CA1 vs PFC
savefile = [savedir 'HP_xcorrmeasures_behscale_PFC']; area = 'PFC'; % PFC vs PFC
%savefile = [savedir 'HP_xcorrmeasures_behscale_CA1']; area = 'CA1'; % CA1 vs CA1





% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days



% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    %     if iscon==0
    %         animals = {'REc,''REd','REe','REf'};
    %         %animals = {'REd','REe','REf'};
    %     else
    animals = {'HPa','HPb'};
    %animals = {'RCb','RCc','RE1'};
    %     end
    %animals = {'REe'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:4'; % Shantanu - I am adding day filter to parse out epoch filter
    runepochfilter = 'isequal($environment, ''wtr1'')'; % Only do 1st w-track. 2 or 1 epochs per day
    sleepepochfilter = 'isequal($type, ''sleep'')'; % Only pre and post sleep marked as sleep
    %presleepepochfilter = 'isequal($environment, ''presleep'')'; % Presleep. 1 epoch per day
    %postsleepepochfilter = 'isequal($environment, ''postsleep'')'; % Postsleep. 1 epoch per day
    
    % Cell filter
    % -----------
    
    %cellpairfilter = {'allcomb','strcmp($tag2, ''CA1Pyr'')','strcmp($tag2, ''PFC'')'}; % This includes all, including silent cells - might be only in sleep
    cellpairfilter = {'allcomb','strcmp($tag2, ''PFC'')','strcmp($tag2, ''PFC'')'};
    %cellpairfilter = {'allcomb','strcmp($tag2, ''CA1Pyr'')','strcmp($tag2, ''CA1Pyr'')'};
    % For more exlcusive choosing, use $tag instead of $tag2
    
    % Time filter
    % -----------
    
    % Could also use getripples for ripples across tetrode. Along with
    % cellthreshold condition as in getpopulation_events2, or sj_HPExpt_ripalign..
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Run time filter. Can also use gethighthetatimes2
    timefilter_place = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6},...
        {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} }; % Dont care about stim here. Get rid of all ripples incl artifacts
    
    % Get ripple times during run - with and without speed criterion. Do as you do in DFSsj_getriprate_run
    timefilterrun_onlyrip_speed = {{'DFTFsj_getvelpos', '(($absvel <= 5))'},...
        {'DFTFsj_getriptimes','($nripples > 0)','tetfilter',riptetfilter,'minthresh',3}};
    timefilterrun_onlyrip = {{'DFTFsj_getriptimes','($nripples > 0)','tetfilter',riptetfilter,'minthresh',3}};
    
    % Sleep Filter - Only Speed criterion
     timefiltersleep = {{'DFTFsj_getvelpos', '(($absvel <= 1))'}};
    
    % Get ripple times in sleep. Ripple on any tet of size > minthres. Stim times removed. No need for speed filter.
    % Use DFTFsj_get2dstate/DFTFsj_getvelpos for speed. Can also use immobility time defined in 2dstate. See also kk_get2dstate
    timefiltersleep_onlyrip = {{'DFTFsj_getriptimes','($nripples > 0)','tetfilter',riptetfilter,'minthresh',3},...
        {'DFTFsj_getvelpos', '(($absvel <= 1))'}};
    %OR use 2dstate instead of getvelpos
    %      timefiltersleep_onlyrip = {{'DFTFsj_getriptimes','($nripples > 0)','tetfilter',riptetfilter,'minthresh',3},...
    %         {'sj_get2dstate', '(($absvel <= 2))'},{'sj_get2dstate', '(($immobilitytime <= 5))'}};
    

    
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    
    % Sleep Corr All - Only speed criterion.
    sleepcorr = createfilter('animal',animals,'days',dayfilter,'epochs',sleepepochfilter,'cellpairs',...
        cellpairfilter,'excludetime', timefiltersleep,'iterator', iterator);
    
    % Sleep Corr - Timefilter: only during ripples
    sleepcorr_rip = createfilter('animal',animals,'days',dayfilter,'epochs',sleepepochfilter,'cellpairs',...
        cellpairfilter,'excludetime', timefiltersleep_onlyrip,'iterator', iterator);
    
    %     presleepcorr = createfilter('animal',animals,'days',dayfilter,'epochs',presleepepochfilter,'cellpairs',...
    %         cellpairfilter,'excludetime', timefiltersleep_onlyrip,'iterator', iterator);
    %     postsleepcorr = createfilter('animal',animals,'days',dayfilter,'epochs',postsleepepochfilter,'cellpairs',...
    %         cellpairfilter,'excludetime', timefiltersleep_onlyrip,'iterator', iterator);
    
    % Run Corr - Filter is speed critetion and no ripples - Beh Time-Scale Correlations
    runcorr = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
        cellpairfilter,'excludetime', timefilter_place,'iterator', iterator);
    
    % Run Corr - Correlations during theta. linvel and noripple criterion. 
    % Can also use gethighthetapower2 to get a more exclusive time filter than above
    runcorr_theta = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
        cellpairfilter,'excludetime', timefilter_place,'iterator', iterator);
    
    % Run Corr - Only ripples during run. Ripple criterion only.
    runcorr_rip = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
        cellpairfilter,'excludetime', timefilterrun_onlyrip,'iterator', iterator);
    
    %     % Run Corr - Only ripples during run at low speed. Ripple criterion and speed criterion.
    %     runcorr_ripsp = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
    %         cellpairfilter,'excludetime', timefilterrun_onlyrip_speed,'iterator', iterator);
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    
    %Can also use DFAsj_calpairxcorr
    
        % Sleep
        % ------
    sleepcorr = setfilterfunction(sleepcorr, 'DFAsj_HPexpt_calcxcorrmeasures', {'spikes'}); % 0.4s win by default
    sleepcorr_rip = setfilterfunction(sleepcorr_rip, 'DFAsj_HPexpt_calcxcorrmeasures', {'spikes'}); % 0.4s win by default
    %presleepcorr = setfilterfunction(presleepcorr, 'DFAsj_calcxcorrmeasures', {'spikes'});
    %postsleepcorr = setfilterfunction(postsleepcorr, 'DFAsj_calcxcorrmeasures', {'spikes'});
    
        % Run
        % ----
    runcorr = setfilterfunction(runcorr, 'DFAsj_HPexpt_calcxcorrmeasures', {'spikes'},'forbeh',1); % 5s win. Behavior time-scale correlations 
    runcorr_theta = setfilterfunction(runcorr_theta, 'DFAsj_HPexpt_calcxcorrmeasures', {'spikes'},'forripples',0); % 1s win
    runcorr_rip = setfilterfunction(runcorr_rip, 'DFAsj_HPexpt_calcxcorrmeasures', {'spikes'}); %0.4s win by default
    
    %runcorr_ripsp = setfilterfunction(runcorr_ripsp, 'DFAsj_calcxcorrmeasures', {'spikes'});
    
    % Run analysis
    % ------------
    sleepcorr = runfilter(sleepcorr);
    sleepcorr_rip = runfilter(sleepcorr_rip);
    %presleepcorr = runfilter(presleepcorr);
    %postsleepcorr = runfilter(postsleepcorr);
    runcorr = runfilter(runcorr);
    runcorr_theta = runfilter(runcorr_theta);
    runcorr_rip = runfilter(runcorr_rip);
    %runcorr_ripsp = runfilter(runcorr_ripsp);
    
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata
        save(savefile);
    end
    
else
    
    load(savefile);
    
end  % end runscript

if ~exist('savedata')
    return
end

% ----------------------------------
% Whether to gather data or to load previously gathered data
% --------------------------------------------------------------------
gatherdata = 0; savegatherdata = 0;
%gatherdatafile = [savedir 'HP_xcorrmeasures_CA1PFC_gather']; area =''; % CA1 vs PFC
gatherdatafile = [savedir 'HP_xcorrmeasures_PFC_gather']; area ='PFC'; % PFC vs PFC
%gatherdatafile = [savedir 'HP_xcorrmeasures_CA1_gather']; area = 'CA1';% CA1 vs CA1




%-----------------------------------------
% Gather Data
% ------------

% To Control Plotting, enter parameters here
if ~isempty(plotanimidx)
    useanim = plotanimidx;
else
    useanim = 1:length(sleepcorr); % Use all animals
end
if ~isempty(plotdays)
    usedays = plotdays;
else
    usedays = [];   % Get for each animal separately
end

if gatherdata
    
    % ---------------------------------
    % Parameters
    % -----------
    %corrwin=0.05; % 2*50ms=100ms for ripples (50 ms on each side). For run-theta, 100ms on each side
    % Increase corrwin for Hipp-PFC
    corrwin=0.1; % 2*100ms=200ms for ripples (50 ms on each side). For run-theta, 300ms on each side 2*300=600ms
    thresh_peakrate = 3; % Originally for place-field overlap. Could use this to parse CA1 cells as well.
    
    
    ripcorrtime=[]; runcorrtime=[];  thetacorrtime=[]; % Get timebase for ripcorr,runcorr, thetacorr only once
    allanimindex_sl=[]; normsmoothcorr_sl = [];
    allanimindex_run=[]; normsmoothcorr_runrip = []; normsmoothcorr_runtheta = []; normsmoothcorr_runall = [];
    cnt_sl=0; cnt_run=0; % This is not no of pairs
    
    % POSSIBLE CHANGE - INCLUDE PAIRS WHERE RIPPLE CORRLN IS 0, SINCE THAT IS A POSSIBLE RESULT.
    % LOW RUN CORRLN => NO CO_SPIKING DURING RIPPLES
    
    
    for an = 1:length(sleepcorr)
        
        % Sleep - Now sleepcorr is allcorr during sleep, sleepcorr_rip is Corrln during ripples
        for i=1:length(sleepcorr(an).output{1}),
            if ~isnan(sleepcorr(an).output{1}(i).normsmoothcorr) % or if Neventscorr~=0
                cnt_sl = cnt_sl+1;
                anim_index_sl{an}(cnt_sl,:)=sleepcorr(an).output{1}(i).index;
                %Only indexes
                animindex_sl=[an sleepcorr(an).output{1}(i).index]; % Put animal index in front
                allanimindex_sl = [allanimindex_sl; animindex_sl]; % Collect all Anim Day Epoch Tet1 Cell1 Tet2 Cell2 Index
                
                % Data
                normsmoothcorr_sl(cnt_sl,:) = sleepcorr(an).output{1}(i).normsmoothcorr;
                Neventscorr_sl(cnt_sl) = sleepcorr(an).output{1}(i).Neventscorr;
                coactivez_sl(cnt_sl) = sleepcorr(an).output{1}(i).coactivez;
                xcorr_sl{cnt_sl} = sleepcorr(an).output{1}(i).corr;
                
                normsmoothcorr_slrip(cnt_sl,:) = sleepcorr_rip(an).output{1}(i).normsmoothcorr;
                Neventscorr_slrip(cnt_sl) = sleepcorr_rip(an).output{1}(i).Neventscorr;
                coactivez_slrip(cnt_sl) = sleepcorr_rip(an).output{1}(i).coactivez;
                xcorr_slrip{cnt_sl} = sleepcorr_rip(an).output{1}(i).corr;
                
                % Time base
                if isempty(ripcorrtime)
                    if isfield(sleepcorr(an).output{1}(i).corr,'time');
                        ripcorrtime =  sleepcorr(an).output{1}(i).corr.time;
                    end
                end
            end
        end
        
        % Run - Corrln during ripples. Corrln during theta. And Corrln during all.
        for i=1:length(runcorr(an).output{1}),
            % None of the 3 run correlations must be undefined for this epoch.
            % Most likely is the ripple correlation
            %if ~isnan(runcorr_rip(an).output{1}(i).normsmoothcorr) & ~isnan(runcorr_theta(an).output{1}(i).normsmoothcorr) & ~isnan(runcorr(an).output{1}(i).normsmoothcorr)
            if (runcorr_rip(an).output{1}(i).Neventscorr > 0) && (runcorr_theta(an).output{1}(i).Neventscorr > 0) && (runcorr(an).output{1}(i).Neventscorr > 0)
                %if ~isnan(runcorr_rip(an).output{1}(i).normsmoothcorr)
                
                cnt_run = cnt_run+1;
                anim_index_run{an}(cnt_run,:)=runcorr(an).output{1}(i).index;
                %Only indexes
                animindex_run=[an runcorr(an).output{1}(i).index]; % Put animal index in front
                allanimindex_run = [allanimindex_run; animindex_run]; % Collect all Anim Day Epoch Tet1 Cell1 Tet2 Cell2 Index
                
                % Data for Run Ripples
                normsmoothcorr_runrip(cnt_run,:) = runcorr_rip(an).output{1}(i).normsmoothcorr;
                Neventscorr_runrip(cnt_run) = runcorr_rip(an).output{1}(i).Neventscorr;
                coactivez_runrip(cnt_run) = runcorr_rip(an).output{1}(i).coactivez;
                xcorr_runrip{cnt_run} = runcorr_rip(an).output{1}(i).corr;
                
                % Data for Run Theta
                normsmoothcorr_runtheta(cnt_run,:) = runcorr_theta(an).output{1}(i).normsmoothcorr;
                Neventscorr_runtheta(cnt_run) = runcorr_theta(an).output{1}(i).Neventscorr;
                coactivez_runtheta(cnt_run) = runcorr_theta(an).output{1}(i).coactivez;
                xcorr_runtheta{cnt_run} = runcorr_theta(an).output{1}(i).corr;
                
                % Data for Run All - Beh Time scale
                normsmoothcorr_runall(cnt_run,:) = runcorr(an).output{1}(i).normsmoothcorr;
                Neventscorr_runall(cnt_run) = runcorr(an).output{1}(i).Neventscorr;
                coactivez_runall(cnt_run) = runcorr(an).output{1}(i).coactivez;
                xcorr_runall{cnt_run} = runcorr(an).output{1}(i).corr;
                
                % Time base for theta run correlations, and beh timescale correlations
                if isempty(runcorrtime)
                    if isfield(runcorr(an).output{1}(i).corr,'time');
                        runcorrtime =  runcorr(an).output{1}(i).corr.time;
                    end
                end
                if isempty(thetacorrtime)
                    if isfield(runcorr_theta(an).output{1}(i).corr,'time');
                        thetacorrtime =  runcorr_theta(an).output{1}(i).corr.time;
                    end
                end
            end
        end
        
    end % end animal
    
    
    % Consolidate single cells across epochs. 2 methods: see DFSsj_getcellinfo
    % Here, I will use combined animal-index. Prevents having to loop over animals
    % ----------------------------------------------------------------------------
    
    % RUN DATA
    % ---------
    runpairoutput = struct;
    dummyindex=allanimindex_run;  % all anim-day-epoch-tet1-cell1-tet2-cell2 indices
    cntpairs=0;
    for i=1:size(allanimindex_run)
        animdaytetcell=allanimindex_run(i,[1 2 4 5 6 7]);
        ind=[];
        while rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7]))~=0          % collect all rows (epochs)
            ind = [ind rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7]))];        % finds the first matching row
            dummyindex(rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7])),:)=[0 0 0 0 0 0 0]; % after adding index, remove the corresponding row
            % so you could find the next one if it exists
        end
        
        % Gather everything for the current cell across epochs
        rip_normsmoothcorr=[]; rip_Nevents=[]; rip_coactivez=[];
        theta_normsmoothcorr=[]; theta_Nevents=[]; theta_coactivez=[];
        all_normsmoothcorr=[]; all_Nevents=[]; all_coactivez=[];
        for r=ind
            rip_normsmoothcorr = [rip_normsmoothcorr; normsmoothcorr_runrip(r,:)];
            rip_Nevents = [rip_Nevents; Neventscorr_runrip(r)];
            rip_coactivez = [rip_coactivez; coactivez_runrip(r)];
            theta_normsmoothcorr = [theta_normsmoothcorr; normsmoothcorr_runtheta(r,:)];
            theta_Nevents = [theta_Nevents; Neventscorr_runtheta(r)];
            theta_coactivez = [theta_coactivez; coactivez_runtheta(r)];
            all_normsmoothcorr = [all_normsmoothcorr; normsmoothcorr_runall(r,:)];
            all_Nevents = [all_Nevents; Neventscorr_runall(r)];
            all_coactivez = [all_coactivez; coactivez_runall(r)];
        end
        
        if ~isempty(rip_normsmoothcorr)
            cntpairs=cntpairs+1;
            runpairoutput_idx(cntpairs,:)=animdaytetcell;
            runpairoutput(cntpairs).index=animdaytetcell; % This is anim-day-tet1-cell1-tet2-cell2. No epoch
            runpairoutput(cntpairs).rip_corr = mean(rip_normsmoothcorr,1); % You can take mean across epochs later also. Then dont do it here.
            runpairoutput(cntpairs).rip_coactivez = mean(rip_coactivez);
            runpairoutput(cntpairs).rip_Nevents = sum(rip_Nevents); % Total no of events across epochs
            runpairoutput(cntpairs).theta_corr = mean(theta_normsmoothcorr,1);
            runpairoutput(cntpairs).theta_coactivez = mean(theta_coactivez);
            runpairoutput(cntpairs).theta_Nevents = sum(theta_Nevents);
            runpairoutput(cntpairs).all_corr = mean(all_normsmoothcorr,1); % You can take mean across epochs later also. Then dont do it here.
            runpairoutput(cntpairs).all_coactivez = mean(all_coactivez);
            runpairoutput(cntpairs).all_Nevents = sum(all_Nevents); % Total no of events across epochs
        end
    end
    
    
    
    % SLEEP DATA
    % -----------
    sleeppairoutput = struct;
    allanimindex_sleep = allanimindex_sl; dummyindex=allanimindex_sl;  % all anim-day-epoch-tet1-cell1-tet2-cell2 indices
    cnt_sleeppairs = 0; cnt_preslpairs=0; cnt_postslpairs=0;
    for i=1:size(allanimindex_sleep)
        animdaytetcell=allanimindex_sleep(i,[1 2 4 5 6 7]);
        ind=[]; ep_ind=[];
        while rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7]))~=0          % collect all rows (epochs)
            ind = [ind rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7]))];        % finds the first matching row
            % Get corresponding epoch index
            ep_ind = [ep_ind allanimindex_sleep(rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7])),3)];
            dummyindex(rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7])),:)=[0 0 0 0 0 0 0]; % after adding index, remove the corresponding row
            % so you could find the next one if it exists
        end
        
        % Gather everything for the current cell across epochs
        % Sleepcorr
        presl_normsmoothcorr=[]; presl_Nevents=[]; presl_coactivez=[];
        postsl_normsmoothcorr=[]; postsl_Nevents=[]; postsl_coactivez=[];
        % Sleepcorrrip
        preslrip_normsmoothcorr=[]; preslrip_Nevents=[]; preslrip_coactivez=[];
        postslrip_normsmoothcorr=[]; postslrip_Nevents=[]; postslrip_coactivez=[];
        
        % Can do pre and post sleep separately, or look for consensus for both epochs.
        % ----------------------------------------------------------------------------
        if length(ind)>1 % Consensus - Proceed only if both presleep and postsleep exist.
            for rcnt=1:length(ind)
                r = ind(rcnt);
                if ep_ind(rcnt)==1 % pre-sleep
                    presl_normsmoothcorr = [presl_normsmoothcorr; normsmoothcorr_sl(r,:)];
                    presl_Nevents = [presl_Nevents; Neventscorr_sl(r)];
                    presl_coactivez = [presl_coactivez; coactivez_sl(r)];
                    
                    preslrip_normsmoothcorr = [preslrip_normsmoothcorr; normsmoothcorr_slrip(r,:)];
                    preslrip_Nevents = [preslrip_Nevents; Neventscorr_slrip(r)];
                    preslrip_coactivez = [preslrip_coactivez; coactivez_slrip(r)];
                    
                else % post-sleep
                    postsl_normsmoothcorr = [postsl_normsmoothcorr; normsmoothcorr_sl(r,:)];
                    postsl_Nevents = [postsl_Nevents; Neventscorr_sl(r)];
                    postsl_coactivez = [postsl_coactivez; coactivez_sl(r)];
                    
                    postslrip_normsmoothcorr = [postslrip_normsmoothcorr; normsmoothcorr_slrip(r,:)];
                    postslrip_Nevents = [postslrip_Nevents; Neventscorr_slrip(r)];
                    postslrip_coactivez = [postslrip_coactivez; coactivez_slrip(r)];
                end
            end
        end
        
        if ~isempty(postsl_normsmoothcorr)
            cnt_sleeppairs=cnt_sleeppairs+1;
            sleeppairoutput_idx(cnt_sleeppairs,:)=animdaytetcell;
            sleeppairoutput(cnt_sleeppairs).index=animdaytetcell; % This is anim-day-tet1-cell1-tet2-cell2. No epoch
            
            % All Sleep
            sleeppairoutput(cnt_sleeppairs).presl_corr = mean(presl_normsmoothcorr,1); % You can take mean across epochs later also. Then dont do it here.
            sleeppairoutput(cnt_sleeppairs).presl_coactivez = mean(presl_coactivez);
            sleeppairoutput(cnt_sleeppairs).presl_Nevents = sum(presl_Nevents); % Total no of events across epochs
            sleeppairoutput(cnt_sleeppairs).postsl_corr = mean(postsl_normsmoothcorr,1);
            sleeppairoutput(cnt_sleeppairs).postsl_coactivez = mean(postsl_coactivez);
            sleeppairoutput(cnt_sleeppairs).postsl_Nevents = sum(postsl_Nevents);
            
            % Sleep Ripples
            sleeppairoutput(cnt_sleeppairs).preslrip_corr = mean(preslrip_normsmoothcorr,1); % You can take mean across epochs later also. Then dont do it here.
            sleeppairoutput(cnt_sleeppairs).preslrip_coactivez = mean(preslrip_coactivez);
            sleeppairoutput(cnt_sleeppairs).preslrip_Nevents = sum(preslrip_Nevents); % Total no of events across epochs
            sleeppairoutput(cnt_sleeppairs).postslrip_corr = mean(postslrip_normsmoothcorr,1);
            sleeppairoutput(cnt_sleeppairs).postslrip_coactivez = mean(postslrip_coactivez);
            sleeppairoutput(cnt_sleeppairs).postslrip_Nevents = sum(postslrip_Nevents);
        end
        
    end
    
    
    % Find corresponding pairs between Run and Sleep - Both have to Exist
    % --------------------------------------------------------------------
    
    cnt_runslpairs = 0; cnt_mismatch=0;
    for i=1:cntpairs % These are Run Pairs. Can also use while loop
        runidx = runpairoutput_idx(i,:);
        match = rowfind(runidx, sleeppairoutput_idx);
        if match~=0
            cnt_runslpairs = cnt_runslpairs+1;
            runsleeppair(cnt_runslpairs).index = runidx;
            
            % Sleep All corr
            runsleeppair(cnt_runslpairs).presl_corr = sleeppairoutput(match).presl_corr;
            runsleeppair(cnt_runslpairs).presl_Nevents = sleeppairoutput(match).presl_Nevents;
            runsleeppair(cnt_runslpairs).presl_coactivez = sleeppairoutput(match).presl_Nevents;
            runsleeppair(cnt_runslpairs).postsl_corr = sleeppairoutput(match).postsl_corr;
            runsleeppair(cnt_runslpairs).postsl_Nevents = sleeppairoutput(match).postsl_Nevents;
            runsleeppair(cnt_runslpairs).postsl_coactivez = sleeppairoutput(match).postsl_Nevents;
           
            % Sleep Ripple corr
            runsleeppair(cnt_runslpairs).preslrip_corr = sleeppairoutput(match).preslrip_corr;
            runsleeppair(cnt_runslpairs).preslrip_Nevents = sleeppairoutput(match).preslrip_Nevents;
            runsleeppair(cnt_runslpairs).preslrip_coactivez = sleeppairoutput(match).preslrip_Nevents;
            runsleeppair(cnt_runslpairs).postslrip_corr = sleeppairoutput(match).postslrip_corr;
            runsleeppair(cnt_runslpairs).postslrip_Nevents = sleeppairoutput(match).postslrip_Nevents;
            runsleeppair(cnt_runslpairs).postslrip_coactivez = sleeppairoutput(match).postslrip_Nevents;
            
            % Run corr
            runsleeppair(cnt_runslpairs).rip_corr = runpairoutput(i).rip_corr;
            runsleeppair(cnt_runslpairs).rip_coactivez = runpairoutput(i).rip_coactivez;
            runsleeppair(cnt_runslpairs).rip_Nevents = runpairoutput(i).rip_Nevents;
            runsleeppair(cnt_runslpairs).theta_corr = runpairoutput(i).theta_corr;
            runsleeppair(cnt_runslpairs).theta_coactivez = runpairoutput(i).theta_coactivez;
            runsleeppair(cnt_runslpairs).theta_Nevents = runpairoutput(i).theta_Nevents;
            runsleeppair(cnt_runslpairs).all_corr = runpairoutput(i).all_corr;
            runsleeppair(cnt_runslpairs).all_coactivez = runpairoutput(i).all_coactivez;
            runsleeppair(cnt_runslpairs).all_Nevents = runpairoutput(i).all_Nevents;
            
        else
            
            cnt_mismatch = cnt_mismatch+1; % How many run pairs were not found in sleep? - 182
            
        end
    end
    
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
    
end  % end gatherdata



% Calculations for Pure Run - Use a condition for N_events as well
% -----------------------------------------------------
bins_rip = find(abs(ripcorrtime)<=corrwin); % Corrln window is Between -100 and 100 ms for ripples
bins_run = find(abs(runcorrtime)<=10*corrwin); % Corrln window is between -1000 and 1000 ms
bins_theta = find(abs(thetacorrtime)<=4*corrwin); % Corrln window is between -400 and 400 ms
thrsev = 10;
cnt_runcorr=0;
for i=1:cntpairs
    rip_Nevents = runpairoutput(i).rip_Nevents; % rip_Nevents is the criterion. Do I need it?
    if rip_Nevents >= thrsev
        cnt_runcorr = cnt_runcorr+1;
        totalcorr_rip(cnt_runcorr) = sum(runpairoutput(i).rip_corr(bins_rip)); % Take total prob in 2*corrwin = -100 to +100ms around 0
        totalcorr_theta(cnt_runcorr) = sum(runpairoutput(i).theta_corr(bins_theta)); % Take total prob in 4*2*corrwin = -400 to +400ms around 0
        totalcorr_all(cnt_runcorr) = sum(runpairoutput(i).all_corr(bins_run)); % Take total prob in 10*2*corrwin = -1000 to +1000ms around 0
        coactiveZ_rip(cnt_runcorr) = runpairoutput(i).rip_coactivez;
        coactiveZ_theta(cnt_runcorr) = runpairoutput(i).theta_coactivez;
        coactiveZ_beh(cnt_runcorr) = runpairoutput(i).all_coactivez;
        % Get the new index as well
        totalcorr_index(cnt_runcorr,:) = runpairoutput(i).index;
        
    end
end




% Calculations for Run and Sleep Consensus - Use a condition for N_events as well
% -----------------------------------------------------

cnt_runslcorr=0;
for i=1:cnt_runslpairs
    rip_Nevents = runsleeppair(i).rip_Nevents; % rip_Nevents is the criterion. Do I need it?
    slrip_Nevents = runsleeppair(i).postsl_Nevents;
    
    if (rip_Nevents >= thrsev) && (slrip_Nevents >= thrsev)
        cnt_runslcorr = cnt_runslcorr+1;
        % Run
        totalrunsl_rip(cnt_runslcorr) = sum(runsleeppair(i).rip_corr(bins_rip)); % Take total prob in 2*corrwin = -100 to +100ms around 0
        totalrunsl_theta(cnt_runslcorr) = sum(runsleeppair(i).theta_corr(bins_theta)); % Take total prob in 4*2*corrwin = -400 to +400ms around 0
        totalrunsl_all(cnt_runslcorr) = sum(runsleeppair(i).all_corr(bins_run)); % Take total prob in 4*2*corrwin = -400 to +400ms around 0
        coactiveZ_rip2(cnt_runslcorr) = runsleeppair(i).rip_coactivez;
        coactiveZ_theta2(cnt_runslcorr) = runsleeppair(i).theta_coactivez;
        % Sleep
        totalrunsl_presl(cnt_runslcorr) = sum(runsleeppair(i).presl_corr(bins_rip));
        totalrunsl_postsl(cnt_runslcorr) = sum(runsleeppair(i).postsl_corr(bins_rip));
        totalrunsl_deltasl(cnt_runslcorr) = totalrunsl_postsl(cnt_runslcorr) - totalrunsl_presl(cnt_runslcorr);
        coactiveZ_presl(cnt_runslcorr) = runsleeppair(i).presl_coactivez;
        coactiveZ_postsl(cnt_runslcorr) = runsleeppair(i).postsl_coactivez;
        % Get the new index as well
        totalcorr_index(cnt_runslcorr,:) = runsleeppair(i).index;
    end
end


figopt1=1;
figdir = '/data25/sjadhav/HPExpt/Figures/Correlation/Egs/';

% Plotting for individual pairs
% ------------------------------
if figopt1==1,
    for i=1:cnt_runslpairs
        rip_Nevents = runsleeppair(i).rip_Nevents;
        theta_Nevents = runsleeppair(i).theta_Nevents;
        all_Nevents = runsleeppair(i).all_Nevents;
        idx = runsleeppair(i).index;
        
        figure; hold on;
        redimscreen_2versubplots	; % 1st is all run. 2nd will be ripple presl, run, and post-sl
        
        subplot(2,1,1); hold on;
        %plot(ripcorrtime, runsleeppair(i).rip_corr,'r','LineWidth',3);
        %plot(thetacorrtime, runsleeppair(i).theta_corr,'b','LineWidth',3);
        plot(runcorrtime, runsleeppair(i).all_corr,'k','LineWidth',3);
        line([0 0], [0 max(runsleeppair(i).all_corr)],'Color',[0.5 0.5 0.5],'LineWidth',2);
        
        set(gca,'XLim',[-4.8 4.8]);
        %set(gca,'XLim',[-0.4 0.4]);
        xlabel('Time (sec)','FontSize',20);
        ylabel('Norm corrln. - Run','FontSize',20);
        if idx(1)==1, pre = 'HPa'; else, pre = 'HPb'; end
        title(sprintf('%s Day%d Tet%d Cell%d, Tet%d Cell%d',...
            pre, idx(2), idx(3), idx(4), idx(5), idx(6)),'FontSize',20);
%         title(sprintf('%s Day%d Tet%d Cell%d Tet%d, Cell%d: Nevrip%d, Nevtheta%d, Nevall%d',...
%             pre, idx(2), idx(3), idx(4), idx(5), idx(6), rip_Nevents, theta_Nevents, all_Nevents ),'FontSize',20);
        
        subplot(2,1,2); hold on;
        %plot(ripcorrtime, runsleeppair(i).rip_corr,'r','LineWidth',3);
        plot(ripcorrtime, runsleeppair(i).postsl_corr,'k','LineWidth',3); % All Corr in post-sleep
        plot(ripcorrtime, runsleeppair(i).postslrip_corr,'b','LineWidth',3); % Ripple Corr in post-sleep
        %plot(ripcorrtime, runsleeppair(i).presl_corr,'k','LineWidth',3);
        line([0 0], [0 max([runsleeppair(i).rip_corr, runsleeppair(i).postsl_corr, runsleeppair(i).presl_corr])],'Color',[0.5 0.5 0.5],'LineWidth',2);
        set(gca,'XLim',[-0.35 0.35]);
        xlabel('Time (sec)','FontSize',20);
        ylabel('Norm corrln. - Sleep All and Sleep Rip','FontSize',20);
        %title(sprintf('Anim % d Day %d HpTet %d HpCell %d PFCtet %d, PFCcell %d Nevrip %d, Nevtheta %d, Nevall %d',...
        %    idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), rip_Nevents, theta_Nevents, all_Nevents ),'FontSize',tfont);
        
        figfile = [figdir,area,'EgCorr_',num2str(i)];
        keyboard;
        
        
        %print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
        
    end
    
end









%*********************************************
%   FIGURES
%********************************************

% ------------------------------
% Figure and Font Sizes

forppr = 0;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Figures/Correlation/';
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

clr = {'b','r','g','c','m','y','k','r'};
novel=1:3;
fam=4:8;
savefig1=0;




%****************************************************
% Overlap during run vs Correlation during run (Obvious: Confirmn) and Corrln during post-sleep
% ---------------------------------------------------

% all_ov = [overlap_ov; nooverlap_ov; midoverlap_ov];  % Overlap value
% allruncorr = [overlap_runallcorr; nooverlap_runallcorr; midoverlap_runallcorr]; % Corrln during all run
% allruncorrpl = [overlap_runplcorr; nooverlap_runplcorr; midoverlap_runplcorr]; % Corrln during place run
% allpostslcorr = [overlap_meancorr_post; nooverlap_meancorr_post; midoverlap_meancorr_post]; % Corrln during post-sleep
% allpreslcorr = [overlap_meancorr_pre; nooverlap_meancorr_pre; midoverlap_meancorr_pre]; % Corrln during pre-sleep
% alldiffslcorr = [allpostslcorr-allpreslcorr];



% ------------------------------------------------------------------------------------------------
% Theta-Corr vs Run-RipCorr: Scatter and Binned
% ------------------------------------------------------------------------------------------------

%Co-active Z Rip run and theta run
% Scatter plot with fit
% ---------------------
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end
[r_thetavsrip,p_thetavsrip] = corrcoef(coactiveZ_rip, coactiveZ_theta);
plot(coactiveZ_theta,coactiveZ_rip,'ro','Markersize',8,'LineWidth',2); % Run on x-axis
%set(gca,'XLim',[-0.04 1.15]);
%set(gca,'YLim',[-0.4 0.7]);
xlabel(['Theta Co-activeZ: ' num2str(4*corrwin*1000) ' ms win'],'FontSize',xfont,'Fontweight','normal');
ylabel(['Ripple Co-activeZ: ' num2str(corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');
text(5,7,['Npairs:' num2str(length(totalcorr_rip))],'FontSize',xfont,'Fontweight','normal');

% Regression
% -----------
[b00,bint00,r00,rint00,stats00] = regress(coactiveZ_rip', [ones(size(coactiveZ_theta')) coactiveZ_theta']);
xpts = 0:0.01:max(coactiveZ_theta);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip

% Do regression after shifting data to make intercept 0
% ------------------------------------------------------
coactiveZ_rip_0 = coactiveZ_rip-mean(coactiveZ_rip);
coactiveZ_theta_0 = coactiveZ_theta-mean(coactiveZ_theta);
[b0,bint0,r0,rint0,stats0] = regress(coactiveZ_rip_0',[ones(size(coactiveZ_theta_0')) coactiveZ_theta_0']);
bfit0 = b0(1)+b0(2)*xpts;

rval = roundn(r_thetavsrip(1,2),-2);
pval = roundn(p_thetavsrip(1,2),-4);
rsquare = roundn(stats0(1),-2);
preg = roundn(stats0(3),-4);


% Shuffling
% ---------
for n=1:1000
    idxs = randperm(length(coactiveZ_rip));
    shuffle = coactiveZ_rip(idxs);
    % Get corrcoeff of shuffle
    [rsh,psh] = corrcoef(coactiveZ_theta, shuffle);
    r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
    % Get regression of shuffle after making intercept 0 / Or Not
    %shuffle_0 = shuffle - mean(shuffle);
    %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(coactiveZ_theta_0')) coactiveZ_theta_0']);
    [bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle', [ones(size(coactiveZ_theta')) coactiveZ_theta']);
    rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
    b_shuffle(n,:) = bsh;
end
prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
% Get regression corresponding to 99 percentile
idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
plot(xpts,bfitsh,'k--','LineWidth',2);  % Theta vs Rip - 99% shuffle line

title(sprintf('Rip vs Theta: r=%g, p=%g; R^2=%g, preg=%g', rval, pval, rsquare, preg),...
    'FontSize',tfont,'Fontweight','normal');
if savefig1==1,
    figfile = [figdir,area,'RipvsThetaCoactiveZ_Scatter2'];
    print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end



% ------------------------------------------------------------------------------------------------
% Theta-Corr vs Run-RipCorr: Scatter and Binned
% ------------------------------------------------------------------------------------------------

% Find examples
egpair = find(totalcorr_rip>1 & totalcorr_theta>6);  %8 - HPa d1t1c2 t15c1
egpairs = find(totalcorr_rip>0.3 & totalcorr_theta>3); % Total 46: 2,6,8,16,29, 1331, 1434, 1460, 1472, 1535, 1545, 1554, 1601, 1604, etc


% Scatter plot with fit
% ---------------------
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end
[r_thetavsrip,p_thetavsrip] = corrcoef(totalcorr_theta, totalcorr_rip);
plot(totalcorr_theta, totalcorr_rip,'ro','Markersize',8,'LineWidth',2); % Run on x-axis
%set(gca,'XLim',[-0.04 1.15]);
%set(gca,'YLim',[-0.4 0.7]);
xlabel(['Theta Corrln: ' num2str(4*corrwin*1000) ' ms win'],'FontSize',xfont,'Fontweight','normal');
ylabel(['Ripple Corrln: ' num2str(corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');
text(4,1.4,['Npairs:' num2str(length(totalcorr_rip))],'FontSize',xfont,'Fontweight','normal');

% Regression
% -----------
[b00,bint00,r00,rint00,stats00] = regress(totalcorr_rip', [ones(size(totalcorr_theta')) totalcorr_theta']);
xpts = 0:0.01:max(totalcorr_theta);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip

% Do regression after shifting data to make intercept 0
% ------------------------------------------------------
totalcorr_rip_0 = totalcorr_rip-mean(totalcorr_rip);
totalcorr_theta_0 = totalcorr_theta-mean(totalcorr_theta);
[b0,bint0,r0,rint0,stats0] = regress(totalcorr_rip_0',[ones(size(totalcorr_theta_0')) totalcorr_theta_0']);
bfit0 = b0(1)+b0(2)*xpts;

rval = roundn(r_thetavsrip(1,2),-2);
pval = roundn(p_thetavsrip(1,2),-4);
rsquare = roundn(stats0(1),-2);
preg = roundn(stats0(3),-4);

% Shuffling
% ---------
for n=1:1000
    idxs = randperm(length(totalcorr_rip));
    shuffle = totalcorr_rip(idxs);
    % Get corrcoeff of shuffle
    [rsh,psh] = corrcoef(totalcorr_theta, shuffle);
    r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
    % Get regression of shuffle after making intercept 0 / Or Not
    %shuffle_0 = shuffle - mean(shuffle);
    %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(totalcorr_theta_0')) totalcorr_theta_0']);
    [bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle', [ones(size(totalcorr_theta')) totalcorr_theta']);
    rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
    b_shuffle(n,:) = bsh;
end
prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
% Get regression corresponding to 99 percentile
idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
plot(xpts,bfitsh,'k--','LineWidth',2);  % Theta vs Rip - 99% shuffle line

title(sprintf('Rip vs Theta: r=%g, p=%g; R^2=%g, preg=%g', rval, pval, rsquare, preg),...
    'FontSize',tfont,'Fontweight','normal');
if savefig1==1,
    figfile = [figdir,area,'RipvsThetaCorr_Scatter2'];
    print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end

% Range Plot
% ------------
% Divide in ranges instead of scatter plot
range = [0:0.5:4.5]; range = [range,7];
totalcorr_rip_range2=[]; totalcorr_rip_range2err=[];
plcorrrangex = range(1:end-1); % x-axis for plot
for i=1:length(range)-1
    rangest=range(i); rangeend=range(i+1);
    rangeidxs = find((totalcorr_theta>=rangest) & (totalcorr_theta<rangeend));
    totalcorr_rip_range2(i) = mean(totalcorr_rip(rangeidxs)); totalcorr_rip_range2err(i) = sem(totalcorr_rip(rangeidxs));
end
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end
errorbar(plcorrrangex,totalcorr_rip_range2,totalcorr_rip_range2err,'r','MarkerSize',16,'LineWidth',3);
set(gca,'XLim',[-0.1 max(plcorrrangex)+0.1]);
%set(gca,'YLim',[-0.04 0.6]);
xlabel(['Theta Corrln: ' num2str(4*corrwin*1000) ' ms win'],'FontSize',xfont,'Fontweight','normal');
ylabel(['Ripple Corrln: ' num2str(corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');
title(sprintf('Rip vs Theta Corrln'),'FontSize',tfont,'Fontweight','normal');
if savefig1==1,
    figfile = [figdir,area,'RipvsThetaCorr_Binned'];
    print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end



% ------------------------------------------------------------------------------------------------
% Theta-Corr vs PostSleep-RipCorr: Scatter and Binned
% ------------------------------------------------------------------------------------------------


% Scatter plot with fit
% ---------------------
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end
[r_thetavspostsl,p_thetavspostsl] = corrcoef(totalrunsl_theta, totalrunsl_postsl);
plot(totalrunsl_theta, totalrunsl_postsl,'ro','Markersize',8,'LineWidth',2); % Run on x-axis
%set(gca,'XLim',[-0.04 1.15]);
set(gca,'YLim',[-0.01 1.8]);
xlabel(['Theta Corrln: ' num2str(4*corrwin*1000) ' ms win'],'FontSize',xfont,'Fontweight','normal');
ylabel(['PostSleep Rip Corrln: ' num2str(corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');
text(4,1.4,['Npairs:' num2str(length(totalrunsl_postsl))],'FontSize',xfont,'Fontweight','normal');

% Regression
% -----------
[b00,bint00,r00,rint00,stats00] = regress(totalrunsl_postsl', [ones(size(totalrunsl_theta')) totalrunsl_theta']);
xpts = 0:0.01:max(totalrunsl_theta);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit0,'k-','LineWidth',4);  % Theta vs postsl

% Do regression after shifting data to make intercept 0
% ------------------------------------------------------
totalrunsl_postsl_0 = totalrunsl_postsl-mean(totalrunsl_postsl);
totalrunsl_theta_0 = totalrunsl_theta-mean(totalrunsl_theta);
[b0,bint0,r0,rint0,stats0] = regress(totalrunsl_postsl_0',[ones(size(totalrunsl_theta_0')) totalrunsl_theta_0']);
bfit0 = b0(1)+b0(2)*xpts;


rval = roundn(r_thetavspostsl(1,2),-2);
pval = roundn(p_thetavspostsl(1,2),-4);
rsquare = roundn(stats0(1),-2);
preg = roundn(stats0(3),-4);

% Shuffling
% ---------
for n=1:1000
    idxs = randperm(length(totalrunsl_postsl));
    shuffle = totalrunsl_postsl(idxs);
    % Get corrcoeff of shuffle
    [rsh,psh] = corrcoef(totalrunsl_theta, shuffle);
    r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
    % Get regression of shuffle after making intercept 0 / Or Not
    %shuffle_0 = shuffle - mean(shuffle);
    %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(totalrunsl_theta_0')) totalrunsl_theta_0']);
    [bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle', [ones(size(totalrunsl_theta')) totalrunsl_theta']);
    rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
    b_shuffle(n,:) = bsh;
end
prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
% Get regression corresponding to 99 percentile
idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
plot(xpts,bfitsh,'k--','LineWidth',2);  % Theta vs postsl - 99% shuffle line

title(sprintf('PostSlvsTheta: r=%g, p=%g; R^2=%g, preg=%g', rval, pval, rsquare, preg),...
    'FontSize',tfont,'Fontweight','normal');
if savefig1==1,
    figfile = [figdir,area,'PostSlRipvsThetaCorr_Scatter'];
    print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end

% Range Plot
% ------------
% Divide in ranges instead of scatter plot
range = [0:0.5:4.5]; range = [range,7];
totalrunsl_postsl_range2=[]; totalrunsl_postsl_range2err=[];
plcorrrangex = range(1:end-1); % x-axis for plot
for i=1:length(range)-1
    rangest=range(i); rangeend=range(i+1);
    rangeidxs = find((totalrunsl_theta>=rangest) & (totalrunsl_theta<rangeend));
    totalrunsl_postsl_range2(i) = mean(totalrunsl_postsl(rangeidxs)); totalrunsl_postsl_range2err(i) = sem(totalrunsl_postsl(rangeidxs));
end
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end
errorbar(plcorrrangex,totalrunsl_postsl_range2,totalrunsl_postsl_range2err,'r','MarkerSize',16,'LineWidth',3);
set(gca,'XLim',[-0.1 max(plcorrrangex)+0.1]);
%set(gca,'YLim',[-0.04 0.6]);
xlabel(['Theta Corrln: ' num2str(4*corrwin*1000) ' ms win'],'FontSize',xfont,'Fontweight','normal');
ylabel(['PostSl Rip Corrln: ' num2str(corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');
title(sprintf('PostSlRip vs Theta Corrln'),'FontSize',tfont,'Fontweight','normal');
if savefig1==1,
    figfile = [figdir,area,'PostSlRipvsThetaCorr_Binned'];
    print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end





% ------------------------------------------------------------------------------------------------
% Theta-Corr vs DeltaSleep-RipCorr: Scatter and Binned
% ------------------------------------------------------------------------------------------------


% Scatter plot with fit
% ---------------------
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end
[r_thetavsdeltasl,p_thetavsdeltasl] = corrcoef(totalrunsl_theta, totalrunsl_deltasl);
plot(totalrunsl_theta, totalrunsl_deltasl,'ro','Markersize',8,'LineWidth',2); % Run on x-axis
%set(gca,'XLim',[-0.04 1.15]);
%set(gca,'YLim',[-0.4 0.7]);
xlabel(['Theta Corrln: ' num2str(4*corrwin*1000) ' ms win'],'FontSize',xfont,'Fontweight','normal');
ylabel(['deltasleep Rip Corrln: ' num2str(corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');
text(4,1.4,['Npairs:' num2str(length(totalrunsl_deltasl))],'FontSize',xfont,'Fontweight','normal');

% Regression
% -----------
[b00,bint00,r00,rint00,stats00] = regress(totalrunsl_deltasl', [ones(size(totalrunsl_theta')) totalrunsl_theta']);
xpts = 0:0.01:max(totalrunsl_theta);
bfit00 = b00(1)+b00(2)*xpts;

% Do regression after shifting data to make intercept 0
% ------------------------------------------------------
totalrunsl_deltasl_0 = totalrunsl_deltasl-mean(totalrunsl_deltasl);
totalrunsl_theta_0 = totalrunsl_theta-mean(totalrunsl_theta);
[b0,bint0,r0,rint0,stats0] = regress(totalrunsl_deltasl_0',[ones(size(totalrunsl_theta_0')) totalrunsl_theta_0']);
bfit0 = b0(1)+b0(2)*xpts;
plot(xpts,bfit0,'k-','LineWidth',4);  % Theta vs deltasl

rval = roundn(r_thetavsdeltasl(1,2),-2);
pval = roundn(p_thetavsdeltasl(1,2),-4);
rsquare = roundn(stats0(1),-2);
preg = roundn(stats0(3),-4);

% Shuffling
% ---------
for n=1:1000
    idxs = randperm(length(totalrunsl_deltasl));
    shuffle = totalrunsl_deltasl(idxs);
    % Get corrcoeff of shuffle
    [rsh,psh] = corrcoef(totalrunsl_theta, shuffle);
    r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
    % Get regression of shuffle after making intercept 0 / Or Not
    shuffle_0 = shuffle - mean(shuffle);
    [bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(totalrunsl_theta_0')) totalrunsl_theta_0']);
    %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle', [ones(size(totalrunsl_theta')) totalrunsl_theta']);
    rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
    b_shuffle(n,:) = bsh;
end
prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
% Get regression corresponding to 99 percentile
idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
plot(xpts,bfitsh,'k--','LineWidth',2);  % Theta vs deltasl - 99% shuffle line

title(sprintf('DeltaSlvsTheta: r=%g, p=%g; R^2=%g, preg=%g', rval, pval, rsquare, preg),...
    'FontSize',tfont,'Fontweight','normal');
if savefig1==1,
    figfile = [figdir,area,'DeltaSlRipvsThetaCorr_Scatter'];
    print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end

% Range Plot
% ------------
% Divide in ranges instead of scatter plot
range = [0:0.5:3.5]; range = [range,7];
totalrunsl_deltasl_range2=[]; totalrunsl_deltasl_range2err=[];
plcorrrangex = range(1:end-1); % x-axis for plot
for i=1:length(range)-1
    rangest=range(i); rangeend=range(i+1);
    rangeidxs = find((totalrunsl_theta>=rangest) & (totalrunsl_theta<rangeend));
    totalrunsl_deltasl_range2(i) = mean(totalrunsl_deltasl(rangeidxs)); totalrunsl_deltasl_range2err(i) = sem(totalrunsl_deltasl(rangeidxs));
end
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end
errorbar(plcorrrangex,totalrunsl_deltasl_range2,totalrunsl_deltasl_range2err,'r','MarkerSize',16,'LineWidth',3);
set(gca,'XLim',[-0.1 max(plcorrrangex)+0.1]);
%set(gca,'YLim',[-0.04 0.6]);
xlabel(['Theta Corrln: ' num2str(4*corrwin*1000) ' ms win'],'FontSize',xfont,'Fontweight','normal');
ylabel(['DeltaSl Rip Corrln: ' num2str(corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');
title(sprintf('DeltaSlRip vs Theta Corrln'),'FontSize',tfont,'Fontweight','normal');
if savefig1==1,
    figfile = [figdir,area,'DeltaSlRipvsThetaCorr_Binned'];
    print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end















%
% %***************************************
% % Delta Corrln for Overlap vs No Overlap
% % --------------------------------------
% figure; hold on;
% if forppr==1
%     redimscreen_figforppr1;
% else
%     redimscreen_figforppt1;
% end
%
% meanov = mean(overlap_deltacorr);
% errov = sem(overlap_deltacorr);
% meanno = nanmean(nooverlap_deltacorr);
% errno = nansem(nooverlap_deltacorr);
%
% plot(1.1,meanno,'ro','MarkerSize',12,'LineWidth',2);
% plot(2.1,meanov,'ro','MarkerSize',12,'LineWidth',2);
% errorbar(1.1,meanno,errno,'r','LineWidth',2);
% errorbar(2.1,meanov,errov,'r','LineWidth',2);
% line([1.1 2.1], [meanno meanov],'Color','r','LineWidth',2);
%
% [hExp,pExp] = ttest2(overlap_deltacorr,nooverlap_deltacorr);
% if hExp==1,
%     mul = sign(mean(overlap_deltacorr));
%     plot(2.1, mean(overlap_deltacorr)+sem(overlap_deltacorr)+0.02, 'r*','MarkerSize',12);
% end
% title([grp,': Reactivn during rest'],'FontSize',tfont,'Fontweight','normal');
% ylabel(['Delta Correlation (Mean Prob in ' num2str(corrwin*1000) 'ms win)'],'FontSize',yfont,'Fontweight','normal');
% set(gca,'XTick',[1:2],'XTickLabel',{'No Overlap';'Overlap'},'FontSize',xfont,'Fontweight','normal');
% text(1.4,0.15,['p = ',num2str(roundn(pExp,-3))],'FontSize',tfont,'Color','r');






% ----------------------------------------------
i=1;
keyboard;

