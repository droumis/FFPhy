% Ver4 : Starting 10Feb2014 - Sync codes with everyone

% Based on "DFSsj_HPexpt_ThetacorrAndRipresp_ver4.m". 
% Get time-lag of corrln peak for theta covariance and ripple periods (+/-500ms)
% Can do for epochs separately, avg the crosscorrlns, and then get the peak, peaklag
% Alternative - in next version, combine data from two epochs, and then do the crosscorrln. corsscov



% --- Notes from original code ------
% Ver3: For ripple corrln, combine the ripplemod across epochs, and then do the corrln.
% So the filter function will simply return the correlation function

% Ver2, Dec 2013 - Implement combining data across epochs 
% Raw Data Files Stay Same - Only Gatherdata Files will change


% For ripple modulated and unmodulated cells, get corrcoeff and coactivez using the ripple response
% Instead of aligning to ripples again, use the saved output oin ripplemod files

% For consolidating across epochs:
% If I run corrcoeff and coactivez using DFAsj_getripresp_corrandcoactz, I will have to take mean when
% consolidating across epochs. Instead, I can run a new functions "DFAsj_gettrialResps", then combine these
% across epochs, and run corrcoeff and coactivez here in this function.


% clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
plotpairs = 0; % pairs
plotpop = 1; %population
cyclemaps = 1;
saveDRfigs = 0;

savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/';
%savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
[y, m, d] = datevec(date);

gatherdata = 0;
savegatherdata = 0;
% 
% mkdir('/mnt/data25/sjadhav/HPExpt/Figures_DR/XCorrgrams_PFCPFC_HPa_day2');
% figdir = '/mnt/data25/sjadhav/HPExpt/Figures_DR/XCorrgrams_PFCPFC_HPa_day2/'; %make these match
% set(0,'defaultaxesfontsize',20);
tfont = 18;
xfont = 16;
yfont = 16;


% PARAMETER
% ---------

% For Ripple Correlation, use either a window around ripples. Parameters define window around ripples
% OR just use ripple times and xcorr. Choose below
% pret and postt window around SWR to use for ripple correlation
rip_pret = 500/1000; rip_postt = 500/1000; % in secs


% Version 4 onwards
% ----------------
% IMP! CA1 (theta modulated only) and PFC ripmod vs ripunmod.
% -----------------------------------------------------------

% PUT MANUAL DATA
val=4;savefile = [savedir 'HPall_PFCPFC_Lags_ThetaAndRip_std3_speed4_ntet2_alldata_2-28-2014']; area = 'PFCPFC'; clr = 'g';
    
% IMP! CA1 (theta modulated only) 
% -------------------------------
%val=2;savefile = [savedir 'HP_CA1_Lags_ThetaAndRip_std3_speed4_ntet2_alldata_2-18-2014']; area = 'CA1thetamod'; clr = 'b';

% IMP! PFC ripmod and PFC ripmod; and PFCall .
% ------------------------------
%val=3;savefile = [savedir 'HP_PFCripmod_Lags_ThetaAndRip_std3_speed4_ntet2_alldata_2-18-2014']; area = 'PFCripmod'; clr = 'c';
%val=4;savefile = [savedir 'HP_PFCall_Lags_ThetaAndRip_std3_speed4_ntet2_alldata_2-18-2014']; area = 'PFCall'; clr = 'c';


savefig1=0;


% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPa'};
%    animals = {'HPa','HPb','HPc'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '2'; % Shantanu - I am adding day filter to parse out epoch filter
    % Either Only do 1st w-track. 2 or 1 epochs per day
    % Or do Wtr1 and Wtr2, 2 epochs per day
    runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'') || isequal($environment, ''ytr'')';
    
    % %Cell filter
    % %-----------
    
    % %IMP! CA1theta-PFCRipmodulated
    % %------------------------
    switch val
        case 1  
            cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'') ','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag2, ''y'')'};
    % %IMP! CA1theta
    % %------------------------
        case 2
            cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'') ','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'') '};
    % %IMP! PFCRipmodulated only
    % %------------------------
        case 3
            cellpairfilter = {'allcomb','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag2, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag2, ''y'')'};
    % %IMP! PFCall
    % %------------------------
        case 4
            cellpairfilter = {'allcomb','strcmp($area, ''PFC'') && ($numspikes > 100)','strcmp($area, ''PFC'') && ($numspikes > 100)'};
    end
    
    
    % Time filter - none.
    % -----------
    riptetfilter = '(isequal($descrip, ''riptet''))';
    timefilter_rip = {{'DFTFsj_getvelpos', '($absvel <= 4)'},{'DFTFsj_getriptimes4','($nripples > 1)','tetfilter',riptetfilter,'minstd',3}};
      
%     timefilter_place = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6},...
%         {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',2} };

    % Use absvel instead of linearvel
    timefilter_place_new = { {'DFTFsj_getvelpos', '($absvel >= 5)'},...
        {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    
%     timefilter_place_new = { {'DFTFsj_getvelpos', '($absvel >= 5)'},{'DFTFsj_getlinstate', '($state ~= -1)', 6},...
%         {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
%     
    % Iterator
    % --------
    iterator = 'singlecellanal';  % Have defined cellpairfilter. Can also use cellpair iterator with cell defn
    
    % Filter creation
    % ----------------
    
    % Ripple Correlation
    % -
    modf = createfilter('animal',animals, 'days', dayfilter, 'epochs',runepochfilter, 'cellpairs',...
        cellpairfilter, 'excludetime', timefilter_rip, 'iterator', iterator);
    
    % Theta Correlation
    % -
    
    %     thetaf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter, 'cellpairs',...
    %         cellpairfilter, 'excludetime', timefilter_place_new, 'iterator', iterator);
    thetaf = createfilter('animal',animals, 'days', dayfilter,'epochs',runepochfilter, 'cellpairs',...
        cellpairfilter, 'excludetime', timefilter_place_new, 'iterator', iterator);

    % Place field overlap
    % -
    % Parameters
    norm_overlap = 1;  % For place-field overlap
    thresh_peakrate = 3; % For place-field overlap
    xrun = createfilter('animal',animals,'days', dayfilter,'epochs',runepochfilter,'cellpairs',...
        cellpairfilter,'iterator', iterator);
    
    
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
   
    % For ripple corrln, CHOOSE
    % a) use window around ripple, OR
    % % For calculating ripplealigned resp from scratch, you will need spikes, ripples, tetinfo, and pos
     modf = setfilterfunction(modf,'DFAsj_getrip_period_xcorr4',{'spikes', 'ripples', 'tetinfo', 'pos'},'dospeed',1,'lowsp_thrs',4,'minrip',2,'minstd',3,'pret',rip_pret,'postt',rip_postt); 
    % b) Use xcorrmeasure with riptimes
   % modf = setfilterfunction(modf, 'DFAsj_HPexpt_calcxcorrmeasures4', {'spikes'});
    
    thetaf = setfilterfunction(thetaf,'DFAsj_getthetacrosscovLAG_timecondition4', {'spikes'},'thrstime',1, 'tmax', 4); % DR 4 seconds on each side to match euston
   
     
  %comgr  xrun = setfilterfunction(xrun, 'DFAsj_calcoverlap4', {'linfields', 'mapfields'},...
    %comgr    'normalize',norm_overlap,'thresh',thresh_peakrate,'minbins',0.5);
    
    % Run analysis
    % ------------
    modf = runfilter(modf);
    thetaf = runfilter(thetaf);
 
    %comgr   xrun = runfilter(xrun);
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata cyclemaps saveDRfigs gatherdata savegatherdata figdir tfont xfont yfont plotpairs plotpop
        save(savefile,'-v7.3');
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


% Version 4 onwards
% Filename of gatherdata. If generating, put current date. If not, then load a file from previous date.
% -----------------------------------------------------------------------------------------------------

% PUT MANUAL DATES


switch val
    % % IMP! - CA1 (theta modulated only) and PFC ripmod vs ripunmod. Also Compare to Sleep Ripmod computed in other Script
    % ------------------------------------------------------------------------
    case 1       
        gatherdatafile = [savedir 'HP_CA1PFC_Lags_ThetaAndRip_std3_speed4_ntet2_alldata_gather_2-18-2014'], area = 'CA1thetamodPFCripmod'; kind = 'thetamodripmod'; state = '';
    case 2
        gatherdatafile = [savedir 'HP_CA1_Lags_ThetaAndRip_std3_speed4_ntet2_alldata_gather_2-18-2014'], area = 'CA1thetamodPFCripUnmod'; kind = 'thetamodripUnmod'; state = '';
   
    case 3
        % % IMP! PFC Ripmod Only
        % -----------------------
        gatherdatafile = [savedir 'HP_PFCripmod_Lags_ThetaAndRip_std3_speed4_ntet2_alldata_gather_2-18-2014'], area = 'PFCripmod'; kind = 'ripmod'; state = '';
        
    case 4      
%         gatherdatafile = [savedir 'HP_PFCall_Lags_ThetaAndRip_std3_speed4_ntet2_alldata_gather_2-18-2014'], area = 'PFCripmod'; kind = 'ripmod'; state = '';
        gatherdatafile = [savedir 'HPa_PFCPFC_day2_Lags_ThetaAndRip_std3_speed4_ntet2_alldata_gather_2-28-2014'], area = 'PFCripmod'; kind = 'ripmod'; state = '';
        
end


if gatherdata
    
    % Parameters if any
    % -----------------
    
    % -------------------------------------------------------------
    
    cnt=0;
    allanimindex=[]; 
    
    allZcrosscov_runtheta=[]; allcrosscov_runtheta_totalcorr=[]; allrawcorr_runtheta=[]; allnormcorr_runtheta=[];
    allZcrosscov_sm_runtheta=[]; allcrosscov_sm_runtheta_totalcorr=[]; allrawcorr_sm_runtheta=[]; allnormcorr_sm_runtheta=[];
    allNeventscorr_runtheta=[];allxcorr_runtheta=[]; allT_runtheta=[]; allp1p2_runtheta=[];
    runcorrtime=[];
    
    allrawcorr_rip=[]; allrawcorr_sm_rip=[];  allnormcorr_rip=[]; allnormcorr_sm_rip=[];
    allNeventscorr_rip=[]; allxcorr_rip=[];
    ripcorrtime=[];
    
    alloverlap=[]; allpeakcomb=[]; alltrajpeakcomb=[];
    alltrajdata1=[]; alltrajdata2=[]; allmapdata1=[]; allmapdata2=[];
    
    corrwin = 0.5; %Window for theta and ripple corrln
    corrwinrip = 0.5;
    
    for an = 1:length(modf)
        for i=1:length(modf(an).output{1})
            cnt=cnt+1;
            anim_index{an}(cnt,:) = modf(an).output{1}(i).index;
            % Only indexes
            animindex=[an modf(an).output{1}(i).index]; % Put animal index in front
            allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet1 Cell1 Tet2 Cell2 Index

            % Data - Ripple period corrln
            % ----------------------------
            allrawcorr_rip(cnt,:) = modf(an).output{1}(i).rawcorr;
            allrawcorr_sm_rip(cnt,:) = modf(an).output{1}(i).rawcorr_sm;
            allnormcorr_rip(cnt,:) = modf(an).output{1}(i).normcorr;
            allnormcorr_sm_rip(cnt,:) = modf(an).output{1}(i).normcorr_sm;
            allNeventscorr_rip(cnt) = modf(an).output{1}(i).Neventscorr;
            allxcorr_rip{cnt} = modf(an).output{1}(i).corrstruct;
            
            % Data - Theta Crosscov
            % ---------------------
            allZcrosscov_runtheta(cnt,:) = thetaf(an).output{1}(i).Zcrosscov;
            allcrosscov_runtheta(cnt,:) = thetaf(an).output{1}(i).crosscov;
            allZcrosscov_sm_runtheta(cnt,:) = thetaf(an).output{1}(i).Zcrosscov_sm;
            allcrosscov_sm_runtheta(cnt,:) = thetaf(an).output{1}(i).crosscov_sm;
            allrawcorr_runtheta(cnt,:) = thetaf(an).output{1}(i).rawcorr;
            allrawcorr_sm_runtheta(cnt,:) = thetaf(an).output{1}(i).rawcorr_sm;
            allnormcorr_runtheta(cnt,:) = thetaf(an).output{1}(i).normcorr;
            allnormcorr_sm_runtheta(cnt,:) = thetaf(an).output{1}(i).normcorr_sm;
            allNeventscorr_runtheta(cnt) = thetaf(an).output{1}(i).Neventscorr;
            allxcorr_runtheta{cnt} = thetaf(an).output{1}(i).corrstruct;
            allT_runtheta(cnt) = thetaf(an).output{1}(i).T;
            allp1p2_runtheta(cnt) = thetaf(an).output{1}(i).p1p2;
            
            % Data - Place field overlap
            % ---------------------
%comgr            alloverlap(cnt) = xrun(an).output{1}(i).overlap;
      %comgr      allpeakcomb(cnt) = xrun(an).output{1}(i).peakcomb;
            %comgr alltrajpeakcomb(cnt) = xrun(an).output{1}(i).trajpeakcomb;
            % trajdata and mapdata are not being returend - too big
            %alltrajdata1{an}{i} = xrun(an).output{1}(i).trajdata1;
            %alltrajdata2{an}{i} = xrun(an).output{1}(i).trajdata2;
            %allmapdata1{an}{i} = xrun(an).output{1}(i).mapdata1;
            %allmapdata2{an}{i} = xrun(an).output{1}(i).mapdata2;
            
            
            %Time base for theta and ripple correlations - only once
            % -----------------------------------------------------
            if isempty(runcorrtime)
                if isfield(thetaf(an).output{1}(i).corrstruct,'time');
                    runcorrtime =  thetaf(an).output{1}(i).corrstruct.time;
                end
                bins_run = find(abs(runcorrtime)<=corrwin); % +/- Corrln window
            end

            if isempty(ripcorrtime)
                if isfield(modf(an).output{1}(i).corrstruct,'time');
                    ripcorrtime =  modf(an).output{1}(i).corrstruct.time;
                end
                bins_rip = find(abs(ripcorrtime)<=corrwinrip); % +/- Corrln window
            end
            
            
            % SHIFT ALL CALCULATIONS TO AFTER DATA COMBINED ACROSS EPOCHS
            % -------------------------------------------------------------           
%             % Calculate a number for theta corr - Total probor pek  in -/+corrwin
%             currthetacorr = allZcrosscov_runtheta(cnt,:);
%             currthetacorr_sm = allZcrosscov_sm_runtheta(cnt,:);
%             % Sum of values in window
%             alltheta_totalcorr(cnt) = nansum(currthetacorr_sm(bins_run))./length(bins_run); % per bin
%             % Peak value in window +/- corrwin
%             alltheta_peakcorr(cnt) = nanmax(currthetacorr_sm(bins_run)); % Already smoothened, or can take +/-3 bins around peak
%             if (~isnan(alltheta_peakcorr(cnt)) && ~isempty(alltheta_peakcorr(cnt)))
%                 alltheta_peaklag_idx(cnt) = min(find(currthetacorr_sm(bins_run) == nanmax(currthetacorr_sm(bins_run)))); % in ms
%                 alltheta_peaklag(cnt) = runcorrtime(bins_run(alltheta_peaklag_idx(cnt)))*1000; %in ms
%             else
%                 alltheta_peaklag_idx(cnt)=0;
%                 alltheta_peaklag(cnt)=0;
%             end
%             % Trough value in window +/- corrwin
%             
%             alltheta_troughcorr(cnt) = nanmin(currthetacorr_sm(bins_run)); % Already smoothened, or can take +/-3 bins around peak
%             if (~isnan(alltheta_troughcorr(cnt)) && ~isempty(alltheta_troughcorr(cnt)))
%                 alltheta_troughlag_idx(cnt) = min(find(currthetacorr_sm(bins_run) == nanmin(currthetacorr_sm(bins_run)))); % in ms
%                 alltheta_troughlag(cnt) = runcorrtime(bins_run(alltheta_troughlag_idx(cnt)))*1000; %in ms
%             else
%                 alltheta_troughlag_idx(cnt)=0;
%                 alltheta_troughlag(cnt)=0;
%             end
%           
%             %alltheta_totalcorr(cnt) = nanmax(currthetacorr(bins_run));
%             %alltheta_peaklag(cnt) = find (currthetacorr(bins_run) == nanmax(currthetacorr(bins_run)));
            
        end
        
    end
    
    
    % CONSOLIDATE ACROSS EPOCHS FOR SAME CELL PAIRS
    % ----------------------------------------------------------
    runpairoutput = struct;
    dummyindex=allanimindex;  % all anim-day-epoch-tet1-cell1-tet2-cell2 indices
    cntpairs=0;
    
    for i=1:size(allanimindex)
        animdaytetcell=allanimindex(i,[1 2 4 5 6 7]);
%         if sum(animdaytetcell(3:end)-[1 2 15 4])==0
%             
%             keyboard
%         end
        ind=[];
        while rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7]))~=0          % collect all rows (epochs)
            ind = [ind rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7]))];        % finds the first matching row
            dummyindex(rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7])),:)=[0 0 0 0 0 0 0]; % after adding index, remove the corresponding row
            % so you could find the next one if it exists
        end
        
        % Gather everything for the current cell across epochs
        % Theta corr variables
        % ---------------------
        allnormcorr_runtheta_epcomb = []; allnormcorr_sm_runtheta_epcomb = []; 
        allNeventscorr_runtheta_epcomb = []; %%alltheta_peakcorr_epcomb = []; alltheta_peaklag_epcomb = [];
        % Ripcorr variables 
        % -----------------
        allnormcorr_rip_epcomb = []; allnormcorr_sm_rip_epcomb = []; 
        allNeventscorr_rip_epcomb = []; 
        % Place Overlap
        % -----------------
        alloverlap_epcomb=[];
        
        for r=ind
            % Theta Corr variables
            % --------------------
            allnormcorr_sm_runtheta_epcomb = [allnormcorr_sm_runtheta_epcomb; allnormcorr_sm_runtheta(r,:)];
            allnormcorr_runtheta_epcomb = [allnormcorr_runtheta_epcomb; allnormcorr_runtheta(r,:)];
            allNeventscorr_runtheta_epcomb = [allNeventscorr_runtheta_epcomb; allNeventscorr_runtheta(r)];
            %%alltheta_peakcorr_epcomb = [alltheta_peakcorr_epcomb; alltheta_peakcorr(r)];
            %%alltheta_peaklag_epcomb = [alltheta_peaklag_epcomb; alltheta_peaklag(r)];
            % Rip Corr variables
            % ------------------
            allnormcorr_sm_rip_epcomb = [allnormcorr_sm_rip_epcomb; allnormcorr_sm_rip(r,:)];
            allnormcorr_rip_epcomb = [allnormcorr_rip_epcomb; allnormcorr_rip(r,:)];
            allNeventscorr_rip_epcomb = [allNeventscorr_rip_epcomb; allNeventscorr_rip(r)];
            % Place Overlap
            % -----------------
%comgr            alloverlap_epcomb = [alloverlap_epcomb; alloverlap(r)];
        end
        
        
        if ~isempty(allnormcorr_runtheta_epcomb)
            cntpairs=cntpairs+1;
            runpairoutput_idx(cntpairs,:)=animdaytetcell;
            runpairoutput(cntpairs).index=animdaytetcell; % This is anim-day-tet1-cell1-tet2-cell2. No epoch
            
            currthetacorr_sm=[]; currpeakcorr=[]; currpeaklag_idx=[]; currpeaklag =[];
            currnormcorr_sm=[]; currpeakripcorr=[]; currpeaklagrip_idx=[]; currpeaklagrip =[];
            
            % Theta corr variables
            % -------------------
            runpairoutput(cntpairs).allnormcorr_sm_runtheta_epcomb = nanmean(allnormcorr_sm_runtheta_epcomb,1);
            runpairoutput(cntpairs).allnormcorr_runtheta_epcomb = nanmean(allnormcorr_runtheta_epcomb,1);   
            runpairoutput(cntpairs).allNeventscorr_runtheta_epcomb = nansum([0;allNeventscorr_runtheta_epcomb]);
            % Get the peakcorr and peak lag from the combined normcorr_sm
            currthetacorr_sm = nanmean(allnormcorr_sm_runtheta_epcomb,1);
            currpeakcorr = nanmax(currthetacorr_sm(bins_run)); % Already smoothened, or can take +/-3 bins around peak
            if (~isnan(currpeakcorr) && ~isempty(currpeakcorr))
                currpeaklag_idx = min(find(currthetacorr_sm(bins_run) == nanmax(currthetacorr_sm(bins_run))));
                currpeaklag = runcorrtime(bins_run(currpeaklag_idx)); %in sec
            else
                currpeaklag_idx=NaN;
                currpeaklag=NaN;
            end         
            runpairoutput(cntpairs).alltheta_peakcorr_epcomb = currpeakcorr;
            runpairoutput(cntpairs).alltheta_peaklag_epcomb = currpeaklag;   
            %%runpairoutput(cntpairs).alltheta_peakcorr_epcomb = nanmean(alltheta_peakcorr_epcomb,1);
            %%runpairoutput(cntpairs).alltheta_peaklag_epcomb = nanmean(alltheta_peaklag_epcomb,1);                    
     
            
            % Rip Corr variables
            % -------------------
            runpairoutput(cntpairs).allnormcorr_rip_epcomb = nanmean(allnormcorr_rip_epcomb,1);  
            runpairoutput(cntpairs).allnormcorr_sm_rip_epcomb = nanmean(allnormcorr_sm_rip_epcomb,1);  
            runpairoutput(cntpairs).allNeventscorr_rip_epcomb = nansum([0;allNeventscorr_rip_epcomb]);
            % Get the peakcorr and peak lag from the combined normcorr_sm
            currnormcorr_sm = nanmean(allnormcorr_sm_rip_epcomb,1);  
            currpeakripcorr = nanmax(currnormcorr_sm(bins_rip)); % Already smoothened, or can take +/-3 bins around peak
            if (~isnan(currpeakripcorr) && ~isempty(currpeakripcorr))
                currpeaklagrip_idx = min(find(currnormcorr_sm(bins_rip) == nanmax(currnormcorr_sm(bins_rip))));
                currpeaklagrip = ripcorrtime(bins_rip(currpeaklagrip_idx)); %in sec
            else
                currpeaklagrip_idx = NaN;
                currpeaklagrip = NaN;
            end         
            runpairoutput(cntpairs).allrip_peakcorr_epcomb = currpeakripcorr;
            runpairoutput(cntpairs).allrip_peaklag_epcomb = currpeaklagrip;               
      
            % Place overlap 
            % --------------------
            runpairoutput(cntpairs).alloverlap_epcomb = nanmean(alloverlap_epcomb);
            
            
            % Save outside of structure format as well
            % Theta corr variables
            % ---------------------
            Sallnormcorr_sm_runtheta_epcomb(cntpairs,:) = nanmean(allnormcorr_sm_runtheta_epcomb,1);
            Sallnormcorr_runtheta_epcomb(cntpairs,:) = nanmean(allnormcorr_runtheta_epcomb,1);
            Salltheta_peakcorr_epcomb(cntpairs) = currpeakcorr;
            Salltheta_peaklag_epcomb(cntpairs) = currpeaklag;
            SallNeventscorr_runtheta_epcomb(cntpairs) = nansum([0;allNeventscorr_runtheta_epcomb]);
            % Rip Corr variables
            % -------------------
            Sallnormcorr_rip_epcomb(cntpairs,:) =  nanmean(allnormcorr_rip_epcomb,1);  
            Sallnormcorr_sm_rip_epcomb(cntpairs,:) =  nanmean(allnormcorr_sm_rip_epcomb,1);  
            Sallrip_peakcorr_epcomb(cntpairs) = currpeakripcorr;
            Sallrip_peaklag_epcomb(cntpairs) = currpeaklagrip;
            SallNeventscorr_rip_epcomb(cntpairs) = nansum([0;allNeventscorr_rip_epcomb]);
             % Place overlap 
            % --------------------
%comgr            Salloverlap_epcomb(cntpairs) = nanmean(alloverlap_epcomb);
        end
    end
    
%     figure; hold on; 
    %comgr hist(Salloverlap_epcomb);
    
    
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
end % end gather data

% Plotting for indiv pairs
% --------------------------
if plotpairs == 1;
    %DR load in the rip sig corr pair indices
%     saveDRfigs = 0;
% load('/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/HP_allPFCCA1sigcorridxs_feb14_DR');
% pairdata = allPFCCA1sigcorridxs;
% for yu = 1:length(pairdata);
%     PFClistindex(yu,:) = pairdata(1,yu).PFCidx(:);
% end

%DR load rip mod to get 'type'
% load('/mnt/data15/gideon/ProcessedData/HP_ripplemod_PFC_alldata_std3_speed4_ntet2_Nspk50_gather_2-19-2014');

    
%    for i=1:cnt       
%         if allp_shuf(i)<0.05   
sigcnt = 0
     for i=1:cntpairs
         % DR only use rip sig corr pairs
%         modind = rowfind(runpairoutput_idx(i,[1 2 5 6]),allripplemod_idx); %find the index of the pfc cell in the ripple mod list
%         if strcmp(allripplemod(modind).type, 'exc') && allripplemod(modind).sig_shuf > 0; %if the pfc cell is signicicantly rip excited 
%             if rowfind(runpairoutput_idx(i,[1 2 5 6]),PFClistindex);
%              pfcidx = rowfind(runpairoutput_idx(i,[1 2 5 6]),PFClistindex);
%              if rowfind(runpairoutput_idx(i,[1 2 3 4]),pairdata(1,pfcidx).CA1sigidxs); %now using only pos corr ca1 cells
                 sigcnt = sigcnt+1;
%                  if sigcnt ==22
%                      keyboard
%                  end
                 
        %if Sallp_epcomb(i)<0.05 
            idx = runpairoutput_idx(i,:) 
%comgr            curroverlap = Salloverlap_epcomb(i);
            %idx = allanimindex(i,:);
            switch idx(1)
                case 1
                    pre ='HPa';
                case 2
                    pre = 'HPb';
                case 3
                    pre = 'HPc';
            end
            
       %comgr     if ~isnan(curroverlap) && curroverlap>0.2 && SallNeventscorr_runtheta_epcomb(i)>=20 && SallNeventscorr_rip_epcomb(i)>=20
                  if SallNeventscorr_runtheta_epcomb(i)>=20 && SallNeventscorr_rip_epcomb(i)>=20
                
                figure; hold on; redimscreen_2versubplots;
                subplot(2,1,1); hold on;
                plot(runcorrtime, Sallnormcorr_sm_runtheta_epcomb(i,:),'r','LineWidth',3);
                plot(Salltheta_peaklag_epcomb(i), Salltheta_peakcorr_epcomb(i),'ko','MarkerSize',20,'LineWidth',2);
                
                
                line([0 0], [min(Sallnormcorr_sm_runtheta_epcomb(i,:)) max(Sallnormcorr_sm_runtheta_epcomb(i,:))],'Color',[0.5 0.5 0.5],'LineWidth',2);
                line([corrwin corrwin], [min(Sallnormcorr_sm_runtheta_epcomb(i,:)) max(Sallnormcorr_sm_runtheta_epcomb(i,:))],'Color',[0.5 0.5 0.5],'LineWidth',2);
                line([-corrwin -corrwin], [min(Sallnormcorr_sm_runtheta_epcomb(i,:)) max(Sallnormcorr_sm_runtheta_epcomb(i,:))],'Color',[0.5 0.5 0.5],'LineWidth',2);
                
               %comgr title(sprintf('%s Day%d Tet%d Cell%d, Tet%d Cell%d Overlap %g',...
               %comgr     pre, idx(2), idx(3), idx(4), idx(5), idx(6),roundn(curroverlap,-2)),'FontSize',20);
                title(sprintf('%s Day%d Tet%d Cell%d, Tet%d Cell%d ',...
                    pre, idx(2), idx(3), idx(4), idx(5), idx(6)),'FontSize',20);
             
               %if allp(i) <0.05, str = '*'; else, str = ''; end
                %text(0.2, 1*max(Sallnormcorr_sm_runtheta_epcomb(i,:)),sprintf('ripcc %0.2f%s',allr(i),str),'FontSize',20);
                set(gca,'XLim',[-4 4]);
                xlabel('Time (sec)','FontSize',20);
                ylabel('Xcorr - RunTheta','FontSize',20);
                
                subplot(2,1,2); hold on;
                plot(ripcorrtime, Sallnormcorr_sm_rip_epcomb(i,:),'b','LineWidth',3);
                plot(Sallrip_peaklag_epcomb(i), Sallrip_peakcorr_epcomb(i),'ko','MarkerSize',20,'LineWidth',2);
                
                corrwinrip = corrwin;
                
                line([0 0], [min(Sallnormcorr_sm_rip_epcomb(i,:)) max(Sallnormcorr_sm_rip_epcomb(i,:))],'Color',[0.5 0.5 0.5],'LineWidth',2);
                line([corrwinrip corrwinrip], [min(Sallnormcorr_sm_rip_epcomb(i,:)) max(Sallnormcorr_sm_rip_epcomb(i,:))],'Color',[0.5 0.5 0.5],'LineWidth',2);
                line([-corrwinrip -corrwinrip], [min(Sallnormcorr_sm_rip_epcomb(i,:)) max(Sallnormcorr_sm_rip_epcomb(i,:))],'Color',[0.5 0.5 0.5],'LineWidth',2);
                
                title(sprintf('Nevs-rip %d ; Nevs-theta %d', SallNeventscorr_rip_epcomb(i),SallNeventscorr_runtheta_epcomb(i)),'FontSize',20);
                
                set(gca,'XLim',[-.5 .5]);
                xlabel('Time (sec)','FontSize',20);
                ylabel('Xcorr - Rip','FontSize',20);
                
                
                %DR save the peaks of rip sig pairs
                ThetaSigPairsAllPeakLags(sigcnt,:) =  Salltheta_peaklag_epcomb(i);
                RipSigPairsAllPeakLags(sigcnt,:) =  Sallrip_peaklag_epcomb(i);
                                
                figfile = [figdir,'PFCPairsXCorr_thetaRip_',num2str(sigcnt)];
                if saveDRfigs==1
                    print('-dpng', figfile, '-r300');
                end
                if ~cyclemaps == 0
                    keyboard;
                end
                close
            end
        %end
%      end %if ca1 sig corr
% end %if pfc sig corr
%      end %if type matches
     end
end


if plotpop ==1;
     
     figure; 
     scatter(Salltheta_peaklag_epcomb, Sallrip_peaklag_epcomb,'MarkerEdgeColor','b','MarkerFaceColor','c','LineWidth',1.5);
     
     %Normalize
     
     for i = 1:length(Sallnormcorr_sm_rip_epcomb(:,1));
         Sallnormcorr_sm_rip_epcomb(i,~Sallnormcorr_sm_rip_epcomb(i,:)) = nan;         
         Sallnormcorr_sm_rip_epcomb(i,:) = (Sallnormcorr_sm_rip_epcomb(i,:)-nanmin(Sallnormcorr_sm_rip_epcomb(i,:)))./(nanmax(Sallnormcorr_sm_rip_epcomb(i,:))-nanmin(Sallnormcorr_sm_rip_epcomb(i,:)));
         Sallnormcorr_sm_rip_epcomb(i,isnan(Sallnormcorr_sm_rip_epcomb(i,:))) = 0;
     end
     
     for i = 1:length(Sallnormcorr_sm_runtheta_epcomb(:,1))   
         Sallnormcorr_sm_runtheta_epcomb(i,~Sallnormcorr_sm_runtheta_epcomb(i,:)) = nan;
         Sallnormcorr_sm_runtheta_epcomb_corrwin(i,bins_rip) = (Sallnormcorr_sm_runtheta_epcomb(i,bins_rip)-nanmin(Sallnormcorr_sm_runtheta_epcomb(i,bins_rip)))./(nanmax(Sallnormcorr_sm_runtheta_epcomb(i,bins_rip))-nanmin(Sallnormcorr_sm_runtheta_epcomb(i,bins_rip)));
         Sallnormcorr_sm_runtheta_epcomb_corrwin(i,isnan(Sallnormcorr_sm_runtheta_epcomb_corrwin(i,:))) = 0;
     end
     
     [B IX] =sort(Sallrip_peaklag_epcomb);
     Sallnormcorr_sm_rip_epcomb = Sallnormcorr_sm_rip_epcomb(IX,:);
     Sallnormcorr_sm_runtheta_epcomb_corrwin = Sallnormcorr_sm_runtheta_epcomb_corrwin(IX,:);
     
     figure;
     subplot(1,2,1);
     hold on
%      imagesc(runcorrtime, [1:length(Sallnormcorr_sm_runtheta_epcomb(:,1))], Sallnormcorr_sm_runtheta_epcomb(IX,:));
     imagesc(runcorrtime, [1:45], Sallnormcorr_sm_runtheta_epcomb_corrwin(IX,:));
%      imagesc(Sallnormcorr_sm_runtheta_epcomb);
set(gca,'YLim',[1 length(Sallnormcorr_sm_runtheta_epcomb_corrwin(:,1))]);
set(gca,'XLim',[-.5 .5]);
set(gca,'YDir','normal');

     subplot(1,2,2);
%      imagesc(ripcorrtime, [1:length(Sallnormcorr_sm_rip_epcomb(:,1))], Sallnormcorr_sm_rip_epcomb);
     imagesc(ripcorrtime, [1:45], Sallnormcorr_sm_rip_epcomb);
     set(gca,'YLim',[1 length(Sallnormcorr_sm_rip_epcomb(:,1))]);
     set(gca,'XLim',[-.5 .5]);
     set(gca,'YDir','normal');
%      imagesc(flipud(Sallnormcorr_sm_rip_epcomb));
     
     
%      rsthetarip = ranksum(Salltheta_peaklag_epcomb,Sallrip_peaklag_epcomb);
%      figure; hold on
%      bar([1], mean(Salltheta_peaklag_epcomb),'FaceColor', 'r', 'EdgeColor', 'none');
%      bar([2], mean(Sallrip_peaklag_epcomb), 'FaceColor', 'b', 'EdgeColor', 'none');
%      errorbar2([1 2], [mean(Salltheta_peaklag_epcomb) mean(Sallrip_peaklag_epcomb)],  [std(Salltheta_peaklag_epcomb) std(Sallrip_peaklag_epcomb)] , 0.3, 'k')
%      set(gca, 'xtick', [1 2], 'xticklabel', {'Theta', 'Rip'})
%      ylabel('Lag')
%      title({'HPa day2 XCorr Peak Lag PFCPairs'; sprintf('RS (%0.5f)',rsthetarip)});

     figfile = [figdir,'Lag_PFCPairsXCorr_thetaRip'];
     if saveDRfigs==1
         print('-dpng', figfile, '-r300');
     end
     if ~cyclemaps == 0
         keyboard;
     end
     close
     
end


return
keyboard

% ------------------
% Population Figures
% ------------------

forppr = 0;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Figures/Jan2014/';
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

if strcmp(state,'sleep'),
    statename = 'Sleep';
else
    statename = 'Run';
end

%clr = 'c'


% Change variable names for epoch-combined data so that you can use the same code to plot as in Version 1
% ----------------------------------------------------------------------------------------------------------
% Theta Corr variables
% allnormcorr_sm_runtheta = Sallnormcorr_sm_runtheta_epcomb;
% allnormcorr_runtheta = Sallnormcorr_runtheta_epcomb;
% alltheta_peakcorr = Salltheta_peakcorr_epcomb;
% alltheta_peaklag = Salltheta_peaklag_epcomb;







keyboard



if 1
    
    % 1a) Rip corrcoeff vs Total thetacorr for SWR Response
    % -----------------------------------------------------
    % Get rid of NaNs / Take Only Significant
    % Save temp
    alltheta_peaklag_tmp = Salltheta_peaklag_epcomb; % Take Peak Lag - thetacorr
    allrip_peaklag_tmp = Sallrip_peaklag_epcomb; % Take Peak Lag - ripcorr
    
   removeidxs = find((SallNeventscorr_rip_epcomb<20)|(SallNeventscorr_runtheta_epcomb<20) ) ;    
   %comgr    removeidxs = find((SallNeventscorr_rip_epcomb<20)|(SallNeventscorr_runtheta_epcomb<20) | (Salloverlap_epcomb<0.2)) ;    
    Sallrip_peaklag_epcomb(removeidxs)=[]; Salltheta_peaklag_epcomb(removeidxs)=[]; 
    
    figure; hold on; redimscreen_figforppt1;
    %set(gcf, 'Position',[205 136 723 446]);
    %xaxis = min(allr):0.1:max(allr);
    plot(Salltheta_peaklag_epcomb, Sallrip_peaklag_epcomb,'k.','MarkerSize',24);
    legend('Theta Xcorr Lag vs Rip Xcorr Lag');
    
    title(sprintf('Peak Lag Compare'),'FontSize',tfont,'Fontweight','normal')
    ylabel('Rip Xcorr Peak Lag ','FontSize',yfont,'Fontweight','normal');
    xlabel('Theta Xcorr Peak Lag','FontSize',xfont,'Fontweight','normal');
%     
%     % Do statistics on this popln plot
%     [r_thetavsrip,p_thetavsrip] = corrcoef(allr, alltheta_peakcorr);
%     [rsig,psig] = corrcoef(allr(sigidx), alltheta_peakcorr(sigidx));
%     [rnosig,pnosig] = corrcoef(allr(nosigidx), alltheta_peakcorr(nosigidx));
%     
%     % Regression
%     % -----------
%     [b00,bint00,r00,rint00,stats00] = regress(allr', [ones(size(alltheta_peakcorr')) alltheta_peakcorr']);
%     xpts = min(alltheta_peakcorr):0.01:max(alltheta_peakcorr);
%     bfit00 = b00(1)+b00(2)*xpts;
%     plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
%     % Do regression after shifting data to make intercept 0
%     % ------------------------------------------------------
%     allr_0 = allr-mean(allr);
%     alltheta_peakcorr_0 = alltheta_peakcorr-mean(alltheta_peakcorr);
%     [b0,bint0,r0,rint0,stats0] = regress(allr_0',[ones(size(alltheta_peakcorr_0')) alltheta_peakcorr_0']);
%     bfit0 = b0(1)+b0(2)*xpts;
%     
%     rval = roundn(r_thetavsrip(1,2),-2);
%     pval = roundn(p_thetavsrip(1,2),-4);
%     rsigval = roundn(rsig(1,2),-2);
%     psigval = roundn(psig(1,2),-4);
%     rnosigval = roundn(rnosig(1,2),-2);
%     pnosigval = roundn(pnosig(1,2),-4);
%     rsquare = roundn(stats0(1),-2);
%     preg = roundn(stats0(3),-4);
%     
%     % Shuffling
%     % ---------
%     for n=1:1000
%         idxs = randperm(length(allr));
%         shuffle = allr(idxs);
%         % Get corrcoeff of shuffle
%         [rsh,psh] = corrcoef(alltheta_peakcorr, shuffle);
%         r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
%         % Get regression of shuffle after making intercept 0 / Or Not
%         %shuffle_0 = shuffle - mean(shuffle);
%         %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(alltheta_peakcorr_0')) alltheta_peakcorr_0']);
%         [bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle', [ones(size(alltheta_peakcorr')) alltheta_peakcorr']);
%         rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
%         b_shuffle(n,:) = bsh;
%     end
%     prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
%     % Get regression corresponding to 99 percentile
%     idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
%     idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
%     bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
%     plot(xpts,bfitsh,'k--','LineWidth',2);  % Theta vs Rip - 99% shuffle line
%     
%     % Add Info to Figure
%     % ------------------
%     set(gca,'XLim',[-4.5 8]); set(gca,'YLim',[-0.2 0.45]);
%     text(-3.8,0.4,['Npairs:' num2str(length(allr))],'FontSize',30,'Fontweight','normal');
%     text(-3.8,0.35,sprintf('R: %0.2f, pval: %0.3f, preg: %0.3f',rval,pval,preg),'FontSize',30,'Fontweight','normal');
%     text(-3.8,0.3,sprintf('Rsig: %0.2f, pval: %0.3f, Nsig: %g',rsigval,psigval,length(sigidx)),'FontSize',30,'Fontweight','normal');
%     text(-3.8,0.25,sprintf('Rnosig: %0.2f, pval: %0.3f',rnosigval,pnosigval),'FontSize',30,'Fontweight','normal');        
%     figfile = [figdir,statename,'_ThetaCovVsRipCorr_',area]
%     if savefig1==1,
%         print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
%     end

    
    
end
    
    







    
    

