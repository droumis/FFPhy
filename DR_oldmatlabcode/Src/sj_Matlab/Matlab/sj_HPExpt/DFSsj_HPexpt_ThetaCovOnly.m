
% Similar to DFSsj_HPexpt_ThetaCorrAndRipResp. 
% Only get ThetaCov and plot for PFC Theta phaselocked vs phase unlocked cells, or PFC ripmod vs PFC ripunmod cells
% Also try including only those Cov that exceed threshold, as in Siapas 2005.

clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 1; % Figure Options - Individual cells

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';

% CA1all files
% ------------
%savefile = [savedir 'HP_thetacov_CA1allPFCripmod']; area = 'CA1allPFCmod'; clr = 'r'; % CA1allPFC
%savefile = [savedir 'HP_thetacov_CA1allPFCripunmod']; area = 'CA1allPFCunmod'; clr = 'b'; % CA1allPFC
%savefile = [savedir 'HP_thetacov_CA1allPFCthetamod']; area = 'CA1allPFCthetamod'; clr = 'r'; % CA1allPFC
%savefile = [savedir 'HP_thetacov_CA1allPFCthetaunmod']; area = 'CA1allPFCthetaunmod'; clr = 'b'; % CA1allPFC

% Compare CA1 (theta modulated only) and PFC (theta modulated vs unmodulated)  
% -------------------------------------------------------------------------
%savefile = [savedir 'HP_thetacov_CA1thetamodPFCthetamod']; area = 'CA1thetamodPFCthetamod'; clr = 'r'; % 
%savefile = [savedir 'HP_thetacov_CA1thetamodPFCthetaunmod']; area = 'CA1thetamodPFCthetaunmod'; clr = 'b'; % 
% %Both unmodulated
%savefile = [savedir 'HP_thetacov_CA1thetaunmodPFCthetaunmod']; area = 'CA1thetaunmodPFCthetaunmod'; clr = 'g'; % PFC

% % Compare CA1 (theta modulated only) and PFC ripandthetabothmod vs unmod  
% ------------------------------------------------------------------------
%savefile = [savedir 'HP_thetacov_CA1thetaPFCboth']; area = 'CA1thetaPFCboth'; clr = 'm'; % 
%savefile = [savedir 'HP_thetacov_CA1thetaPFCunboth']; area = 'CA1thetaPFCunboth'; clr = 'c'; % 

% %IMP! CA1 (theta modulated only) and PFC ripmod vs ripunmod. Also Compare to Sleep Ripmod computed in other Script
% ------------------------------------------------------------------------
savefile = [savedir 'HP_thetacov_CA1thetamodPFCripmod']; area = 'CA1thetamodPFCripmod'; clr = 'r';
%savefile = [savedir 'HP_thetacov_CA1thetamodPFCripunmod']; area = 'CA1thetamodPFCripunmod'; clr = 'b';
% ---------

% IMP! For Sleep, get the propagted tag from the sleep to the run session. This way, only cells defined in sleep will
% be taken for the covariance computation. No need for doing a index match between run and sleep later.
% ---------
%savefile = [savedir 'HP_thetacov_CA1thetamodPFCsleepripmod']; area = 'CA1thetamodPFCsleepripmod'; clr = 'c'; %post-sleep only
%savefile = [savedir 'HP_thetacov_CA1thetamodPFCsleepripunmod']; area = 'CA1thetamodPFCsleepripunmod'; clr = 'g'; %post-sleep only

% % Theta and Ripple Mod mix
% --------------------------
% PFC only theta mod, but not ripple mod
% savefile = [savedir 'HP_thetacov_CA1thetamodPFConlythetamod']; area = 'CA1thetamodPFConlythetamod'; clr = 'k'; % 
% PFC only ripple mod, but not theta mod
%savefile = [savedir 'HP_thetacov_CA1thetamodPFConlyripmod']; area = 'CA1thetamodPFConlyripmod'; clr = 'g'; % 

% %IMP! CA1 ThetaMod vs PFC all - for sorting by ripmodln later. Run and Sleep
% ------------------------------------------------------------------------------
%savefile = [savedir 'HP_thetacov_CA1thetaPFCall']; area = 'CA1thetaPFCall'; clr = 'k'; % 
% All defined cells in postsleep: yes or no for postsleepripmodtag
%savefile = [savedir 'HP_thetacov_CA1thetamodPFCsleepall']; area = 'CA1thetamodPFCsleepall'; clr = 'k'; % 



savefig1=0;


% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPa','HPb','HPc'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    % Either Only do 1st w-track. 2 or 1 epochs per day
    % Or do Wtr1 and Wtr2, 2 epochs per day
    runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
    
    % %Cell filter
    % %-----------
    
     % CA1all calculation
    % -------------------
    % %CA1all-PFC theta modulated
    % %------------------------
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100)','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($thetamodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100)','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($thetamodtag, ''n'')'};
    
    % %CA1all-PFC ripple modulated vs unmod
    % %------------------------------------
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100)','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100)','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag, ''n'')'};
    

    % %CA1 theta mod - PFC both vs unboth
    % %---------------------------------
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($thetamodtag, ''y'') && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($thetamodtag, ''n'') && strcmp($ripmodtag, ''n'')'};
    
    
    % %CA1 theta mod  and PFC theta modulated vs unmodulated
    % %------------------------------------------------------------
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($thetamodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($thetamodtag, ''n'')'};
    % %Both unmodulated
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''n'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($thetamodtag, ''n'')'};
    
    
    % %IMP! CA1 (theta modulated only) and PFC ripmod vs ripunmod. Also Compare to Sleep Ripmod 
    % %---------------------------------------------------------
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag, ''n'')'};        
    
    % %IMP! CA1 (theta modulated only) and PFC Post-Sleep ripmod vs ripunmod.
    % %---------------------------------------------------------
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($postsleepripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($postsleepripmodtag, ''n'')'};

    % % Theta and Ripple Mod mix
    % --------------------------
    % PFC only theta mod, but not ripple mod
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($thetamodtag, ''y'') && strcmp($ripmodtag, ''n'')'};
    % PFC only ripple mod, but not theta mod
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($thetamodtag, ''n'') && strcmp($ripmodtag, ''y'')'};

    % %IMP! CA1 ThetaMod vs PFC all - for sorting by ripmodln later. Run and Sleep
    % ----------------------------------------------------------
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100)'};
    % All defined cells in postsleep: yes or no for postsleepripmodtag
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && (strcmp($postsleepripmodtag, ''y'') || strcmp($postsleepripmodtag, ''n'') )'};


    
   
    
    % Time filter - none.
    % -----------
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    timefilter_place = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6},...
        {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',2} };
    
    % Iterator
    % --------
    iterator = 'singlecellanal';  % Have defined cellpairfilter. Can also use cellpair iterator with cell defn
    
    % Filter creation
    % ----------------
    
    
    thetaf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter, 'cellpairs',...
        cellpairfilter, 'excludetime', timefilter_place, 'iterator', iterator);
    
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------

    
    %thetaf = setfilterfunction(thetaf,'DFAsj_calcxcorrmeasures', {'spikes'},'forripples',0);
    %thetaf = setfilterfunction(thetaf,'DFAsj_getthetacovariogram', {'spikes'}); % Shuffle corrected
    thetaf = setfilterfunction(thetaf,'DFAsj_getthetacrosscov', {'spikes'}); % Firing rate corected, as in Siapas, 2005, and Wierzynski 2009
    
    % Run analysis
    % ------------
    thetaf = runfilter(thetaf);
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata
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
gatherdata = 0; savegatherdata = 0;

% CA1all files
% ----
%gatherdatafile = [savedir 'HP_thetacov_CA1allPFCripmod_gather']; area = 'CA1allPFC'; kind = 'allvsripmod'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1allPFCripunmod_gather']; area = 'CA1allPFC'; kind = 'allvsripunmod'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1allPFCthetamod_gather']; area = 'CA1allPFC'; kind = 'allvsthetamod'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1allPFCthetaunmod_gather']; area = 'CA1allPFC'; kind = 'allvsthetaunmod'; state = '';

% Compare CA1 (theta modulated only) and PFC (theta modulated vs unmodulated)  
% -------------------------------------------------------------------------
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCthetamod_gather']; area = 'CA1thetamodPFCthetamod'; kind = 'thetamod'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCthetaunmod_gather']; area = 'CA1thetamodPFCthetaunmod'; kind = 'thetaunmod'; state = '';
% %Both unmodulated
%gatherdatafile = [savedir 'HP_thetacov_CA1thetaunmodPFCthetaunmod_gather']; area = 'CA1thetaunmodPFCthetaunmod'; kind = 'unmod'; state = '';


% % Compare CA1 (theta modulated only) and PFC ripandtheta both mod vs unmod  
% ------------------------------------------------------------------------
%gatherdatafile = [savedir 'HP_thetacov_CA1thetaPFCboth_gather']; area = 'CA1thetaPFCboth'; kind = 'thetavsboth'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1thetaPFCunboth_gather']; area = 'CA1thetaPFCunboth'; kind = 'thetavsunboth'; state = '';

% % IMP! - CA1 (theta modulated only) and PFC ripmod vs ripunmod. Also Compare to Sleep Ripmod computed in other Script
% ------------------------------------------------------------------------
gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCripmod_gather']; area = 'CA1thetamodPFCripmod'; kind = 'thetamodripmod'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCripunmod_gather']; area = 'CA1thetamodPFCripunmod'; kind = 'thetamodripunmod'; state = '';
% ------

% IMP! -Sleep RipMod
% ------
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCsleepripmod_gather']; area = 'CA1thetamodPFCsleepripmod'; kind = 'thetamodsleepripmod'; state = 'sleep';
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCsleepripunmod_gather']; area = 'CA1thetamodPFCsleepripunmod'; kind = 'thetamodsleepripunmod'; state = 'sleep';


% % Theta and Ripple Mod mix
% --------------------------
% PFC only theta mod, but not ripple mod
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFConlythetamod_gather']; area = 'CA1thetamodPFConlythetamod'; kind = 'thetamod'; state = '';
% PFC only ripple mod, but not theta mod
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFConlyripmod_gather']; area = 'CA1thetamodPFConlyripmod'; kind = 'thetamod'; state = '';

% % IMP! - CA1 ThetaMod vs PFC all - for soarting by ripmodln later. Run and Sleep
% -------------------------------------------------------------------------
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCall_gather']; area = 'CA1thetamodPFCall'; kind = 'thetamod'; state = '';
% All defined cells in postsleep: yes or no for postsleepripmodtag
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCsleepall_gather']; area = 'CA1thetamodPFCsleepall'; kind = 'thetamod'; state = 'sleep';


thrs = 3.66; % thrs = 2.58 (p=0.01) to thrs = 3.66 (p = 0.01/40)

if gatherdata
    
    % Parameters if any
    % -----------------
    
    % -------------------------------------------------------------
    
    cnt=0; 
    allanimindex=[];
    allZcrosscov_runtheta=[]; allcrosscov_runtheta_totalcorr=[]; allrawcorr_runtheta=[];
    allZcrosscov_sm_runtheta=[]; allcrosscov_sm_runtheta_totalcorr=[]; allrawcorr_sm_runtheta=[];
    allNeventscorr_runtheta=[];allxcorr_runtheta=[]; allT_runtheta=[]; allp1p2_runtheta=[];
    runcorrtime=[];
    corrwin = 0.2; %Window for theta corrln
    
    for an = 1:length(thetaf)
        for i=1:length(thetaf(an).output{1})
                cnt=cnt+1;
                anim_index{an}(cnt,:) = thetaf(an).output{1}(i).index;
                % Only indexes
                animindex=[an thetaf(an).output{1}(i).index]; % Put animal index in front
                allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
                
                % Data - Theta
                allZcrosscov_runtheta(cnt,:) = thetaf(an).output{1}(i).Zcrosscov;
                allcrosscov_runtheta(cnt,:) = thetaf(an).output{1}(i).crosscov;
                allrawcorr_runtheta(cnt,:) = thetaf(an).output{1}(i).rawcorr;
                allZcrosscov_sm_runtheta(cnt,:) = thetaf(an).output{1}(i).Zcrosscov_sm;
                allcrosscov_sm_runtheta(cnt,:) = thetaf(an).output{1}(i).crosscov_sm;
                allrawcorr_sm_runtheta(cnt,:) = thetaf(an).output{1}(i).rawcorr_sm;
                allNeventscorr_runtheta(cnt) = thetaf(an).output{1}(i).Neventscorr;
                allxcorr_runtheta{cnt} = thetaf(an).output{1}(i).corr;
                allT_runtheta(cnt) = thetaf(an).output{1}(i).T;
                allp1p2_runtheta(cnt) = thetaf(an).output{1}(i).p1p2;
                
                %Time base for theta correlations - only once
                if isempty(runcorrtime)
                    if isfield(thetaf(an).output{1}(i).corr,'time');
                        runcorrtime =  thetaf(an).output{1}(i).corr.time;
                    end                   
                    bins_run = find(abs(runcorrtime)<=corrwin); % +/- Corrln window                     
                end
                
                % Calculate a number for theta corr - Total prob or peak in -/+corrwin
                currthetacorr = allZcrosscov_runtheta(cnt,:);             
                currthetacorr_sm = allZcrosscov_sm_runtheta(cnt,:); 
                % Sum of values in window
                alltheta_totalcorr(cnt) = nansum(currthetacorr_sm(bins_run))./length(bins_run); % per bin
                % Peak value in window +/- corrwin
                alltheta_peakcorr(cnt) = nanmax(currthetacorr_sm(bins_run)); % Already smoothened, or can take +/-3 bins around peak
                if (~isnan(alltheta_peakcorr(cnt)) && ~isempty(alltheta_peakcorr(cnt)))
                    alltheta_peaklag_idx(cnt) = min(find(currthetacorr_sm(bins_run) == nanmax(currthetacorr_sm(bins_run)))); % in ms
                    alltheta_peaklag(cnt) = runcorrtime(bins_run(alltheta_peaklag_idx(cnt)))*1000; %in ms
                else
                    alltheta_peaklag_idx(cnt)=0;
                    alltheta_peaklag(cnt)=0;
                end
                % Trough value in window +/- corrwin
               
                alltheta_troughcorr(cnt) = nanmin(currthetacorr_sm(bins_run)); % Already smoothened, or can take +/-3 bins around peak
                if (~isnan(alltheta_troughcorr(cnt)) && ~isempty(alltheta_troughcorr(cnt)))
                    alltheta_troughlag_idx(cnt) = min(find(currthetacorr_sm(bins_run) == nanmin(currthetacorr_sm(bins_run)))); % in ms
                    alltheta_troughlag(cnt) = runcorrtime(bins_run(alltheta_troughlag_idx(cnt)))*1000; %in ms
                else
                    alltheta_troughlag_idx(cnt)=0;
                    alltheta_troughlag(cnt)=0;
                end
                     
                
                %alltheta_totalcorr(cnt) = nanmax(currthetacorr_sm(bins_run));
                %alltheta_peaklag(cnt) = find (currthetacorr(bins_run) == nanmax(currthetacorr(bins_run)));

        end
        
    end
    
   
   % NEED TO GET FOR EACH PFC cell, corresponding CA1-PFC pairs. Then get norm measure for that cell.
   
   Qdata = struct;
   cntc = 0; 
   dummyindex=allanimindex;  % all anim-day-epoch-tet1-cell1-tet2-cell2 indices
   for i=1:size(allanimindex)
        animdayeptet2cell2=allanimindex(i,[1 2 3 6 7]); % Unique PFC cell, including epoch       
        ind=[];
        while rowfind(animdayeptet2cell2,dummyindex(:,[1 2 3 6 7]))~=0          % collect all rows: matching CA1cells
            ind = [ind rowfind(animdayeptet2cell2,dummyindex(:,[1 2 3 6 7]))];        % finds the first matching row
            dummyindex(rowfind(animdayeptet2cell2,dummyindex(:,[1 2 3 6 7])),:)=[0 0 0 0 0 0 0]; % after adding index, remove the corresponding row
            % so you could find the next one if it exists
        end    
        
        % Gather the corresponding crosscov for current PFC cell
        currZsm = []; currZ = []; currpeak = []; currtrough = []; 
        for r=ind
            currZsm = [currZsm; allZcrosscov_sm_runtheta(r,:)]; 
            currZ = [currZ; allZcrosscov_runtheta(r,:)];
            currpeak = [currpeak; alltheta_peakcorr(r)]; 
            currtrough = [currtrough; alltheta_troughcorr(r)]; 
        end
        
        if ~isempty(ind)
            cntc = cntc+1;
            Qdata(cntc).index = animdayeptet2cell2;
            Qdata(cntc).Zsm = currZsm;
            Qdata(cntc).Z = currZ;
            Qdata(cntc).peak = currpeak;
            Qdata(cntc).trough = currtrough;
        end    
   end

   % Get mean cc for unique PFC cells
   Qmean = []; Qmeansig = []; Qmeansigonly = [];
   cntsig = 0; % How many PFC cells had at least 1 significant interaction?
   for i=1:cntc
       currZsm = Qdata(i).Zsm;
       currpeak = Qdata(i).peak; currtrough = Qdata(i).trough;
       
       npairs = size(currZsm,1);
       % All pairs
       Qmean(i,:) = nansum(currZsm,1)./(sqrt(npairs));
       
       % Only sig pairs
       %sigidx = find((currpeak >= thrs) | (currtrough <= -2*thrs));
       sigidx = find(currpeak >= thrs);
       if ~isempty(sigidx)
           cntsig = cntsig+1;
           nsigpairs = length(sigidx);         
           Qmeansig(i,:) = nansum(currZsm(sigidx,:),1)./(sqrt(nsigpairs));
           Qmeansigonly(cntsig,:) = nansum(currZsm(sigidx,:),1)./(sqrt(nsigpairs)); % Discarding cells with no sign CCs for check. Not correct way. 
       else
           Qmeansig(i,:) = zeros(size(runcorrtime));
       end
   end
    
   cntc
   cntsig
   
    
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
end % end gather data

figdir = '/data25/sjadhav/HPExpt/Figures/03Nov/';
set(0,'defaultaxesfontsize',20);
tfont = 20;
xfont = 20;
yfont = 20;
% Plotting for indiv pairs
% --------------------------
if 1
    for i=1:cnt      
        
        %if alltheta_peakcorr(i)>=2 || alltheta_troughcorr(i)<=-2
            
            idx = allanimindex(i,:);
            switch idx(1)
                case 1
                    pre ='HPa';
                case 2
                    pre = 'HPb';
                case 3
                    pre = 'HPc';
            end
            
            figure; hold on;
            plot(runcorrtime, allZcrosscov_sm_runtheta(i,:),'r','LineWidth',3);
            line([0 0], [min(allZcrosscov_sm_runtheta(i,:)) max(allZcrosscov_sm_runtheta(i,:))],'Color',[0.5 0.5 0.5],'LineWidth',2);
            line([0.2 0.2], [min(allZcrosscov_sm_runtheta(i,:)) max(allZcrosscov_sm_runtheta(i,:))],'Color',[0.5 0.5 0.5],'LineWidth',1);
            line([-0.2 -0.2], [min(allZcrosscov_sm_runtheta(i,:)) max(allZcrosscov_sm_runtheta(i,:))],'Color',[0.5 0.5 0.5],'LineWidth',1);
            
            title(sprintf('%s Day%d Ep %d Tet%d Cell%d, Tet%d Cell%d',...
                pre, idx(2), idx(3), idx(4), idx(5), idx(6), idx(7)),'FontSize',20);
            set(gca,'XLim',[-0.4 0.4]);
            xlabel('Time (sec)','FontSize',20);
            ylabel('Std. CrossCov - Run','FontSize',20);
            
            keyboard;
        %end
    end  
end


% Correct way is to get all unique PFC cells, and then get all PFC-CA1 pairs for that cell. 
% Then check for significance and get 
% mean stand. ccov = sum of sig. cov./(sqrt(no. of sig pairs), or
% mean stand. ccov = sum of cov./(sqrt(no. of pairs) 


% ------------------
% Population Figures
% ------------------

forppr = 0;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Figures/03Nov/';
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


if 1    
    
        
    % Plot Qmean sig 
    % -----------
    
    %Qmeansigpop = nanmean(Qmeansig,1);
    
    %Qmeansigpop = nanmean(Qmean,1);
    Qmeansigpop = nanmean(allZcrosscov_sm_runtheta,1);
    
    figure(3); hold on;
    plot(runcorrtime, Qmeansigpop,[clr '-'],'LineWidth',3);
    low = min(Qmeansigpop); high = max(Qmeansigpop);

    
    %title(sprintf('Mean Std. CrossCov - Theta %s Cells',kind),'FontSize',20);
    ylabel(sprintf('Mean Std. CrossCov'),'FontSize',24);
    set(gca,'XLim',[-0.4 0.4]);
    xlabel('Time (sec)','FontSize',24);
    title('Std. Mean CrossCov - Sig only','FontSize',24);
    
    low = -0.2; high = 1.3;
    line([0 0], [low high],'Color',[0.5 0.5 0.5],'LineWidth',2);
    line([0.2 0.2], [low high],'Color',[0.5 0.5 0.5],'LineWidth',1);
    line([-0.2 -0.2], [low high],'Color',[0.5 0.5 0.5],'LineWidth',1);
    
    %figfile = [figdir,area,'_',x,'_CrossCov'];
    figfile = [figdir,'Run_CrossCov_Summ2'];
    if savefig1==1,
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end

    
    
    % Plot Qmean sigonly 
    % -----------
    
%     Qmeansigpop = nanmean(Qmeansigonly,1);
%     
%     figure(4); hold on;
%     plot(runcorrtime, Qmeansigpop,[clr '--'],'LineWidth',3);
%     low = min(Qmeansigpop); high = max(Qmeansigpop);
% 
%     
%     %title(sprintf('Mean Std. CrossCov - Theta %s Cells',kind),'FontSize',20);
%     ylabel(sprintf('Mean Std. CrossCov'),'FontSize',24);
%     set(gca,'XLim',[-0.4 0.4]);
%     xlabel('Time (sec)','FontSize',24);
%     title('Std. Mean CrossCov - Sig only with discard','FontSize',24);
%     
%     line([0 0], [low high],'Color',[0.5 0.5 0.5],'LineWidth',2);
%     line([0.2 0.2], [low high],'Color',[0.5 0.5 0.5],'LineWidth',1);
%     line([-0.2 -0.2], [low high],'Color',[0.5 0.5 0.5],'LineWidth',1);
%     
%     figfile = [figdir,area,'_Both_CrossCov'];
%     if savefig1==1,
%         print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
%     end
%     
    





    
%     
%     % Plot Qmean 
%     % -----------
%     
%     Qmeanpop = nanmean(Qmean,1);
%     
%     figure; hold on;
%     plot(runcorrtime, Qmeanpop,'b','LineWidth',3);
%     low = min(Qmeanpop); high = max(Qmeanpop);
%     line([0 0], [low high],'Color',[0.5 0.5 0.5],'LineWidth',2);
%     line([0.2 0.2], [low high],'Color',[0.5 0.5 0.5],'LineWidth',1);
%     line([-0.2 -0.2], [low high],'Color',[0.5 0.5 0.5],'LineWidth',1);
%     
%     title(sprintf('Mean Std. CrossCov - Theta %s Cells',...
%         kind),'FontSize',20)
%     set(gca,'XLim',[-0.4 0.4]);
%     xlabel('Time (sec)','FontSize',20);
%     ylabel('Std. Mean CrossCov','FontSize',20);
%     

    
    
    % Plot Mean Standardized Cross-Cov for Entire Popln
    % -----------------------------------------------------
    
    %     runthetaZsm = nanmean(allZcrosscov_sm_runtheta,1);
    %     runthetaZ = nanmean(allZcrosscov_runtheta,1);
    %     %runthetaripZ = nanmean(allZcrosscov_sm_runthetarip,1);
    %     figure; hold on;
    %     plot(runcorrtime, runthetaZ,'r--','LineWidth',3);
    %     plot(runcorrtime, runthetaZsm,'r','LineWidth',3);
    %     line([0 0], [min(runthetaZ) max(runthetaZ)],'Color',[0.5 0.5 0.5],'LineWidth',2);
    %     line([0.2 0.2], [min(runthetaZ) max(runthetaZ)],'Color',[0.5 0.5 0.5],'LineWidth',1);
    %     line([-0.2 -0.2], [min(runthetaZ) max(runthetaZ)],'Color',[0.5 0.5 0.5],'LineWidth',1);
    %
    %     title(sprintf('Mean Std. CrossCov - Theta %s Cells',...
    %         kind),'FontSize',20)
    %     set(gca,'XLim',[-0.4 0.4]);
    %     xlabel('Time (sec)','FontSize',20);
    %     ylabel('Std. CrossCov','FontSize',20);
    %     legend('Theta','ThetaSm');
    %
    
    % Plot Mean Standardized Cross-Cov for Sig Popln
    % -----------------------------------------------------
    
%     sigcc = [];
%     for i=1:cnt        
%         %if alltheta_peakcorr(i)>=2 || alltheta_troughcorr(i)<=-2
%         if alltheta_peakcorr(i)>=2 
%             sigcc = [sigcc; allZcrosscov_sm_runtheta(i,:)];    
%         end 
%     end
%     
%     popsigcc = nanmean(sigcc,1);
%     normpopsigcc = popsigcc./(sqrt(size(sigcc,1)));
%     %figure; hold on;
%     plot(runcorrtime, normpopsigcc,'r','LineWidth',3);
%     line([0 0], [min(normpopsigcc) max(normpopsigcc)],'Color',[0.5 0.5 0.5],'LineWidth',2);
%     line([0.2 0.2], [min(normpopsigcc) max(normpopsigcc)],'Color',[0.5 0.5 0.5],'LineWidth',1);
%     line([-0.2 -0.2], [min(normpopsigcc) max(normpopsigcc)],'Color',[0.5 0.5 0.5],'LineWidth',1);
%     title(sprintf('Mean SIG Std. CrossCov - Theta %s Cells',kind),'FontSize',20)
%     set(gca,'XLim',[-0.4 0.4]);
%     xlabel('Time (sec)','FontSize',20);
%     ylabel('Std. CrossCov','FontSize',20);    
end












% ------------------------------------------------------------------
% COMBINING PLOTS ACROSS FILES
% ------------------------------------------------------------------


keyboard;

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
figdir = '/data25/sjadhav/HPExpt/Figures/RippleMod/Popln/';

% Define area
area = 'CA1allPFC'; state ='sleep'; % state = '';, or state = 'sleep';

if strcmp(state,'sleep'),
    statename = 'Sleep';
else
    statename = 'Run';
end


% Modulated Units
load([savedir 'HP_ripmod',state,'_corrandcoactz_',area,'_gather'])
% Corr Fig
figure(1); hold on; redimscreen_figforppt1;
set(gcf, 'Position',[205 136 723 446]);
%xaxis = min(allr):0.1:max(allr);
xaxis = -1:0.05:1;
h = histc(allr,xaxis); normh = h./max(h);
plot(xaxis,normh,clr,'Linewidth',3);

% % Coactive Z Fig
% figure(2); hold on; redimscreen_figforppt1;
% set(gcf, 'Position',[205 136 723 446]);
% xaxis = min(allcoactivez):0.5:max(allcoactivez);
% h = histc(allcoactivez,xaxis); normh = h./max(h);
% plot(xaxis,normh,clr,'Linewidth',3);

allrm = allr; allpm = allp; allcozm = allcoactivez;


% UnModulated Units
load([savedir 'HP_ripunmod',state,'_corrandcoactz_',area,'_gather'])
% Corr Fig
figure(1); hold on;
xaxis = -1:0.05:1;
h = histc(allr,xaxis); normh = h./max(h);
plot(xaxis,normh,[clr '--'],'Linewidth',3);

% % Coactive Z Fig
% figure(2); hold on; redimscreen_figforppt1;
% set(gcf, 'Position',[205 136 723 446]);
% xaxis = min(allcoactivez):0.5:max(allcoactivez);
% h = histc(allcoactivez,xaxis); normh = h./max(h);
% plot(xaxis,normh,[clr '--'],'Linewidth',3);

legend('Rip Mod','Rip Unmod');
title(sprintf('%s %s - units: Corr Coeff Hist', area, statename),'FontSize',tfont,'Fontweight','normal')
xlabel('Corr Coeff','FontSize',xfont,'Fontweight','normal');
ylabel('Fraction of cells','FontSize',yfont,'Fontweight','normal');

corrRateRipMod = nanmean(allpm < 0.05), corrRateRipUnMod = nanmean(allp < 0.05),
[r_modunmod p_modunmod] = ttest2(allpm<0.05,allp<0.05)

set(gca,'XLim',[-0.2 0.25]);
text(0.07,0.7,sprintf('Corr Rate Mod: %0.2f',corrRateRipMod),'FontSize',30,'Fontweight','normal');
text(0.07,0.6,sprintf('Corr Rate Unmod: %0.2f',corrRateRipUnMod),'FontSize',30,'Fontweight','normal');
text(0.07,0.5,sprintf('Diff Sig: %0.3f',p_modunmod),'FontSize',30,'Fontweight','normal');

figfile = [figdir,area,'_',statename,'_RippleModvsUnmod_CorrCoeffHist']
if savefig1==1,
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end











