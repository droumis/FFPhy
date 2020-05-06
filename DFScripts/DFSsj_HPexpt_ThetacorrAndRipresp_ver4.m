% Ver5: CHange to averaging covariance acroos epochs, and then finding the peak
% Ver4 : Starting 10Feb2014 - Sync codes with everyone

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


clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - Individual cells

%savedir = '/data15/gideon/ProcessedData/';
savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
[y, m, d] = datevec(date);

% Version 4 onwards
% ----------------

if runscript==1
    % IMP! CA1 (theta modulated only) and PFC ripmod vs ripunmod.
    % -----------------------------------------------------------
    val=1;savefile = [savedir 'HP_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_',num2str(m),'-',num2str(d),'-',num2str(y)]; area = 'CA1thetamodPFCripmod'; clr = 'r';
    %val=2;savefile = [savedir 'HP_ThetacovAndRipUnmodcorr_std3_speed4_ntet2_alldata_',num2str(m),'-',num2str(d),'-',num2str(y)]; area = 'CA1thetamodPFCripUnmod'; clr = 'b';
    % IMP! PFC ripmod and PFC ripmod; and PFCall .
    % ------------------------------
    %val=3;savefile = [savedir 'HP_PFCripmod_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_',num2str(m),'-',num2str(d),'-',num2str(y)]; area = 'PFCripmod'; clr = 'c';
    %val=4;savefile = [savedir 'HP_PFCall_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_',num2str(m),'-',num2str(d),'-',num2str(y)]; area = 'PFCripmod'; clr = 'c';
else
    val=1;savefile = [savedir 'HP_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_2-12-2014']; area = 'CA1thetamodPFCripmod'; clr = 'r';
    %val=2; savefile = [savedir 'HP_ThetacovAndRipUnmodcorr_std3_speed4_ntet2_alldata_2-12-2014',]; area = 'CA1thetamodPFCripUnmod'; clr = 'b';
    %val=3; savefile = [savedir 'HP_PFCripmod_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_2-12-2014'];
    %val=4; savefile = [savedir 'HP_PFCall_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_2-12-2014'];
end



% % Pre version4
% ---------------
% ---------------
% IMP! CA1 (theta modulated only) and PFC ripmod vs ripunmod.
% -----------------------------------------------------------
%val=1;savefile = [savedir 'HP_ThetacovAndRipmodcorr_newtimefilter_alldata3_feb14']; area = 'CA1thetamodPFCripmod'; clr = 'r';
%val=2;savefile = [savedir 'HP_ThetacovAndRipUnmodcorr_newtimefilter_alldata3']; area = 'CA1thetamodPFCripUnmod'; clr = 'b';
% IMP! PFC ripmod and PFC ripmod .
% ------------------------------
%val=3;savefile = [savedir 'HP_PFC_ThetacovAndRipmodcorr_newtimefilter_alldata2']; area = 'PFCripmod'; clr = 'c';
% savefile = [savedir 'HP_thetaandripmodcorr_corrandcoactz_PFC']; area = 'PFC'; clr = 'b'; % PFC
% savefile = [savedir 'HP_thetaandripunmodcorr_corrandcoactz_PFC']; area = 'PFC'; clr = 'b'; % PFC
% savefile = [savedir 'HP_thetaandripmodcorr_corrandcoactz_CA1']; area = 'CA1'; clr = 'r';
% savefile = [savedir 'HP_thetaandripunmodcorr_corrandcoactz_CA1']; area = 'CA1'; clr = 'r';
% savefile = [savedir 'HP_thetaandripmodcorr_corrandcoactz_CA1PFC']; area = 'CA1PFC'; clr = 'c';
% savefile = [savedir 'HP_thetaandripunmodcorr_corrandcoactz_CA1PFC']; area = 'CA1PFC'; clr = 'c';
% savefile = [savedir 'HP_thetaandripmodcorr_corrandcoactz_CA1allPFC']; area = 'CA1allPFC'; clr = 'c';
% savefile = [savedir 'HP_thetaandripunmodcorr_corrandcoactz_CA1allPFC']; area = 'CA1allPFC'; clr = 'c';



savefig1=0;


% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    %animals = {'HPa','HPb','HPc','Nadal'};
    animals = {'HPa','HPb','HPc','Ndl'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    %dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    % Either Only do 1st w-track. 2 or 1 epochs per day
    % Or do Wtr1 and Wtr2, 2 epochs per day
    runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'') || isequal($environment, ''ytr'')';
    
    % %Cell filter
    % %-----------
    
    % %IMP! CA1theta-PFCRipmodulated
    % %------------------------
    switch val
        case 1  
            cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag2, ''y'')'};
        case 2
            cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($thetamodtag, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag2, ''n'')'};
    % %IMP! PFCRipmodulated only
    % %------------------------
        case 3
            cellpairfilter = {'allcomb','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag2, ''y'')','strcmp($area, ''PFC'') && ($numspikes > 100) && strcmp($ripmodtag2, ''y'')'};
    % %IMP! PFCall
    % %------------------------
        case 4
            cellpairfilter = {'allcomb','strcmp($area, ''PFC'') && ($numspikes > 100)','strcmp($area, ''PFC'') && ($numspikes > 100)'};
    end

    % % Pre version4
    % Misc - OLD
    % -----
    % %PFC
    % %----
    % cellpairfilter = {'allcomb','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')'}; % Ripple mod
    %cellpairfilter = {'allcomb','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')'}; % Ripple unmod
    % %CA1
    % %----
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($ripmodtag, ''y'')','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($ripmodtag, ''n'')','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($ripmodtag, ''n'')'};
    % %CA1-PFC
    % %--------
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($ripmodtag, ''y'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($numspikes > 100) && strcmp($ripmodtag, ''n'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')'};
    % %CA1all-PFCmodulated
    % %------------------------
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($meanrate > 0.2)','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) && ($meanrate > 0.2)','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')'};
    
    
    
    
    
    % Time filter - none.
    % -----------
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
%     timefilter_place = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6},...
%         {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',2} };


    % Version 4 onwards: 
    % ------------------------------------------------------------------------
    % Use absvel instead of linearvel
    timefilter_place_new = { {'DFTFsj_getvelpos', '(($absvel >= 5))'},...
        {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    
    % Iterator
    % --------
    iterator = 'singlecellanal';  % Have defined cellpairfilter. Can also use cellpair iterator with cell defn
    
    % Filter creation
    % ----------------
    
    % Ripp corrcoeff and coactivez
%     modf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter, 'cellpairs',...
%         cellpairfilter, 'iterator', iterator);
    modf = createfilter('animal',animals,'epochs',runepochfilter, 'cellpairs',...
        cellpairfilter, 'iterator', iterator);
    
%     thetaf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter, 'cellpairs',...
%         cellpairfilter, 'excludetime', timefilter_place, 'iterator', iterator);
    
%     thetaf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter, 'cellpairs',...
%         cellpairfilter, 'excludetime', timefilter_place_new, 'iterator', iterator);
     thetaf = createfilter('animal',animals,'epochs',runepochfilter, 'cellpairs',...
        cellpairfilter, 'excludetime', timefilter_place_new, 'iterator', iterator);
    
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    
    % % Need only the ripplemod structure for all of this.
    % % For calculating ripplealigned resp from scratch, you will need spikes, ripples, tetinfo, and pos
    modf = setfilterfunction(modf,'DFAsj_getripmoddata4',{'ripplemod'}); % ONLY RETURN RIPPLEMOD/TRIALRESPS DATA
    %modf = setfilterfunction(modf,'DFAsj_getripresp_corrandcoactz',{'ripplemod'}); %
    %modf = setfilterfunction(modf,'DFAsj_getripresp_corrandcoactz',{'ripplemod','cellinfo','spikes', 'ripples', 'tetinfo', 'pos'}); %
    
    %thetaf = setfilterfunction(thetaf,'DFAsj_calcxcorrmeasures', {'spikes'},'forripples',0);
    %thetaf = setfilterfunction(thetaf,'DFAsj_getthetacovariogram', {'spikes'});
    %thetaf = setfilterfunction(thetaf,'DFAsj_getthetacrosscov', {'spikes'});
    thetaf = setfilterfunction(thetaf,'DFAsj_getthetacrosscov_timecondition4', {'spikes'},'thrstime',1); % With includetime condition

    
    
    % Run analysis
    % ------------
    modf = runfilter(modf);
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
[y, m, d] = datevec(date);

% Version 4 onwards
% Filename of gatherdata. If generating, put current date. If not, then load a file from previous date.
% -----------------------------------------------------------------------------------------------------


% % IMP! - CA1 (theta modulated only) and PFC ripmod vs ripunmod. Also Compare to Sleep Ripmod computed in other Script
% ------------------------------------------------------------------------
switch val
    case 1  
        if gatherdata
            gatherdatafile = [savedir 'HP_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_gather_',num2str(m),'-',num2str(d),'-',num2str(y)], area = 'CA1thetamodPFCripmod'; kind = 'thetamodripmod'; state = '';
            %keyboard;
        else
            gatherdatafile = [savedir 'HP_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_gather_2-12-2014'], area = 'CA1thetamodPFCripmod'; kind = 'thetamodripmod'; state = '';
            %keyboard;
        end   
    case 2
        if gatherdata
            gatherdatafile = [savedir 'HP_ThetacovAndRipUnmodcorr_std3_speed4_ntet2_alldata_gather_',num2str(m),'-',num2str(d),'-',num2str(y)], area = 'CA1thetamodPFCripUnmod'; kind = 'thetamodripUnmod'; state = '';
        else
            gatherdatafile = [savedir 'HP_ThetacovAndRipUnmodcorr_std3_speed4_ntet2_alldata_gather_2-12-2014'], area = 'CA1thetamodPFCripUnmod'; kind = 'thetamodripUnmod'; state = '';
        end
% % IMP! PFC Ripmod Only
% -----------------------
    case 3
        if gatherdata
            gatherdatafile = [savedir 'HP_PFCripmod_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_gather_',num2str(m),'-',num2str(d),'-',num2str(y)], area = 'PFCripmod'; kind = 'ripmod'; state = ''; 
        else
            gatherdatafile = [savedir 'HP_PFCripmod_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_gather_2-12-2014'], area = 'PFCripmod'; kind = 'ripmod'; state = '';            
        end
    case 4
        if gatherdata
            gatherdatafile = [savedir 'HP_PFCall_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_gather_',num2str(m),'-',num2str(d),'-',num2str(y)], area = 'PFCripmod'; kind = 'ripmod'; state = ''; 
        else
            gatherdatafile = [savedir 'HP_PFCall_ThetacovAndRipmodcorr_std3_speed4_ntet2_alldata_gather_2-12-2014'], area = 'PFCripmod'; kind = 'ripmod'; state = '';            
        end
end




% Pre version4
% ------------
% IMP! - EPOCH COMBINED FILES
% -----------------------------------
% CA1 (theta modulated only) and PFC ripmod vs ripunmod.
% -----------------------------------------------------------
%gatherdatafile = [savedir 'HP_ThetacovAndRipmodcorr_newtimefilter_gather_epcomb']; area = 'CA1thetamodPFCripmod'; kind = 'thetamodripmod'; state = '';
% SAVED FILES FROM VER1 WITHOUT/WITH EPOCH COMBINE
% % IMP! - CA1 (theta modulated only) and PFC ripmod vs ripunmod. Also Compare to Sleep Ripmod computed in other Script
% ------------------------------------------------------------------------
% switch val
% case 1
% %gatherdatafile = [savedir 'HP_ThetacovAndRipmodcorr_newtimefilter_alldata_gather3']; area = 'CA1thetamodPFCripmod'; kind = 'thetamodripmod'; state = '';
% gatherdatafile = [savedir 'HP_ThetacovAndRipmodcorr_newtimefilter_alldata_gather3_feb14']; area = 'CA1thetamodPFCripmod'; kind = 'thetamodripmod'; state = '';
% case 2
% gatherdatafile = [savedir 'HP_ThetacovAndRipUnmodcorr_newtimefilter_alldata_gather3']; area = 'CA1thetamodPFCripUnmod'; kind = 'thetamodripUnmod'; state = '';
% % % IMP! PFC Ripmod Only
% % -----------------------
% case 3
% gatherdatafile = [savedir 'HP_PFC_ThetacovAndRipmodcorr_newtimefilter_alldata_gather3']; area = 'PFCripmod'; kind = 'ripmod'; state = '';
% end
% OLD
%gatherdatafile = [savedir 'HP_thetaandripmodcorr_corrandcoactz_PFC_gather']; area = 'PFC'; kind = 'mod'; state = '';
% gatherdatafile = [savedir 'HP_thetaandripunmodcorr_corrandcoactz_PFC_gather']; area = 'PFC'; kind = 'unmod'; state = '';
% gatherdatafile = [savedir 'HP_thetaandripmodcorr_corrandcoactz_CA1_gather']; area = 'CA1'; kind = 'mod'; state = '';
% gatherdatafile = [savedir 'HP_thetaandripunmodcorr_corrandcoactz_CA1_gather']; area = 'CA1'; kind = 'unmod'; state = '';
% gatherdatafile = [savedir 'HP_thetaandripmodcorr_corrandcoactz_CA1PFC_gather']; area = 'CA1PFC'; kind = 'mod'; state = '';
% gatherdatafile = [savedir 'HP_thetaandripunmodcorr_corrandcoactz_CA1PFC_gather']; area = 'CA1PFC'; kind = 'unmod'; state = '';
%gatherdatafile = [savedir 'HP_thetaandripmodcorr_corrandcoactz_CA1allPFC_gather']; area = 'CA1allPFC'; kind = 'mod'; state = '';
%gatherdatafile = [savedir 'HP_thetaandripunmodcorr_corrandcoactz_CA1allPFC_gather']; area = 'CA1allPFC'; kind = 'unmod'; state = '';


if gatherdata
    
    % Parameters if any
    % -----------------
    
    % -------------------------------------------------------------
    
    cnt=0;
    allanimindex=[]; 
    allZcrosscov_runtheta=[]; allcrosscov_runtheta_totalcorr=[]; allrawcorr_runtheta=[];
    allZcrosscov_sm_runtheta=[]; allcrosscov_sm_runtheta_totalcorr=[]; allrawcorr_sm_runtheta=[];
    allNeventscorr_runtheta=[];allxcorr_runtheta=[]; allT_runtheta=[]; allp1p2_runtheta=[];
    
    all_nsimul=[]; all_nsimul_rdm=[]; 
    all_trialResps1 = []; all_trialResps2 = []; all_trialResps1_rdm = []; all_trialResps2_rdm = [];
    runcorrtime=[];
    corrwin = 0.2; %Window for theta corrln
    
    for an = 1:length(modf)
        for i=1:length(modf(an).output{1})
            cnt=cnt+1;
            anim_index{an}(cnt,:) = modf(an).output{1}(i).index;
            % Only indexes
            animindex=[an modf(an).output{1}(i).index]; % Put animal index in front
            allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet1 Cell1 Tet2 Cell2 Index
            % Data - Ripmoddata
            all_trialResps1{cnt} = modf(an).output{1}(i).trialResps1;
            all_trialResps2{cnt} = modf(an).output{1}(i).trialResps2;
            all_trialResps1rdm{cnt} = modf(an).output{1}(i).trialResps1_rdm;
            all_trialResps2rdm{cnt} = modf(an).output{1}(i).trialResps2_rdm;
            all_nsimul(cnt) = modf(an).output{1}(i).nsimul;
            all_nsimulrdm(cnt) = modf(an).output{1}(i).nsimul_rdm;
            
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
        % --------------------
        allZcrosscov_runtheta_epcomb = []; allZcrosscov_sm_runtheta_epcomb = []; 
        allNeventscorr_runtheta_epcomb = []; %alltheta_peakcorr_epcomb = []; alltheta_peaklag_epcomb = [];
        % Ripcorr variables 
        % --------------------
        all_nsimul_epcomb=0; all_nsimulrdm_epcomb=0;
        all_trialResps1_epcomb=[]; all_trialResps2_epcomb=[]; all_trialResps1rdm_epcomb=[]; all_trialResps2rdm_epcomb=[];
        r=[]; p=[]; r_rdm=[]; p_rdm=[];
        
        for r=inde
            % Theta Corr variables
            % --------------------
            allZcrosscov_sm_runtheta_epcomb = [allZcrosscov_sm_runtheta_epcomb; allZcrosscov_sm_runtheta(r,:)];
            allZcrosscov_runtheta_epcomb = [allZcrosscov_runtheta_epcomb; allZcrosscov_runtheta(r,:)];
            allNeventscorr_runtheta_epcomb = [allNeventscorr_runtheta_epcomb; allNeventscorr_runtheta(r)];
            %alltheta_peakcorr_epcomb = [alltheta_peakcorr_epcomb; alltheta_peakcorr(r)];
            %alltheta_peaklag_epcomb = [alltheta_peaklag_epcomb; alltheta_peaklag(r)];
            % Rip Corr variables
            % --------------------
            all_trialResps1_epcomb = [all_trialResps1_epcomb; all_trialResps1{r}]; 
            all_trialResps2_epcomb = [all_trialResps2_epcomb; all_trialResps2{r}];
            all_trialResps1rdm_epcomb = [all_trialResps1rdm_epcomb; all_trialResps1rdm{r}]; 
            all_trialResps2rdm_epcomb = [all_trialResps2rdm_epcomb; all_trialResps2rdm{r}]; 
            all_nsimul_epcomb = all_nsimul_epcomb + all_nsimul(r);
            all_nsimulrdm_epcomb = all_nsimulrdm_epcomb + all_nsimulrdm(r);
        end
        
        % Calculate corrln for combined epoch data
        [r, p] = corrcoef(all_trialResps1_epcomb,all_trialResps2_epcomb);
        [r_rdm, p_rdm] = corrcoef(all_trialResps1rdm_epcomb,all_trialResps2rdm_epcomb);
        r = r(1,2); p = p(1,2);
        r_rdm = r_rdm(1,2); p_rdm = p_rdm(1,2);
        
        
        if ~isempty(allZcrosscov_runtheta_epcomb)
            cntpairs=cntpairs+1;
            runpairoutput_idx(cntpairs,:)=animdaytetcell;
            runpairoutput(cntpairs).index=animdaytetcell; % This is anim-day-tet1-cell1-tet2-cell2. No epoch
            
            currthetacorr_sm=[]; currpeakcorr=[]; currpeaklag_idx=[]; currpeaklag =[];            
            
            % Theta corr variables
            % --------------------
            runpairoutput(cntpairs).allZcrosscov_sm_runtheta_epcomb = nanmean(allZcrosscov_sm_runtheta_epcomb,1);
            runpairoutput(cntpairs).allZcrosscov_runtheta_epcomb = nanmean(allZcrosscov_runtheta_epcomb,1);
            runpairoutput(cntpairs).allNeventscorr_runtheta_epcomb = nansum([0;allNeventscorr_runtheta_epcomb]);
            % Get the peakcorr and peak lag from the combined normcorr_sm  
            currthetacorr_sm = nanmean(allZcrosscov_sm_runtheta_epcomb,1);
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
            %runpairoutput(cntpairs).alltheta_peakcorr_epcomb = nanmean(alltheta_peakcorr_epcomb,1);
            %runpairoutput(cntpairs).alltheta_peaklag_epcomb = nanmean(alltheta_peaklag_epcomb,1);
           
            % Rip Corr variables
            % --------------------
            runpairoutput(cntpairs).all_nsimul_epcomb = all_nsimul_epcomb;
            runpairoutput(cntpairs).all_nsimulrdm_epcomb = all_nsimulrdm_epcomb;
            runpairoutput(cntpairs).all_trialResps1_epcomb = all_trialResps1_epcomb;
            runpairoutput(cntpairs).all_trialResps2_epcomb = all_trialResps2_epcomb;
            runpairoutput(cntpairs).all_trialResps1rdm_epcomb = all_trialResps1rdm_epcomb;
            runpairoutput(cntpairs).all_trialResps2rdm_epcomb = all_trialResps2rdm_epcomb;
            runpairoutput(cntpairs).allr_epcomb = r;
            runpairoutput(cntpairs).allp_epcomb = p;
            runpairoutput(cntpairs).allr_rdm_epcomb = r_rdm;
            runpairoutput(cntpairs).allp_rdm_epcomb = p_rdm;
            
            % Save outside of structure format
            % Theta corr variables
            % --------------------
            SallZcrosscov_sm_runtheta_epcomb(cntpairs,:) = nanmean(allZcrosscov_sm_runtheta_epcomb,1);
            SallZcrosscov_runtheta_epcomb(cntpairs,:) = nanmean(allZcrosscov_runtheta_epcomb,1);
            Salltheta_peakcorr_epcomb(cntpairs) = currpeakcorr;
            Salltheta_peaklag_epcomb(cntpairs) = currpeaklag;
            SallNeventscorr_runtheta_epcomb(cntpairs) = nansum([0;allNeventscorr_runtheta_epcomb]);
            % Rip Corr variables
            % --------------------
            Sallr_epcomb(cntpairs) = r;
            Sallr_rdm_epcomb(cntpairs) = r_rdm;
            Sallp_epcomb(cntpairs) = p;
            Sallp_rdm_epcomb(cntpairs) = p_rdm;
            Sall_nsimul_epcomb(cntpairs) = all_nsimul_epcomb;
            Sall_nsimulrdm_epcomb(cntpairs) = all_nsimul_epcomb;
        end
    end
    
    
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
end % end gather data

figdir = '/data25/sjadhav/HPExpt/Figures/Correlation/Egs/';
set(0,'defaultaxesfontsize',20);
tfont = 20;
xfont = 20;
yfont = 20;
% Plotting for indiv pairs
% --------------------------
if 0
%    for i=1:cnt       
%         if allp_shuf(i)<0.05   
     for i=1:cntpairs       
        if Sallp_epcomb(i)<0.05 
            idx = runpairoutput_idx(i,:) 
            %idx = allanimindex(i,:);
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
            line([100 100], [min(allZcrosscov_sm_runtheta(i,:)) max(allZcrosscov_sm_runtheta(i,:))],'Color',[0.5 0.5 0.5],'LineWidth',2);
            line([-100 -100], [min(allZcrosscov_sm_runtheta(i,:)) max(allZcrosscov_sm_runtheta(i,:))],'Color',[0.5 0.5 0.5],'LineWidth',2);
            
            title(sprintf('%s Day%d Ep %d Tet%d Cell%d, Tet%d Cell%d',...
                pre, idx(2), idx(3), idx(4), idx(5), idx(6), idx(7)),'FontSize',20)
            if allp(i) <0.05, str = '*'; else, str = ''; end
            text(0.2, 1*max(allZcrosscov_sm_runtheta(i,:)),sprintf('ripcc %0.2f%s',allr(i),str),'FontSize',20);
            set(gca,'XLim',[-0.4 0.4]);
            xlabel('Time (sec)','FontSize',20);
            ylabel('Std. CrossCov - Run','FontSize',20);
            
            keyboard;
        end
    end
end



corrRateRipples = nanmean(Sallp_epcomb < 0.05)
corrRateRdm = nanmean(Sallp_rdm_epcomb < 0.05)
[r_realrdm, p_realrdm] = ttest(Sallp_epcomb<0.05,Sallp_rdm_epcomb<0.05)




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
allZcrosscov_sm_runtheta = SallZcrosscov_sm_runtheta_epcomb;
allZcrosscov_runtheta = SallZcrosscov_runtheta_epcomb;
alltheta_peakcorr = Salltheta_peakcorr_epcomb;
alltheta_peaklag = Salltheta_peaklag_epcomb;
% Rip Corr variables
allr = Sallr_epcomb;
allr_rdm = Sallr_rdm_epcomb;
allp = Sallp_epcomb;
allp_rdm = Sallp_rdm_epcomb;
all_nsimul = Sall_nsimul_epcomb;
all_nsimulrdm = Sall_nsimulrdm_epcomb;


% Plot Mean Standardized Cross-Cov for Entire Popln
% -----------------------------------------------------

if 0
    runthetaZsm = nanmean(allZcrosscov_sm_runtheta,1);
    runthetaZ = nanmean(allZcrosscov_runtheta,1);
    %runthetaripZ = nanmean(allZcrosscov_sm_runthetarip,1);
    figure; hold on;
    plot(runcorrtime, runthetaZ,'r--','LineWidth',3);
    plot(runcorrtime, runthetaZsm,'r','LineWidth',3);
    line([0 0], [min(runthetaZ) max(runthetaZ)],'Color',[0.5 0.5 0.5],'LineWidth',2);
    line([0.1 0.1], [min(runthetaZ) max(runthetaZ)],'Color',[0.5 0.5 0.5],'LineWidth',1);
    line([-0.1 -0.1], [min(runthetaZ) max(runthetaZ)],'Color',[0.5 0.5 0.5],'LineWidth',1);
    
    title(sprintf('Mean Cov - Rip %s Cells',...
        kind),'FontSize',20)
    set(gca,'XLim',[-0.4 0.4]);
    xlabel('Time (sec)','FontSize',20);
    ylabel('Std. CrossCov','FontSize',20);
    legend('Theta','ThetaSm');
    
    
    % Plot Mean Standardized Cross-Cov for Sig Ripcorr
    % -----------------------------------------------------
    
    nonsig = find (allp >= 0.05);
    tallZcrosscov_runtheta = allZcrosscov_runtheta; tallZcrosscov_sm_runtheta = allZcrosscov_sm_runtheta;
    tallZcrosscov_runtheta(nonsig,:)=[]; tallZcrosscov_sm_runtheta(nonsig,:)=[];
    runthetaZsm = nanmean(tallZcrosscov_sm_runtheta,1);
    runthetaZ = nanmean(tallZcrosscov_runtheta,1);
    %runthetaripZ = nanmean(allZcrosscov_sm_runthetarip,1);
    figure; hold on;
    plot(runcorrtime, runthetaZ,'r--','LineWidth',3);
    plot(runcorrtime, runthetaZsm,'r','LineWidth',3);
    line([0 0], [min(runthetaZ) max(runthetaZ)],'Color',[0.5 0.5 0.5],'LineWidth',2);
    line([0.1 0.1], [min(runthetaZ) max(runthetaZ)],'Color',[0.5 0.5 0.5],'LineWidth',1);
    line([-0.1 -0.1], [min(runthetaZ) max(runthetaZ)],'Color',[0.5 0.5 0.5],'LineWidth',1);
    
    title(sprintf('Mean Std. CrossCov - Rip %s Cells. Sig RipCorr Only',...
        kind),'FontSize',20)
    set(gca,'XLim',[-0.4 0.4]);
    xlabel('Time (sec)','FontSize',20);
    ylabel('Std. CrossCov','FontSize',20);
    legend('Theta','ThetaSm');
    
end





if 1
    
    % 1a) Rip corrcoeff vs Total thetacorr for SWR Response
    % -----------------------------------------------------
    % Get rid of NaNs / Take Only Significant
    % Save temp
    alltheta_peakcorr_tmp = alltheta_peakcorr; % Take Peak Corrln
    allr_tmp = allr; allp_tmp = allp;
    
    rnan = find(isnan(allr) | isnan(alltheta_peakcorr));
    %tnan = find(isnan(alltheta_peakcorr));
    nocc = find(all_nsimul<5);
    removeidxs = union(rnan,nocc);
    
    
    allr(removeidxs)=[]; alltheta_peakcorr(removeidxs)=[]; allp(removeidxs)=[];
    sigidx = find(allp<0.05); nosigidx = find(allp>=0.05);
    %allr = allr(sigidx); alltheta_peakcorr = alltheta_peakcorr(sigidx);
    
    figure; hold on; redimscreen_figforppt1;
    %set(gcf, 'Position',[205 136 723 446]);
    %xaxis = min(allr):0.1:max(allr);
    plot(alltheta_peakcorr, allr,'k.','MarkerSize',24);
    plot(alltheta_peakcorr(sigidx), allr(sigidx),'r.','MarkerSize',24);
    legend('Theta Cov vs SWR CorrCoeff');
    
    title(sprintf('%s %s units', area, statename),'FontSize',tfont,'Fontweight','normal')
    title(sprintf('CA1-PFC SWR modulated pairs'),'FontSize',tfont,'Fontweight','normal')
    ylabel('SWR Response Correlation','FontSize',yfont,'Fontweight','normal');
    xlabel('Peak Theta Covariance','FontSize',xfont,'Fontweight','normal');
    
    % Do statistics on this popln plot
    [r_thetavsrip,p_thetavsrip] = corrcoef(allr, alltheta_peakcorr);
    [rsig,psig] = corrcoef(allr(sigidx), alltheta_peakcorr(sigidx));
    [rnosig,pnosig] = corrcoef(allr(nosigidx), alltheta_peakcorr(nosigidx));
    
    % Regression
    % -----------
    [b00,bint00,r00,rint00,stats00] = regress(allr', [ones(size(alltheta_peakcorr')) alltheta_peakcorr']);
    xpts = min(alltheta_peakcorr):0.01:max(alltheta_peakcorr);
    bfit00 = b00(1)+b00(2)*xpts;
    plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
    % Do regression after shifting data to make intercept 0
    % ------------------------------------------------------
    allr_0 = allr-mean(allr);
    alltheta_peakcorr_0 = alltheta_peakcorr-mean(alltheta_peakcorr);
    [b0,bint0,r0,rint0,stats0] = regress(allr_0',[ones(size(alltheta_peakcorr_0')) alltheta_peakcorr_0']);
    bfit0 = b0(1)+b0(2)*xpts;
    
    rval = roundn(r_thetavsrip(1,2),-2);
    pval = roundn(p_thetavsrip(1,2),-4);
    rsigval = roundn(rsig(1,2),-2);
    psigval = roundn(psig(1,2),-4);
    rnosigval = roundn(rnosig(1,2),-2);
    pnosigval = roundn(pnosig(1,2),-4);
    rsquare = roundn(stats0(1),-2);
    preg = roundn(stats0(3),-4);
    
    % Shuffling
    % ---------
    for n=1:1000
        idxs = randperm(length(allr));
        shuffle = allr(idxs);
        % Get corrcoeff of shuffle
        [rsh,psh] = corrcoef(alltheta_peakcorr, shuffle);
        r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
        % Get regression of shuffle after making intercept 0 / Or Not
        %shuffle_0 = shuffle - mean(shuffle);
        %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(alltheta_peakcorr_0')) alltheta_peakcorr_0']);
        [bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle', [ones(size(alltheta_peakcorr')) alltheta_peakcorr']);
        rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
        b_shuffle(n,:) = bsh;
    end
    prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
    % Get regression corresponding to 99 percentile
    idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
    idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
    bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
    plot(xpts,bfitsh,'k--','LineWidth',2);  % Theta vs Rip - 99% shuffle line
    
    % Add Info to Figure
    % ------------------
    set(gca,'XLim',[-4.5 8]); set(gca,'YLim',[-0.2 0.45]);
    text(-3.8,0.4,['Npairs:' num2str(length(allr))],'FontSize',30,'Fontweight','normal');
    text(-3.8,0.35,sprintf('R: %0.2f, pval: %0.3f, preg: %0.3f',rval,pval,preg),'FontSize',30,'Fontweight','normal');
    text(-3.8,0.3,sprintf('Rsig: %0.2f, pval: %0.3f, Nsig: %g',rsigval,psigval,length(sigidx)),'FontSize',30,'Fontweight','normal');
    text(-3.8,0.25,sprintf('Rnosig: %0.2f, pval: %0.3f',rnosigval,pnosigval),'FontSize',30,'Fontweight','normal');

        
    figfile = [figdir,statename,'_ThetaCovVsRipCorr_',area]
    if savefig1==1,
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
    
    
    %RESTOR VALUES AFTER PLOT
    alltheta_peakcorr = alltheta_peakcorr_tmp;
    allr=allr_tmp;
    allp = allp_tmp;
    
    
    
    % 1b) Rip corrcoeff vs Total thetacorr for SWR BACKGROUND
    % -----------------------------------------------------
    % Get rid of NaNs / ONLY SIG CORR
    alltheta_peakcorr_tmp = alltheta_peakcorr; % Take Peak Corrln
    allr_rdm_tmp = allr_rdm; allp_rdmtmp = allp_rdm;
    
    rnan = find(isnan(allr_rdm) | isnan(alltheta_peakcorr));
    %tnan = find(isnan(alltheta_peakcorr));
    nocc = find(all_nsimulrdm<5);
    removeidxs = union(rnan,nocc);

    allr_rdm(removeidxs)=[]; alltheta_peakcorr(removeidxs)=[]; allp_rdm(removeidxs)=[];
    sigidx = find(allp_rdm<0.05); nosigidx = find(allp_rdm>=0.05);
    %allr_rdm = allr_rdm(sigidx); alltheta_peakcorr = alltheta_peakcorr(sigidx);
    
    figure; hold on; redimscreen_figforppt1;
    %set(gcf, 'Position',[205 136 723 446]);
    %xaxis = min(allr_rdm):0.1:max(allr_rdm);
    plot(alltheta_peakcorr, allr_rdm,['k.'],'MarkerSize',24);
    plot(alltheta_peakcorr(sigidx), allr_rdm(sigidx),['r.'],'MarkerSize',24);
    legend('Theta Cov vs SWR BCKGND CorrCoeff');
    
    title(sprintf('%s %s units', area, statename),'FontSize',tfont,'Fontweight','normal')
    title(sprintf('CA1-PFC SWR modulated pairs'),'FontSize',tfont,'Fontweight','normal')
    ylabel('SWR BCKGND Correlation','FontSize',yfont,'Fontweight','normal');
    xlabel('Peak Theta Covariance','FontSize',xfont,'Fontweight','normal');
    
    % Do stattistics on this popln plot
    [r_thetavsrip,p_thetavsrip] = corrcoef(allr_rdm, alltheta_peakcorr);
    [rsig,psig] = corrcoef(allr_rdm(sigidx), alltheta_peakcorr(sigidx));
    [rnosig,pnosig] = corrcoef(allr_rdm(nosigidx), alltheta_peakcorr(nosigidx));
    
    % Regression
    % -----------
    [b00,bint00,r00,rint00,stats00] = regress(allr_rdm', [ones(size(alltheta_peakcorr')) alltheta_peakcorr']);
    xpts = min(alltheta_peakcorr):0.01:max(alltheta_peakcorr);
    bfit00 = b00(1)+b00(2)*xpts;
    plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
    % Do regression after shifting data to make intercept 0
    % ------------------------------------------------------
    allr_rdm_0 = allr_rdm-mean(allr_rdm);
    alltheta_peakcorr_0 = alltheta_peakcorr-mean(alltheta_peakcorr);
    [b0,bint0,r0,rint0,stats0] = regress(allr_rdm_0',[ones(size(alltheta_peakcorr_0')) alltheta_peakcorr_0']);
    bfit0 = b0(1)+b0(2)*xpts;
    
    rval = roundn(r_thetavsrip(1,2),-2);
    pval = roundn(p_thetavsrip(1,2),-4);
    rsigval = roundn(rsig(1,2),-2);
    psigval = roundn(psig(1,2),-4);
    rnosigval = roundn(rnosig(1,2),-2);
    pnosigval = roundn(pnosig(1,2),-4);
    rsquare = roundn(stats0(1),-2);
    preg = roundn(stats0(3),-4);
    
    % Shuffling
    % ---------
    for n=1:1000
        idxs = randperm(length(allr_rdm));
        shuffle = allr_rdm(idxs);
        % Get corrcoeff of shuffle
        [rsh,psh] = corrcoef(alltheta_peakcorr, shuffle);
        r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
        % Get regression of shuffle after making intercept 0 / Or Not
        %shuffle_0 = shuffle - mean(shuffle);
        %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(alltheta_peakcorr_0')) alltheta_peakcorr_0']);
        [bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle', [ones(size(alltheta_peakcorr')) alltheta_peakcorr']);
        rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
        b_shuffle(n,:) = bsh;
    end
    prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
    % Get regression corresponding to 99 percentile
    idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
    idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
    bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
    plot(xpts,bfitsh,'k--','LineWidth',2);  % Theta vs Rip - 99% shuffle line
    
    % Add Info to Figure
    % ------------------
    set(gca,'XLim',[-4.5 8]); set(gca,'YLim',[-0.2 0.45]);
    text(-3.8,0.4,['Npairs:' num2str(length(allr_rdm))],'FontSize',30,'Fontweight','normal');
    text(-3.8,0.35,sprintf('R: %0.2f, pval: %0.3f, preg: %0.3f',rval,pval,preg),'FontSize',30,'Fontweight','normal');
    text(-3.8,0.3,sprintf('Rsig: %0.2f, pval: %0.3f, Nsig: %g', rsigval,psigval,length(sigidx)),'FontSize',30,'Fontweight','normal');
    text(-3.8,0.25,sprintf('Rnosig: %0.2f, pval: %0.3f',rnosigval,pnosigval),'FontSize',30,'Fontweight','normal');

    
    figfile = [figdir,statename,'_ThetaCovVsRipBCKCorr_',area]
    if savefig1==1,
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
end


alltheta_peakcorr = alltheta_peakcorr_tmp;



if 0
    % 1) Rip Coactivez vs Total thetacorr
    % -----------------------------------------------------
    % Get rid of NaNs
    rnan = find(isnan(allcoactivez));
    tnan = find(isnan(alltheta_peakcorr));
    removeidxs = [rnan,tnan];
    allcoactivez(removeidxs)=[]; alltheta_peakcorr(removeidxs)=[];
    
    figure; hold on; redimscreen_figforppt1;
    %set(gcf, 'Position',[205 136 723 446]);
    %xaxis = min(allcoactivez):0.1:max(allcoactivez);
    plot(alltheta_peakcorr, allcoactivez,[clr '.'],'MarkerSize',20);
    legend('Theta corr vs Rip CoactiveZ');
    
    title(sprintf('%s %s - Rip%s units', area, statename, kind),'FontSize',tfont,'Fontweight','normal')
    ylabel('Rip Resp CoactiveZ','FontSize',xfont,'Fontweight','normal');
    xlabel('Total Theta Corr','FontSize',yfont,'Fontweight','normal');
    
    % Do stattistics on this popln plot
    [r_thetavsrip,p_thetavsrip] = corrcoef(allcoactivez, alltheta_peakcorr);
    
    % Regression
    % -----------
    [b00,bint00,r00,rint00,stats00] = regress(allcoactivez', [ones(size(alltheta_peakcorr')) alltheta_peakcorr']);
    xpts = min(alltheta_peakcorr):0.01:max(alltheta_peakcorr);
    bfit00 = b00(1)+b00(2)*xpts;
    plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
    % Do regression after shifting data to make intercept 0
    % ------------------------------------------------------
    allcoactivez_0 = allcoactivez-mean(allcoactivez);
    alltheta_peakcorr_0 = alltheta_peakcorr-mean(alltheta_peakcorr);
    [b0,bint0,r0,rint0,stats0] = regress(allcoactivez_0',[ones(size(alltheta_peakcorr_0')) alltheta_peakcorr_0']);
    bfit0 = b0(1)+b0(2)*xpts;
    
    rval = roundn(r_thetavsrip(1,2),-2);
    pval = roundn(p_thetavsrip(1,2),-4);
    rsquare = roundn(stats0(1),-2);
    preg = roundn(stats0(3),-4);
    
    % Shuffling
    % ---------
    for n=1:1000
        idxs = randperm(length(allcoactivez));
        shuffle = allcoactivez(idxs);
        % Get corrcoeff of shuffle
        [rsh,psh] = corrcoef(alltheta_peakcorr, shuffle);
        r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
        % Get regression of shuffle after making intercept 0 / Or Not
        %shuffle_0 = shuffle - mean(shuffle);
        %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(alltheta_peakcorr_0')) alltheta_peakcorr_0']);
        [bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle', [ones(size(alltheta_peakcorr')) alltheta_peakcorr']);
        rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
        b_shuffle(n,:) = bsh;
    end
    prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
    % Get regression corresponding to 99 percentile
    idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
    idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
    bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
    plot(xpts,bfitsh,'k--','LineWidth',2);  % Theta vs Rip - 99% shuffle line
    
    % Add Info to Figure
    % ------------------
    set(gca,'XLim',[-4.5 8]); set(gca,'YLim',[-4 5]);
    text(4,5,['Npairs:' num2str(length(allcoactivez))],'FontSize',30,'Fontweight','normal');
    text(4,5,sprintf('R: %0.2f, pval: %0.3f, preg: %0.3f',rval,pval,preg),'FontSize',30,'Fontweight','normal');
    
    figfile = [figdir,area,'_',statename,'_ThetaCorrVsRip',kind,'_CoactiveZ']
    if savefig1==1,
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
end


   alltheta_peakcorr = alltheta_peakcorr_tmp;






% ------------------------------------------------------------------
% COMBINING PLOTS ACROSS FILES - DO MANUALLY
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











