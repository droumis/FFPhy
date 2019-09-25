
% Same as no sleep file - theta corr and rip corr coeff or coactivez
% Need to look across run and sleep epochs

clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - Individual cells

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';

% Sleep
% ------
% savefile = [savedir 'HP_thetaandripmodcorrsleep_corrandcoactz_PFC']; area = 'PFC'; clr = 'b'; % PFC
% savefile = [savedir 'HP_thetaandripunmodcorrsleep_corrandcoactz_PFC']; area = 'PFC'; clr = 'b'; % PFC
% savefile = [savedir 'HP_thetaandripmodcorrsleep_corrandcoactz_CA1']; area = 'CA1'; clr = 'r'; 
% savefile = [savedir 'HP_thetaandripunmodcorrsleep_corrandcoactz_CA1']; area = 'CA1'; clr = 'r'; 
% savefile = [savedir 'HP_thetaandripmodcorrsleep_corrandcoactz_CA1PFC']; area = 'CA1PFC'; clr = 'c'; 
% savefile = [savedir 'HP_thetaandripunmodcorrsleep_corrandcoactz_CA1PFC']; area = 'CA1PFC'; clr = 'c'; 

 savefile = [savedir 'HP_thetaandripmodcorrsleep_corrandcoactz_CA1allPFC']; area = 'CA1allPFC'; clr = 'c'; 
% savefile = [savedir 'HP_thetaandripunmodcorrsleep_corrandcoactz_CA1allPFC']; area = 'CA1allPFC'; clr = 'c';

savefig1=0;


% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPa','HPb'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    % Either Only do 1st w-track. 2 or 1 epochs per day
    % Or do Wtr1 and Wtr2, 2 epochs per day
    runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')';
    %sleepepochfilter = 'isequal($type, ''sleep'')'; % Only pre and post sleep marked as sleep
    sleepepochfilter = 'isequal($environment, ''postsleep'')'; % Only pre and post sleep marked as sleep
    
    % %Cell filter
    % %-----------
    % %PFC
    % %----
    %cellpairfilter = {'allcomb','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')'}; % Ripple mod
    %cellpairfilter = {'allcomb','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')'}; % Ripple unmod
    % %CA1
    % %----
    %cellpairfilter = {'allcomb','strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') && ($numspikes > 100) && strcmp($ripmodtag, ''y'')','strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') && ($numspikes > 100) && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') && strcmp($ripmodtag, ''n'')','strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') && strcmp($ripmodtag, ''n'')'};
    % %CA1-PFC
    % %--------
    %cellpairfilter = {'allcomb','strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') && ($numspikes > 100) && strcmp($ripmodtag, ''y'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') && ($numspikes > 100) && strcmp($ripmodtag, ''n'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')'};
    
    % %CA1all-PFCmodulated
    % %------------------------
     cellpairfilter = {'allcomb','strcmp($area, ''CA1'') || strcmp($area, ''iCA1'') && ($meanrate > 0.2)','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''y'')'};
    %cellpairfilter = {'allcomb','strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')','strcmp($area, ''PFC'') && strcmp($ripmodtag, ''n'')'};
    
    % Time filter - none.
    % -----------
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    timefilter_place = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6},...
        {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    
    % Iterator
    % --------
    iterator = 'singlecellanal';  % Have defined cellpairfilter. Can also use cellpair iterator with cell defn
    
    % Filter creation
    % ----------------
    
    % Ripp corrcoeff and coactivez
    modf = createfilter('animal',animals,'days',dayfilter,'epochs',sleepepochfilter, 'cellpairs',...
        cellpairfilter, 'iterator', iterator);
    
    thetaf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter, 'cellpairs',...
        cellpairfilter, 'excludetime', timefilter_place, 'iterator', iterator);
    
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    
    % % Need only the ripplemod structure for all of this.
    % % For calculating ripplealigned resp from scratch, you will need spikes, ripples, tetinfo, and pos
    modf = setfilterfunction(modf,'DFAsj_getripresp_corrandcoactz',{'ripplemod'}); %
    %modf = setfilterfunction(modf,'DFAsj_getripresp_corrandcoactz',{'ripplemod','cellinfo','spikes', 'ripples', 'tetinfo', 'pos'}); %
    
    %thetaf = setfilterfunction(thetaf,'DFAsj_calcxcorrmeasures', {'spikes'},'forripples',0);
    %thetaf = setfilterfunction(thetaf,'DFAsj_getthetacovariogram', {'spikes'});
    thetaf = setfilterfunction(thetaf,'DFAsj_getthetacrosscov', {'spikes'});
    
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
gatherdatafile = [savedir 'HP_thetaandripmodcorrsleep_corrandcoactz_CA1allPFC_gather']; area = 'CA1allPFC'; kind = 'mod'; state = 'sleep';
%gatherdatafile = [savedir 'HP_thetaandripunmodcorrsleep_corrandcoactz_CA1allPFC_gather']; area = 'CA1allPFC'; kind = 'unmod'; state = 'sleep';

% gatherdatafile = [savedir 'HP_thetaandripmodcorrsleep_corrandcoactz_PFC_gather']; area = 'PFC'; kind = 'mod'; state = 'sleep';
% gatherdatafile = [savedir 'HP_thetaandripunmodcorrsleep_corrandcoactz_PFC_gather']; area = 'PFC'; kind = 'unmod'; state = 'sleep';
% gatherdatafile = [savedir 'HP_thetaandripmodcorrsleep_corrandcoactz_CA1_gather']; area = 'CA1'; kind = 'mod'; state = 'sleep';
% gatherdatafile = [savedir 'HP_thetaandripunmodcorrsleep_corrandcoactz_CA1_gather']; area = 'CA1'; kind = 'unmod'; state = 'sleep';
% gatherdatafile = [savedir 'HP_thetaandripmodcorrsleep_corrandcoactz_CA1PFC_gather']; area = 'CA1PFC'; kind = 'mod'; state = 'sleep';
% gatherdatafile = [savedir 'HP_thetaandripunmodcorrsleep_corrandcoactz_CA1PFC_gather']; area = 'CA1PFC'; kind = 'unmod'; state = 'sleep';


if gatherdata
    
    % Parameters if any
    % -----------------
    
    % -------------------------------------------------------------
    
    cnt=0; 
    allanimindex=[]; allr=[]; allp = []; allr_rdm=[]; allp_rdm = []; allr_shuf=[]; allp_shuf = []; allcoactivez=[]; allcoactivez_rdm=[]; allcoactivez_shuf=[]; 
    allZcrosscov_runtheta=[]; allcrosscov_runtheta_totalcorr=[]; allrawcorr_runtheta=[];
    allZcrosscov_sm_runtheta=[]; allcrosscov_sm_runtheta_totalcorr=[]; allrawcorr_sm_runtheta=[];
    allNeventscorr_runtheta=[];allxcorr_runtheta=[]; allT_runtheta=[]; allp1p2_runtheta=[];
    runcorrtime=[];
    corrwin = 0.2; %Window for theta corrln
    
    % Rip Corr - In Sleep
    for an = 1:length(modf)
        for i=1:length(modf(an).output{1})
                cnt=cnt+1;
                anim_index{an}(cnt,:) = modf(an).output{1}(i).index;
                % Only indexes
                animindex=[an modf(an).output{1}(i).index]; % Put animal index in front
                allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
                % Data - Rip corr coeff and coact z
                allr(cnt) = modf(an).output{1}(i).r; 
                allr_rdm(cnt) = modf(an).output{1}(i).r_rdm;
                allr_shuf(cnt) = modf(an).output{1}(i).r_shuf;
                allp(cnt) = modf(an).output{1}(i).p; 
                allp_rdm(cnt) = modf(an).output{1}(i).p_rdm;
                allp_shuf(cnt) = modf(an).output{1}(i).p_shuf;
                
                allcoactivez(cnt) = modf(an).output{1}(i).coactivez; 
                allcoactivez_rdm(cnt) = modf(an).output{1}(i).coactivez_rdm;
                allcoactivez_shuf(cnt) = modf(an).output{1}(i).coactivez_shuf;
    
        end      
    end
    
    cnt=0; allanimindex_run=[];
    % Theta Corr - In RUn
    for an = 1:length(thetaf)
        for i=1:length(thetaf(an).output{1})
            cnt=cnt+1;
            anim_index_run{an}(cnt,:) = thetaf(an).output{1}(i).index;
            % Only indexes
            animindex_run=[an thetaf(an).output{1}(i).index]; % Put animal index in front
            allanimindex_run = [allanimindex_run; animindex_run]; % Collect all Anim Day Epoch Tet Cell Index
            
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
            
            % Calculate a number for theta corr - Total prob in -/+corrwin
            currthetacorr = allZcrosscov_runtheta(cnt,:);
            alltheta_totalcorr(cnt) = nansum(currthetacorr(bins_run))./length(bins_run); % per bin
            %alltheta_totalcorr(cnt) = nanmax(currthetacorr(bins_run));
            %alltheta_peaklag(cnt) = find (currthetacorr(bins_run) == nanmax(currthetacorr(bins_run)));
        end
    end
    
    % NEED TO CONSOLIDATE ACROSS EPOCHS TO COMPARE RUN AND SLEEP
    % ------------------------------------------------------------------------------
    
    % SLEEP DATA
    % ---------
    sleeppairoutput = struct;
    dummyindex=allanimindex;  % all anim-day-epoch-tet1-cell1-tet2-cell2 indices
    cntpairs=0;
    
    for i=1:size(allanimindex)
        animdaytetcell=allanimindex(i,[1 2 4 5 6 7]);
        ind=[];
        while rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7]))~=0          % collect all rows (epochs)
            ind = [ind rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7]))];        % finds the first matching row
            dummyindex(rowfind(animdaytetcell,dummyindex(:,[1 2 4 5 6 7])),:)=[0 0 0 0 0 0 0]; % after adding index, remove the corresponding row
            % so you could find the next one if it exists
        end
        
        % Gather everything for the current cell across epochs
        allrs=[]; allcoactivezs=[]; 
        for r=ind
            allrs = [allrs; allr(r)];
            allcoactivezs = [allcoactivezs; allcoactivez(r)];
        end
        
        if ~isempty(allrs)
            cntpairs=cntpairs+1;
            sleeppairoutput_idx(cntpairs,:)=animdaytetcell;
            sleeppairoutput(cntpairs).index=animdaytetcell; % This is anim-day-tet1-cell1-tet2-cell2. No epoch
            sleeppairoutput(cntpairs).allr = nanmean(allrs);
            sleeppairoutput(cntpairs).allcoactivez = nanmean(allcoactivezs);
        end
    end
    
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
            dummyindex(rowfind(animdaytetcell,dummyindex(:,[1 2 3 4 5 6 7])),:)=[0 0 0 0 0 0 0]; % after adding index, remove the corresponding row
            % so you could find the next one if it exists
        end
        
        % Gather everything for the current cell across epochs
        theta_normsmoothcorr=[]; theta_Nevents=[]; theta_totalcorr = [];
        for r=ind
            theta_normsmoothcorr = [theta_normsmoothcorr; allZcrosscov_sm_runtheta(r,:)];
            theta_Nevents = [theta_Nevents; allNeventscorr_runtheta(r)];
            theta_totalcorr = [theta_totalcorr; alltheta_totalcorr(r)];
        end
        
        if ~isempty(theta_totalcorr)
            cntpairs=cntpairs+1;
            runpairoutput_idx(cntpairs,:)=animdaytetcell;
            runpairoutput(cntpairs).index=animdaytetcell; % This is anim-day-tet1-cell1-tet2-cell2. No epoch
            runpairoutput(cntpairs).theta_corr = mean(theta_normsmoothcorr,1);
            runpairoutput(cntpairs).theta_totalcorr = nanmean(theta_totalcorr);
            runpairoutput(cntpairs).theta_Nevents = nansum(theta_Nevents);
        end
    end
    
 
    
    
    % Find corresponding Cell Pairs between Run and Sleep - Both have to Exist
    % --------------------------------------------------------------------
    
    cnt_runslpairs = 0; cnt_mismatch=0; 
    runsleeppair=struct; allfidx=[]; allftotalthetacorr=[]; allfr=[]; allfallcoactivez=[];
    for i=1:cntpairs % These are Run Pairs. Can also use while loop
        runidx = runpairoutput_idx(i,:);
        match = rowfind(runidx, sleeppairoutput_idx);
        if match~=0
            cnt_runslpairs = cnt_runslpairs+1;
            runsleeppair(cnt_runslpairs).index = runidx;
            runsleeppair(cnt_runslpairs).theta_totalcorr = runpairoutput(i).theta_totalcorr;
            runsleeppair(cnt_runslpairs).allr = sleeppairoutput(match).allr;
            runsleeppair(cnt_runslpairs).allcoactivez = sleeppairoutput(match).allcoactivez; 
            
            allfidx(cnt_runslpairs,:) = runidx; 
            allftheta_totalcorr(cnt_runslpairs) = runpairoutput(i).theta_totalcorr; 
            allfr(cnt_runslpairs) = sleeppairoutput(match).allr;
            allfallcoactivez(cnt_runslpairs) = sleeppairoutput(match).allcoactivez;  
        else
            
            cnt_mismatch = cnt_mismatch+1; % How many run pairs were not found in sleep? 
            
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

%corrRateRipples = nanmean(allp < 0.05)
%corrRateRdm = nanmean(allp_rdm < 0.05)
%corrRateShuf = nanmean(allp_shuf < 0.05)
%[r_realrdm p_realrdm] = ttest(allp<0.05,allp_rdm<0.05)
%[r_realshuf p_realshuf] = ttest(allp<0.05,allp_shuf<0.05)
%[r_rdmshuf p_rdmshuf] = ttest(allp_rdm<0.05,allp_shuf<0.05)

% ------------------
% Population Figures
% ------------------

forppr = 0;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Figures/RippleMod/Popln/';
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


allr = allfr;
allcoactivez =  allfallcoactivez;
alltheta_totalcorr = allftheta_totalcorr;


if 1   
      
    % 1) Rip corrcoeff vs Total thetacorr
    % -----------------------------------------------------
    % Get rid of NaNs
    alltheta_totalcorr_tmp = alltheta_totalcorr;
    rnan = find(isnan(allr));
    tnan = find(isnan(alltheta_totalcorr));
    removeidxs = [rnan,tnan];
    allr(removeidxs)=[]; alltheta_totalcorr(removeidxs)=[];
    
    figure; hold on; redimscreen_figforppt1;
    %set(gcf, 'Position',[205 136 723 446]);
    %xaxis = min(allr):0.1:max(allr);
    plot(alltheta_totalcorr, allr,[clr '.'],'MarkerSize',20);
    legend('Theta corr vs Rip Corrcoeff');

    title(sprintf('%s %s - Rip%s units', area, statename, kind),'FontSize',tfont,'Fontweight','normal')
    ylabel('Rip Resp Corr Coeff','FontSize',xfont,'Fontweight','normal');
    xlabel('Total Theta Corr','FontSize',yfont,'Fontweight','normal');

    % Do stattistics on this popln plot
    [r_thetavsrip,p_thetavsrip] = corrcoef(allr, alltheta_totalcorr);

    % Regression
    % -----------
    [b00,bint00,r00,rint00,stats00] = regress(allr', [ones(size(alltheta_totalcorr')) alltheta_totalcorr']);
    xpts = min(alltheta_totalcorr):0.01:max(alltheta_totalcorr);
    bfit00 = b00(1)+b00(2)*xpts;
    plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
    % Do regression after shifting data to make intercept 0
    % ------------------------------------------------------
    allr_0 = allr-mean(allr);
    alltheta_totalcorr_0 = alltheta_totalcorr-mean(alltheta_totalcorr);
    [b0,bint0,r0,rint0,stats0] = regress(allr_0',[ones(size(alltheta_totalcorr_0')) alltheta_totalcorr_0']);
    bfit0 = b0(1)+b0(2)*xpts;
    
    rval = roundn(r_thetavsrip(1,2),-2);
    pval = roundn(p_thetavsrip(1,2),-4);
    rsquare = roundn(stats0(1),-2);
    preg = roundn(stats0(3),-4);
    
    % Shuffling
    % ---------
    for n=1:1000
        idxs = randperm(length(allr));
        shuffle = allr(idxs);
        % Get corrcoeff of shuffle
        [rsh,psh] = corrcoef(alltheta_totalcorr, shuffle);
        r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
        % Get regression of shuffle after making intercept 0 / Or Not
        %shuffle_0 = shuffle - mean(shuffle);
        %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(alltheta_totalcorr_0')) alltheta_totalcorr_0']);
        [bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle', [ones(size(alltheta_totalcorr')) alltheta_totalcorr']);
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
    %set(gca,'XLim',[0 5]);
    text(4,0.3,['Npairs:' num2str(length(allr))],'FontSize',30,'Fontweight','normal');
    text(4,0.2,sprintf('R: %0.2f, pval: %0.3f, preg: %0.3f',rval,pval,preg),'FontSize',30,'Fontweight','normal');
        
    figfile = [figdir,area,'_',statename,'_ThetaCorrVsRip',kind,'_CorrCoeff']
    if savefig1==1,
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
end


alltheta_totalcorr = alltheta_totalcorr_tmp;


if 1
     % 1) Rip Coactivez vs Total thetacorr
    % -----------------------------------------------------
    % Get rid of NaNs
    rnan = find(isnan(allcoactivez));
    tnan = find(isnan(alltheta_totalcorr));
    removeidxs = [rnan,tnan];
    allcoactivez(removeidxs)=[]; alltheta_totalcorr(removeidxs)=[];
    
    figure; hold on; redimscreen_figforppt1;
    %set(gcf, 'Position',[205 136 723 446]);
    %xaxis = min(allcoactivez):0.1:max(allcoactivez);
    plot(alltheta_totalcorr, allcoactivez,[clr '.'],'MarkerSize',20);
    legend('Theta corr vs Rip CoactiveZ');

    title(sprintf('%s %s - Rip%s units', area, statename, kind),'FontSize',tfont,'Fontweight','normal')
    ylabel('Rip Resp CoactiveZ','FontSize',xfont,'Fontweight','normal');
    xlabel('Total Theta Corr','FontSize',yfont,'Fontweight','normal');

    % Do stattistics on this popln plot
    [r_thetavsrip,p_thetavsrip] = corrcoef(allcoactivez, alltheta_totalcorr);

    % Regression
    % -----------
    [b00,bint00,r00,rint00,stats00] = regress(allcoactivez', [ones(size(alltheta_totalcorr')) alltheta_totalcorr']);
    xpts = min(alltheta_totalcorr):0.01:max(alltheta_totalcorr);
    bfit00 = b00(1)+b00(2)*xpts;
    plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
    % Do regression after shifting data to make intercept 0
    % ------------------------------------------------------
    allcoactivez_0 = allcoactivez-mean(allcoactivez);
    alltheta_totalcorr_0 = alltheta_totalcorr-mean(alltheta_totalcorr);
    [b0,bint0,r0,rint0,stats0] = regress(allcoactivez_0',[ones(size(alltheta_totalcorr_0')) alltheta_totalcorr_0']);
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
        [rsh,psh] = corrcoef(alltheta_totalcorr, shuffle);
        r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
        % Get regression of shuffle after making intercept 0 / Or Not
        %shuffle_0 = shuffle - mean(shuffle);
        %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(alltheta_totalcorr_0')) alltheta_totalcorr_0']);
        [bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle', [ones(size(alltheta_totalcorr')) alltheta_totalcorr']);
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
    %set(gca,'XLim',[0 5]);
    text(4,5,['Npairs:' num2str(length(allcoactivez))],'FontSize',30,'Fontweight','normal');
    text(4,5,sprintf('R: %0.2f, pval: %0.3f, preg: %0.3f',rval,pval,preg),'FontSize',30,'Fontweight','normal');
        
    figfile = [figdir,area,'_',statename,'_ThetaCorrVsRip',kind,'_CoactiveZ']
    if savefig1==1,
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
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











