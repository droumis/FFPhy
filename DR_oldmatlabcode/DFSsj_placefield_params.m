
% From placefield1
% Get parameters for placefields for exp and con groups

clear; close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 1; % Figure Options - Plot individual cell place plots

savedir = '/mnt/data25/sjadhav/RippleInterruption/ProcessedData/';
savefile = [savedir 'PlaceFieldParamsn'];
%savefile = [savedir 'REGrp_PlaceFieldParamsn'];  % n is after std is changed for mapfields
%savefile = [savedir 'RCGrp_PlaceFieldParamsn'];

% Parameters
minabsvel = 3;  % cm/sec - Most conservative for runs and place fields
minlinvel = 5;

% Plot options
plotanimidx = []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days

binsize=2; % for linearized trajectories


% If runscript, run Datafilter and save data
if runscript == 1
    
    
    %Animal selection
    %-----------------------------------------------------
    Expanimals = {'REc,''REd','REe','REf'};
    Conanimals = {'RCa','RCb','RCc','RCd'};
    %Expanimals = {'REd','REe','REf'};
    %Conanimals = {'RCb','RCc','RE1'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    epochfilter = 'isequal($type, ''run'')';
    epochfiltersl = 'isequal($type, ''sleep'')'; % for mean rate
    
    % Cell filter
    % -----------
    placecellfilter = 'strcmp($tag, ''PyrSR'')';
    
    % Time filter
    % -----------
    
    % Either use tetlist for riple detection for ripfilter,
    % or use ripple tet filter based on tag in tetinfo marking elec as ripple tet
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Abs linear velocity(thrs 5)/ velocity(thrs 3) time filter
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getvelpos', '(($absvel >= 3))'}, {'DFTFsj_getriptimes','($nripples == 0)',tetlist,'minthresh',2} }
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getvelpos', '(($absvel >= 3))'}, {'DFTFsj_getriptimes','($nripples == 0)',[],'tetfilter',riptetfilter,'minthresh',2} }
    % Linear velocity time filter implemented in getlinstate
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getriptimes','($nripples == 0)',tetlist,'minthresh',2} };
    timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',2} };
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    Expplaf = createfilter('animal',Expanimals,'days',dayfilter,'epochs',epochfilter,'cells',placecellfilter,'excludetime', timefilter,'iterator', iterator);
    Conplaf = createfilter('animal',Conanimals,'days',dayfilter,'epochs',epochfilter,'cells',placecellfilter,'excludetime', timefilter,'iterator', iterator);
    
    Expslef = createfilter('animal',Expanimals,'days',dayfilter,'epochs',epochfiltersl,'cells',placecellfilter,'iterator', iterator);
    Conslef = createfilter('animal',Conanimals,'days',dayfilter,'epochs',epochfiltersl,'cells',placecellfilter,'iterator', iterator);
    
    % Set analysis function
    % ----------------------
    % combine utilities from calctotalmeanrate and DFAsj_calcoverlap/sj_calcoverlap. Really dont need linpos
    % Relevant stuff also in DFAsj_peakdistance_linpos and DFAsj_peakdistance_traj
    Expplaf = setfilterfunction(Expplaf, 'DFAsj_getplacefieldparams', {'spikes', 'linfields','mapfields'}, 'binsize', 2);
    Conplaf = setfilterfunction(Conplaf, 'DFAsj_getplacefieldparams', {'spikes', 'linfields','mapfields'}, 'binsize', 2);
    
    Expslef = setfilterfunction(Expslef, 'DFAsj_calctotalmeanrate', {'spikes'}); % for sleep
    Conslef = setfilterfunction(Conslef, 'DFAsj_calctotalmeanrate', {'spikes'});
    
    % Run analysis
    % ------------
    Expplaf = runfilter(Expplaf);  
    Conplaf = runfilter(Conplaf);  
    
    Expslef = runfilter(Expslef);  
    Conslef = runfilter(Conslef);  
    
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

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/PlaceFields/Params/';
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

% Saving processed data and figures
savefig1 = 1;
%saveprocdatafile = [savedir,'PlaceFieldParamData'];
%saveprocdata = 0; saveproctag = 'exp';

% ---------------------------------------

str=['Con';'Exp'];
for g = 1:size(str,1)   % Do Exp and Con groups separately
    
  % Anim
    if ~isempty(plotanimidx)
        useanim = plotanimidx;
    else
        useanim = eval(['1:length(',str(g,:),'plaf);']); 
    end 
    
    % Days
    if ~isempty(plotdays)
        days = plotdays;
    else
        days = eval(['unique(',str(g,:),'plaf(1).epochs{1}(:,1));']);   % all days get from struct
    end
 
    allepochs = eval(['unique(',str(g,:),'plaf(1).epochs{1}(:,2));']); 
    totanim = length(useanim);
 
    % saving across animals and days. Can keep track of animals and days as well if I want to
    cntallcells=0; cntallcellssl=0;
    for anidx=1:totanim % Across animals
        an=useanim(anidx);
        
        % For run
        currindex=[]; 
        eval(['nidxs_an = length(',str(g,:),'plaf(an).output{1});']);
        for i=1:nidxs_an
            eval(['currindex(i,:) = ',str(g,:),'plaf(an).output{1}(i).index;']);
        end
        ep2s = find(currindex(:,2)==2); 
        ep4s = find(currindex(:,2)==4); % should have same number of both
        
        for n=1:length(ep2s)
            cntallcells = cntallcells+1;
            % Epoch2
            eval([str(g,:),'meanrate(cntallcells,1) =',str(g,:),'plaf(an).output{1}(ep2s(n)).meanrate;']);
            eval([str(g,:),'peakrate(cntallcells,1) =',str(g,:),'plaf(an).output{1}(ep2s(n)).peakrate;']);
            eval([str(g,:),'peakrate(cntallcells,1) =',str(g,:),'plaf(an).output{1}(ep2s(n)).peakrate;']);
            
            eval([str(g,:),'total_fracunder1(cntallcells,1) =',str(g,:),'plaf(an).output{1}(ep2s(n)).total_fracunder1;']);
            eval([str(g,:),'total_fracunder3(cntallcells,1) =',str(g,:),'plaf(an).output{1}(ep2s(n)).total_fracunder3;']);
            eval([str(g,:),'total_lthunder1(cntallcells,1) =',str(g,:),'plaf(an).output{1}(ep2s(n)).total_lthunder1;']);
            eval([str(g,:),'total_lthunder3(cntallcells,1) =',str(g,:),'plaf(an).output{1}(ep2s(n)).total_lthunder3;']);
            eval([str(g,:),'total_trajlth(cntallcells,1) =',str(g,:),'plaf(an).output{1}(ep2s(n)).total_trajlth;']);
            
            eval([str(g,:),'fracunder1(cntallcells,1,:) =',str(g,:),'plaf(an).output{1}(ep2s(n)).fracunder1;']);
            eval([str(g,:),'fracunder3(cntallcells,1,:) =',str(g,:),'plaf(an).output{1}(ep2s(n)).fracunder3;']);
            eval([str(g,:),'lthunder1(cntallcells,1,:) =',str(g,:),'plaf(an).output{1}(ep2s(n)).lthunder1;']);
            eval([str(g,:),'lthunder3(cntallcells,1,:) =',str(g,:),'plaf(an).output{1}(ep2s(n)).lthunder3;']);
            eval([str(g,:),'trajlth(cntallcells,1,:) =',str(g,:),'plaf(an).output{1}(ep2s(n)).trajlth;']);
            
            eval([str(g,:),'trajdata{cntallcells}{1} =',str(g,:),'plaf(an).output{1}(ep2s(n)).trajdata;']);
            eval([str(g,:),'mapdata{cntallcells}{1} =',str(g,:),'plaf(an).output{1}(ep2s(n)).mapdata;']);
            
            % Epoch4
            eval([str(g,:),'meanrate(cntallcells,2) =',str(g,:),'plaf(an).output{1}(ep4s(n)).meanrate;']);
            eval([str(g,:),'peakrate(cntallcells,2) =',str(g,:),'plaf(an).output{1}(ep4s(n)).peakrate;']);
            eval([str(g,:),'peakrate(cntallcells,2) =',str(g,:),'plaf(an).output{1}(ep4s(n)).peakrate;']);
            
            eval([str(g,:),'total_fracunder1(cntallcells,2) =',str(g,:),'plaf(an).output{1}(ep4s(n)).total_fracunder1;']);
            eval([str(g,:),'total_fracunder3(cntallcells,2) =',str(g,:),'plaf(an).output{1}(ep4s(n)).total_fracunder3;']);
            eval([str(g,:),'total_lthunder1(cntallcells,2) =',str(g,:),'plaf(an).output{1}(ep4s(n)).total_lthunder1;']);
            eval([str(g,:),'total_lthunder3(cntallcells,2) =',str(g,:),'plaf(an).output{1}(ep4s(n)).total_lthunder3;']);
            eval([str(g,:),'total_trajlth(cntallcells,2) =',str(g,:),'plaf(an).output{1}(ep4s(n)).total_trajlth;']);
            
            eval([str(g,:),'fracunder1(cntallcells,2,:) =',str(g,:),'plaf(an).output{1}(ep4s(n)).fracunder1;']);
            eval([str(g,:),'fracunder3(cntallcells,2,:) =',str(g,:),'plaf(an).output{1}(ep4s(n)).fracunder3;']);
            eval([str(g,:),'lthunder1(cntallcells,2,:) =',str(g,:),'plaf(an).output{1}(ep4s(n)).lthunder1;']);
            eval([str(g,:),'lthunder3(cntallcells,2,:) =',str(g,:),'plaf(an).output{1}(ep4s(n)).lthunder3;']);
            eval([str(g,:),'trajlth(cntallcells,2,:) =',str(g,:),'plaf(an).output{1}(ep4s(n)).trajlth;']);
            
            eval([str(g,:),'trajdata{cntallcells}{1} =',str(g,:),'plaf(an).output{1}(ep4s(n)).trajdata;']);
            eval([str(g,:),'mapdata{cntallcells}{1} =',str(g,:),'plaf(an).output{1}(ep4s(n)).mapdata;']);
            
        end % end ep2s
        
        % For sleep - mean rate
        currindexsl=[];
        eval(['nidxs_ansl = length(',str(g,:),'slef(an).output{1});']);
        for i=1:nidxs_ansl
            eval(['currindexsl(i,:) = ',str(g,:),'slef(an).output{1}(i).index;']);
        end
        ep1s = find(currindexsl(:,2)==1);
        ep3s = find(currindexsl(:,2)==3);
        ep5s = find(currindexsl(:,2)==5); % should have same number for all epochs
        
        for n=1:length(ep1s)
            cntallcellssl = cntallcellssl+1;
            % Epoch1
            eval([str(g,:),'meanratesl(cntallcellssl,1) =',str(g,:),'slef(an).output{1}(ep1s(n)).meanrate;']);
            % Epoch3
            eval([str(g,:),'meanratesl(cntallcellssl,2) =',str(g,:),'slef(an).output{1}(ep3s(n)).meanrate;']);
            % Epoch5
            eval([str(g,:),'meanratesl(cntallcellssl,3) =',str(g,:),'slef(an).output{1}(ep5s(n)).meanrate;']);
        end
        
        
        
        
        %           Wrong
        %         for d=1:length(days)
        %             currday = days(d);
        %             for ep = 1:length(allepochs)
        %                 currep = allepochs(ep);
        %
        %                 index = eval(['find( (',str(g,:),'plaf(an).epochs{1}(:,1)==currday) & (',str(g,:),'plaf(an).epochs{1}(:,2)==currep) );']);
        %                 d,ep,index
        %                 eval([str(g,:),'meanrate(d,ep) =',str(g,:),'plaf(an).output{1}(index).meanrate;']);
        %                 eval([str(g,:),'peakrate(d,ep) =',str(g,:),'plaf(an).output{1}(index).peakrate;']);
        %                 eval([str(g,:),'peaktraj(d,ep) =',str(g,:),'plaf(an).output{1}(index).peaktraj(1);']); % Take only 1 of there are multiple
        %
        %                 eval([str(g,:),'total_fracunder1(d,ep) =',str(g,:),'plaf(an).output{1}(index).total_fracunder1;']);
        %                 eval([str(g,:),'total_fracunder3(d,ep) =',str(g,:),'plaf(an).output{1}(index).total_fracunder3;']);
        %                 eval([str(g,:),'total_lthunder1(d,ep) =',str(g,:),'plaf(an).output{1}(index).total_lthunder1;']);
        %                 eval([str(g,:),'total_lthunder3(d,ep) =',str(g,:),'plaf(an).output{1}(index).total_lthunder3;']);
        %                 eval([str(g,:),'total_trajlth(d,ep) =',str(g,:),'plaf(an).output{1}(index).total_trajlth;']);
        %
        %                 eval([str(g,:),'fracunder1(d,ep,:) =',str(g,:),'plaf(an).output{1}(index).fracunder1;']);
        %                 eval([str(g,:),'fracunder3(d,ep,:) =',str(g,:),'plaf(an).output{1}(index).fracunder3;']);
        %                 eval([str(g,:),'lthunder1(d,ep,:) =',str(g,:),'plaf(an).output{1}(index).lthunder1;']);
        %                 eval([str(g,:),'lthunder3(d,ep,:) =',str(g,:),'plaf(an).output{1}(index).lthunder3;']);
        %                 eval([str(g,:),'trajlth(d,ep,:) =',str(g,:),'plaf(an).output{1}(index).trajlth;']);
        %
        %                 eval([str(g,:),'trajdata{an}{d}{ep} =',str(g,:),'plaf(an).output{1}(index).trajdata;']);
        %                 eval([str(g,:),'mapdata{an}{d}{ep} =',str(g,:),'plaf(an).output{1}(index).mapdata;']);
        %
        %             end % ep
        %         end % day
        
        
    end % anim
end % grp



%****************************************
% Figures -
%*******************************************

% Mean across epochs: Mean rates / Peak rates, etc

Expmeanratesl = mean(Expmeanratesl,2); % for sleep
Expmeanrate = mean(Expmeanrate,2); % Hz - for runs only
Exppeakrate = mean(Exppeakrate,2); % Hz - for runs only
Exptotal_lthunder1 = mean(Exptotal_lthunder1,2)*binsize; % in cm
Exptotal_lthunder3 = mean(Exptotal_lthunder3,2)*binsize;
Exptotal_trajlth = mean(Exptotal_trajlth,2)*binsize; % Total Traj Lth - length of all 4 trajs. 
Exptotal_fracunder1 = mean(Exptotal_fracunder1,2);
Exptotal_fracunder3 = mean(Exptotal_fracunder3,2);

Conmeanratesl = mean(Conmeanratesl,2); % for sleep
Conmeanrate = mean(Conmeanrate,2);
Conpeakrate = mean(Conpeakrate,2);
Contotal_lthunder1 = mean(Contotal_lthunder1,2)*binsize;
Contotal_lthunder3 = mean(Contotal_lthunder3,2)*binsize;
Contotal_trajlth = mean(Contotal_trajlth,2)*binsize; % Total Traj Lth - length of all 4 trajs. 
Contotal_fracunder1 = mean(Contotal_fracunder1,2);
Contotal_fracunder3 = mean(Contotal_fracunder3,2);



if figopt1 == 1
    
    
    % Total Lth under 1
    % ------------------
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    bar(1,mean(Contotal_lthunder1),'b');
    bar(2,mean(Exptotal_lthunder1),'r');
    errorbar(1,mean(Contotal_lthunder1),sem(Contotal_lthunder1),'k','LineWidth',3);
    errorbar(2,mean(Exptotal_lthunder1),sem(Exptotal_lthunder1),'k','LineWidth',3);
    [h_lthunder1,p_lthunder1] = ttest2(Contotal_lthunder1,Exptotal_lthunder1); %h=0, p=0.74
    title('Total length under 1Hz (all 4 traj)');
    ylabel('Place Field lth (cm)');
    set(gca,'XTick',[1 2],'XTickLabel',{'Con';'Exp'});
    if savefig1==1,
        figfile = [figdir,'PlaceFieldLth_1Hz'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
            
    % Proportion Track Active under 1
    % --------------------------------
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    bar(1,mean(Contotal_fracunder1),'b');
    bar(2,mean(Exptotal_fracunder1),'r');
    errorbar(1,mean(Contotal_fracunder1),sem(Contotal_fracunder1),'k','LineWidth',3);
    errorbar(2,mean(Exptotal_fracunder1),sem(Exptotal_fracunder1),'k','LineWidth',3);
    [h_fracunder1,p_fracunder1] = ttest2(Contotal_fracunder1,Exptotal_fracunder1); %h=0, p=0.76
    title('Proportion track active under 1Hz (all 4 traj)');
    ylabel('Proportion track active');
    set(gca,'XTick',[1 2],'XTickLabel',{'Con';'Exp'});
    if savefig1==1,
        figfile = [figdir,'PlaceFieldFrac_1Hz'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
    % Total Lth under 3
    % ------------------
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    bar(1,mean(Contotal_lthunder3),'b');
    bar(2,mean(Exptotal_lthunder3),'r');
    errorbar(1,mean(Contotal_lthunder3),sem(Contotal_lthunder3),'k','LineWidth',3);
    errorbar(2,mean(Exptotal_lthunder3),sem(Exptotal_lthunder3),'k','LineWidth',3);
    [h_lthunder3,p_lthunder3] = ttest2(Contotal_lthunder3,Exptotal_lthunder3); %h0=, p=0.75
    title('Total length under 3Hz (all 4 traj)');
    ylabel('Place Field lth (cm)');
    set(gca,'XTick',[1 2],'XTickLabel',{'Con';'Exp'});
    if savefig1==1,
        figfile = [figdir,'PlaceFieldLth_3Hz'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    % Proportion Track Active under 3
    % --------------------------------
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    bar(1,mean(Contotal_fracunder3),'b');
    bar(2,mean(Exptotal_fracunder3),'r');
    errorbar(1,mean(Contotal_fracunder3),sem(Contotal_fracunder3),'k','LineWidth',3);
    errorbar(2,mean(Exptotal_fracunder3),sem(Exptotal_fracunder3),'k','LineWidth',3);
    [h_fracunder3,p_fracunder3] = ttest2(Contotal_fracunder3,Exptotal_fracunder3); %h=0, p=0.77
    title('Proportion track active under 3Hz (all 4 traj)');
    ylabel('Proportion track active');
    set(gca,'XTick',[1 2],'XTickLabel',{'Con';'Exp'});
    if savefig1==1,
        figfile = [figdir,'PlaceFieldFrac_3Hz'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
        
        
    % Traj Lth 
    % ------------------
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    bar(1,mean(Contotal_trajlth),'b');
    bar(2,mean(Exptotal_trajlth),'r');
    errorbar(1,mean(Contotal_trajlth),sem(Contotal_trajlth),'k','LineWidth',3);
    errorbar(2,mean(Exptotal_trajlth),sem(Exptotal_trajlth),'k','LineWidth',3);
    [h_trajlth,p_trajlth] = ttest2(Contotal_trajlth,Exptotal_trajlth); %h=0, p=0.77
    title('Total traj length (all 4 traj)');
    ylabel('Trajectory lth (cm)');
    set(gca,'XTick',[1 2],'XTickLabel',{'Con';'Exp'});
    if savefig1==1,
        figfile = [figdir,'TrajLths'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    % Peak rates
    % ----------
    rem=find(Conpeakrate<3); Conpeakrate(rem)=[]; Contotal_fracunder1(rem)=[]; Contotal_fracunder3(rem)=[];
    rem=find(Exppeakrate<3); Exppeakrate(rem)=[]; Exptotal_fracunder1(rem)=[]; Exptotal_fracunder3(rem)=[];
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    hCon=hist(Conpeakrate,[0:5:50]);
    hExp=hist(Exppeakrate,[0:5:50]);
    plot([0:5:50],hCon,'b.-','Linewidth',4,'MarkerSize',24);
    plot([0:5:50],hExp,'r.-','Linewidth',4,'MarkerSize',24);
    [h_peakrate,p_peakrate] = ttest2(Conpeakrate,Exppeakrate); %h=0, p=0.9325
    [h_peakratedist,p_peakratedist] = kstest2(hCon,hExp);
    title('Peak rate distr');
    ylabel('No of neurons');
    xlabel('Peak rate (Hz)');
    if savefig1==1,
        figfile = [figdir,'PeakrateDistr'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    % Peak rate cdf
    cdf_hCon = cumsum(hCon)./max(cumsum(hCon));
    cdf_hExp = cumsum(hExp)./max(cumsum(hExp));
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    plot([0:5:50],cdf_hCon,'bd-','Linewidth',3,'MarkerSize',10);
    plot([0:5:50],cdf_hExp,'ro-','Linewidth',3,'MarkerSize',10);
    [h_peakratecdf,p_peakratecdf] = kstest2(cdf_hCon,cdf_hExp); % p = 0.9852
    [p_peakratecdf_ranksum,h_peakratecdf_ranksum] = ranksum(cdf_hCon,cdf_hExp); % p = 0.8688
    title('Peak rate distr');
    ylabel('Cumulative Fraction');
    xlabel('Peak rate (Hz)');
    if savefig1==1,
        figfile = [figdir,'PeakrateCDF'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    % Proportion Track active under 1 distributions
    % ---------------------------------------------
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    hCon=hist(Contotal_fracunder1,[0:0.05:1]);
    hExp=hist(Exptotal_fracunder1,[0:0.05:1]);
    plot([0:0.05:1],hCon,'b.-','Linewidth',4,'MarkerSize',24);
    plot([0:0.05:1],hExp,'r.-','Linewidth',4,'MarkerSize',24);
    [h_propact,p_propact] = ttest2(Contotal_fracunder1,Exptotal_fracunder1); %h=0, p=0.7162
    [h_propactdist,p_propactdist] = kstest2(hCon,hExp);
    title('Prop Track Active distr');
    ylabel('No of neurons');
    xlabel('Prop under 1Hz');
    if savefig1==1,
        figfile = [figdir,'PropTrackActiveDistr'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    % Peak rate cdf
    cdf_hCon = cumsum(hCon)./max(cumsum(hCon));
    cdf_hExp = cumsum(hExp)./max(cumsum(hExp));
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    plot([0:0.05:1],cdf_hCon,'bd-','Linewidth',3,'MarkerSize',10);
    plot([0:0.05:1],cdf_hExp,'ro-','Linewidth',3,'MarkerSize',10);
    [h_propactcdf,p_propactcdf] = kstest2(cdf_hCon,cdf_hExp); % p = 0.9999
    [p_propactcdf_ranksum,h_propactcdf_ranksum] = ranksum(cdf_hCon,cdf_hExp); % p = 0.9571
    title('Proportion Track Active');
    ylabel('Cumulative Fraction');
    xlabel('Proportion under 1Hz');
    set(gca,'XLim',[0 0.6])
    if savefig1==1,
        figfile = [figdir,'PropTrackActiveCDF'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
    
    % Proportion Track Active under 1 vs Peak Rate
    % ---------------------------------------------
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    plot(Contotal_fracunder1,Conpeakrate,'bd','MarkerSize',8,'Linewidth',2);
    plot(Exptotal_fracunder1,Exppeakrate,'ro','MarkerSize',8,'Linewidth',2);
    
    % Regression
    % Option: shift data to have 0 intercept or not
    
    %Exppeakrate0=Exppeakrate-mean(Exppeakrate);
    %Exptotal_fracunder10=Exptotal_fracunder1-mean(Exptotal_fracunder1);
    [be,binte,re,rinte,statse] = regress(Exppeakrate,[ones(size(Exptotal_fracunder1)) Exptotal_fracunder1]);
    xpts = 0:0.05:max(Exptotal_fracunder1);
    bfite = be(1)+be(2)*xpts;
    plot(xpts,bfite,'r-','LineWidth',4); 
    
    Conpeakrate0=Conpeakrate-mean(Conpeakrate);
    Contotal_fracunder10=Contotal_fracunder1-mean(Contotal_fracunder1);
    [bc,bintc,rc,rintc,statsc] = regress(Conpeakrate,[ones(size(Contotal_fracunder1)) Contotal_fracunder1]);
    xpts = 0:0.05:max(Contotal_fracunder1);
    bfitc = bc(1)+bc(2)*xpts;
    plot(xpts,bfitc,'b-','LineWidth',4);  
   
    
    title('Peak rate vs prop track active');
    ylabel('Peak Rate (Hz)');
    xlabel('Proportion track active');
    if savefig1==1,
        figfile = [figdir,'PeakRateVsPropTrackActive1Hz'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
end % end figopt



keyboard;



















