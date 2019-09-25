
% From DFSsj_getriprate_run. Just choose one session to plot. So pick anim and then check over days. 
% No errors on ripple size curves

% Allows Plotting of each group separately
% DIO is controlled for

% This will do run, in parallel to DFSsj_getripraten which does sleep

clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 1; % Figure Options for Ripple Rate
figopt2 = 1; % Figure Options for Ripple Size - all days
dosepdays=1; % Plot sep days for rip size

ntet=1;
savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/';

% New: Use DFTF_getstimtimes to filter out stim times (set all those times to 0) and then use DFAsj_getriprate
% Filtering out stimtimes: also done in DFTFsj_getriptimes_nostim
%savefile = [savedir 'RippleRate_All_run_eg']; % Speed criterion of 10cm/s. 
savefile = [savedir 'RippleRate_All_run_eg2']; % No Speed criterion. 

% Plot options - like DFSsj_xcorrmeasures1, but modified. Pick only 1 animal for each group
 plotanimidx_Con =  []; % To pick animals for plotting 
 plotanimidx_Exp = [];
 plotdays_Con = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days
 plotdays_Exp = [];

% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    Expanimals = {'REc';'REd';'REe';'REf'};
    Conanimals = {'RCa';'RCb';'RCc';'RCd'};
    
    %Filter creation
    %--------------------------------------------------------
    
    % epoch filter
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    epochfilter{1} = ['isequal($type, ''run'')'];
    
    % Tet filter for ripple detection
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Filter out stimtimes
    % Exp
%     timefilter = {{'DFTFsj_getvelpos', '(($absvel <= 10))'},...
%         {'DFTFsj_getstimtimes','($nstim == 0)','tetfilter',riptetfilter,'timewin',0.2}}; % Default timewin is 0.1=100ms
     timefilter = {{'DFTFsj_getstimtimes','($nstim == 0)','tetfilter',riptetfilter,'timewin',0.15}}; % Default timewin is 0.1=100ms
    
    % Con
%     timefilter_con = {{'DFTFsj_getvelpos', '(($absvel <= 10))'},...
%         {'DFTFsj_getstimtimes','($nstim == 0)','tetfilter',riptetfilter,'timewin1',0.05,'timewin2',0.1}}; %asymmfilt
    timefilter_con = {{'DFTFsj_getstimtimes','($nstim == 0)','tetfilter',riptetfilter,'timewin1',0.075,'timewin2',0.15}}; % Default timewin is 0.1=100ms
    
    % iterator
    iterator = 'multitetrodeanal'; % / iterator = eeganal;
    
    % filter creation
    Expripf = createfilter('animal',Expanimals,'days',dayfilter,'epochs',epochfilter,'eegtetrodes',riptetfilter,'excludetime', timefilter,'iterator', iterator);
    Conripf = createfilter('animal',Conanimals,'days',dayfilter,'epochs',epochfilter,'eegtetrodes',riptetfilter,'excludetime', timefilter_con,'iterator', iterator);
 
    % set analysis function
    % Call get riprate after filtering out all stim times
    Expripf = setfilterfunction(Expripf, 'DFAsj_getriprate', {'ripplesep1'}, 'numtetrodes', ntet,'minthresh',3);
    Conripf = setfilterfunction(Conripf, 'DFAsj_getriprate', {'ripplesep1'}, 'numtetrodes', ntet,'minthresh',3);
    
    % run analysis
    Expripf = runfilter(Expripf);  % Ripple rate
    Conripf = runfilter(Conripf);
    
    %--------------------- Finished Filter Function Run -------------------
    
    disp('Finished running filter');
    
    if savedata == 1
        clear figopt1 figopt2 figoptsize dosepdays runscript savedata plotanimidx_Con plotanimidx_Exp plotdays_Con plotdays_Exp
        save(savefile);
    end
    
else
    
    load(savefile);
    
end  % end runscript

if ~exist('savedata')
    return
end


%--------------------- End Run Script -------------------


% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/RippleRate/Run/';
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

% ---------------------------------------
% % To Control Plotting, enter parameters here
% if ~isempty(plotdays)
%     usedays = plotdays;
% else
%     usedays = [];   % get from struct
% end
% -------------------------------------------

% Extract Data and Get ripple rate
str=['Con';'Exp'];
allepochs = unique(Expripf(1).epochs{1}(:,2));   % Run epochs 2 and 4

for g = 1:size(str,1)   % Do Exp and Con groups separately
    
    % Anim
    eval(['plotanimidx = plotanimidx_',str(g,:),';']);
    if ~isempty(plotanimidx)
        useanim = plotanimidx;
    else
        useanim = 3; % default
    end 
    if length(useanim)>1
        error('Enter only 1 animal in each group for plotting')
    end
    
    % Days
    eval(['plotdays = plotdays_',str(g,:),';']);
    if ~isempty(plotdays)
        days = plotdays;
    else
        days = eval(['unique(',str(g,:),'ripf(1).epochs{1}(:,1));']);   % all days get from struct
    end
   
    % Across animals: Inititalize to store ripplerate_ep, ripplepercenttime_ep, and totalstilltime_ep;
    for ep = 1:length(allepochs)
        currep = allepochs(ep);
        eval([str(g,:),'riprate_ep',num2str(currep),'=[];']);
        eval([str(g,:),'ripper_ep',num2str(currep),'=[];']);
        eval([str(g,:),'stilltime_ep',num2str(currep),'=[];']);
    end

    anidx=1;
    an=useanim(anidx);
    for d=1:length(days)
        currday = days(d);
        for ep = 1:length(allepochs)
            currep = allepochs(ep);
            
            index = eval(['find( (',str(g,:),'ripf(an).epochs{1}(:,1)==currday) & (',str(g,:),'ripf(an).epochs{1}(:,2)==currep) );']);
            eval([str(g,:),'riprate(d,ep) =',str(g,:),'ripf(an).output{1}(index).rip(1);']);
            eval([str(g,:),'ripper(d,ep) =',str(g,:),'ripf(an).output{1}(index).rip(2);']);
            eval([str(g,:),'stilltime(d,ep) =',str(g,:),'ripf(an).output{1}(index).rip(3);']);
            
            eval([str(g,:),'riprate_ep',num2str(currep),'(d) =',str(g,:),'ripf(an).output{1}(index).rip(1);']);
            eval([str(g,:),'ripper_ep',num2str(currep),'(d) =',str(g,:),'ripf(an).output{1}(index).rip(2);']);
            eval([str(g,:),'stilltime_ep',num2str(currep),'(d) =',str(g,:),'ripf(an).output{1}(index).rip(3);']);
            
            % Ripple  Size
            eval([str(g,:),'ripntet_ep',num2str(currep),'{d} =',str(g,:),'ripf(an).output{1}(index).ripntet;']);
            eval([str(g,:),'ripsize_ep',num2str(currep),'{d} =',str(g,:),'ripf(an).output{1}(index).ripsize;']);
            eval([str(g,:),'ripbaseline_ep',num2str(currep),'{d} =',str(g,:),'ripf(an).output{1}(index).rip_baseline;']); % Single number
            eval([str(g,:),'ripstd_ep',num2str(currep),'{d} =',str(g,:),'ripf(an).output{1}(index).rip_std;']); % Single number
            eval([str(g,:),'ripthresh_ep',num2str(currep),'{d} =',str(g,:),'ripf(an).output{1}(index).rip_thresh;']); % Single number
            
        end
    end
    
    
    % ***************   Ripple Rate  ******************
    % *************************************************
    % Calculate for ripple rate - not necessary. Done later during  plots as well
    runs = [2,4];
 
    % ***************   Ripple Size  ******************
    % *************************************************
    
    % Histogram Edges
    ripsizeedges = 3:0.5:12;
    % Make vector to combine across days: ripsize. Keep epochs separate for now. Easy to combine later
    for ep=1:length(allepochs)
        currep=allepochs(ep);
        eval([str(g,:),'ripsize_ep',num2str(currep),'_alldays=[];']); % combine across days. Keep epochs separate.
    end
    
    totanim=1;
    % Loop over each day and get separately
    for n = 1:length(days)
        currday = days(n);
        for ep=1:length(allepochs) % loop over the two run epochs
            
            currep=allepochs(ep);
            % Get current data
            eval(['currripsize =',str(g,:),'ripsize_ep',num2str(currep),'{currday};']);
            % ripple parameters
            eval(['curr_ripbaseline =',str(g,:),'ripbaseline_ep',num2str(currep),'{currday};']);
            eval(['curr_ripstd =',str(g,:),'ripstd_ep',num2str(currep),'{currday};']);
            
            
            % Make histogram for current day and epoch - ripsize
            h = histc(currripsize,ripsizeedges);
            h = cumsum(h); h = h./max(h);  % Cumulative proportion for current day and epoch
            % Store. You will take mean and sem across animals for plotting for day
            if length(find(isnan(h)))>0 % you get allnans when no ripp. Exp-Anim1Ep2-Days5and6/Anim3D3Ep2
                h(1)=0; h(2:length(ripsizeedges))=1;
            end
            eval([str(g,:),'ripsizehist_day_ep',num2str(currep),'(currday,:) = h;']);
            
            % Store vectors
            % Put raw values in vector for current day from all animals - Stats for each day based on raw stds
            eval([str(g,:),'ripsize_ep',num2str(currep),'_day{',num2str(currday),'}=currripsize;']);
            % Put raw values in vector for all days and all animals
            eval([str(g,:),'ripsize_ep',num2str(currep),'_alldays=[',str(g,:),'ripsize_ep',num2str(currep),'_alldays, currripsize];']);
            
        end % end loop over run epochs
    end  % end loop over days for ripple size
    
end   % end str = Exp and Con



%****************************************
% Figures -
%*******************************************


% Can combine epoch2 and epoch4 data or plot them separate
epcomb=0;
savefig1=0;

if figopt1 == 1
    
    % ********************** Ripple Rate ****************************
    % ***************************************************************
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    %  1) Bar Plot - Mean Ripple Rates
    
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    % Epoch 2 and 4 combined
    %-----------------------
    if (epcomb==1 | epcomb==0)
        bar(1,mean([Conriprate_ep2(:); Conriprate_ep4(:)]),'b');
        bar(2,mean([Expriprate_ep2(:); Expriprate_ep4(:)]),'r');
        errorbar(1,mean([Conriprate_ep2(:); Conriprate_ep4(:)]),sem([Conriprate_ep2(:); Conriprate_ep4(:)]),'k','LineWidth',2);
        errorbar(2,mean([Expriprate_ep2(:); Expriprate_ep4(:)]),sem([Expriprate_ep2(:); Expriprate_ep4(:)]),'k','LineWidth',2);
        % ----------- Stats -------------
        %Ndays  pts in each group: each day is a measure
        % Con against Exp epoch 2 and 4 combined
        [p_allriprate,h_allriprate] = ranksum([Conriprate_ep2(:); Conriprate_ep4(:)],[Expriprate_ep2(:); Expriprate_ep4(:)]);
        % Plot midway between the two bars
        if h_allriprate==1,
            mul = sign(mean([Conriprate_ep2(:); Conriprate_ep4(:)]));
            plot(1.5, (mean([Conriprate_ep2(:); Conriprate_ep4(:)]) + sem([Conriprate_ep2(:); Conriprate_ep4(:)]) + ...
                1.1*mul*(sem([Conriprate_ep2(:); Conriprate_ep4(:)]))), 'r*','MarkerSize',12);
        end
        % -------------------------
        title('Ripple Rate during Run: All Days');
        ylabel('Ripple Rate (Hz)');
        set(gca,'XTick',[1:2],'XTickLabel',{'Con';'Exp'},'FontSize',xfont,'Fontweight','normal');
        set(gca,'YLim',[0 1.3]);
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        if savefig1==1,
            figfile = [figdir,'RippleRateRunEg_EpComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
        
%     else
%         
%         %Epoch 2 and 4 separate
%         %-----------------------
%         runs = [2,4];
%         for ep=1:length(runs)
%             currep=runs(ep);
%             eval(['allmean_riprate(',num2str(ep),',:) = [Conallmean_riprateep',num2str(currep),';  Expallmean_riprateep',num2str(currep),'*0.7];'])
%         end
%         bar(allmean_riprate,'grouped');
%         errorbar(0.8,Conallmean_riprateep2,Conallerr_riprateep2,'k');
%         errorbar(1.2,Expallmean_riprateep2*0.7,Expallerr_riprateep2,'k');
%         errorbar(1.8,Conallmean_riprateep4,Conallerr_riprateep4,'k');
%         errorbar(2.2,Expallmean_riprateep4*0.7,Expallerr_riprateep4,'k');
%         set(gca,'XTick',[1:2],'XTickLabel',{'Ep2';'Ep4'},'FontSize',xfont,'Fontweight','normal');
%         %----------- Stats -------------
%         %Ndays  pts in each group: each day is a measure
%         %Con against Exp for each epoch
%         [p22_allriprate,h22_allriprate] = ranksum(Conriprate_ep2(:),Expriprate_ep2(:));
%         [p44_allriprate,h44_allriprate] = ranksum(Conriprate_ep4(:),Expriprate_ep4(:));
%         
%         %Plot midway between the two graphs
%         if h22_allriprate==1,
%             mul = sign(Conallmean_riprateep2);
%             plot(1, (Conallmean_riprateep2+Conallerr_riprateep2+1.1*mul*Conallerr_riprateep2), 'r*','MarkerSize',12);
%         end
%         if h44_allriprate==1,
%             mul = sign(Conallmean_riprateep4);
%             plot(2, (Conallmean_riprateep4+Conallerr_riprateep4+1.1*mul*Conallerr_riprateep4), 'r*','MarkerSize',12);
%         end
%         
%          if savefig1==1,
%             figfile = [figdir,'RippleRateRunEg_EpSep'];
%             print('-dpdf', figfile);
%             print('-djpeg', figfile);
%             saveas(gcf,figfile,'fig');
%         end
        
    end % end epcomb  
    
end % end figopt1





% ********************** Ripple Size ****************************
% ***************************************************************

if figopt2 == 1
    
    ripsizeedgesn =[3,ripsizeedges];
    
    rem=[];
    Conripsizehist_day_ep2(rem,:)=[]; Conripsizehist_day_ep4(rem,:)=[];
    Expripsizehist_day_ep2(rem,:)=[]; Expripsizehist_day_ep4(rem,:)=[];
    % Combine across days
    %-------------------------
    % Ep2
    Conmeanhist2 = mean(Conripsizehist_day_ep2); % Mean across days
    Conerrhist2 = sem(Conripsizehist_day_ep2);
    Expmeanhist2 = mean(Expripsizehist_day_ep2);
    Experrhist2 = sem(Expripsizehist_day_ep2);
    % Ep4
    Conmeanhist4 = mean(Conripsizehist_day_ep4);
    Conerrhist4 = sem(Conripsizehist_day_ep4);
    Expmeanhist4 = mean(Expripsizehist_day_ep4);
    Experrhist4 = sem(Expripsizehist_day_ep4);
    % Combine across epochs
    Conmeanhist = mean([Conmeanhist2;Conmeanhist4]);
    Conerrhist = sem([Conerrhist2;Conerrhist4]);
    Expmeanhist = mean([Expmeanhist2;Expmeanhist4]);
    Experrhist = sem([Experrhist2;Experrhist4]);
    % Vectors of histogram for current day - for stats
    Conhistvec_ep2_alldays = Conripsizehist_day_ep2(:);
    Exphistvec_ep2_alldays = Expripsizehist_day_ep2(:);
    Conhistvec_ep4_alldays = Conripsizehist_day_ep4(:);
    Exphistvec_ep4_alldays = Expripsizehist_day_ep4(:);
    
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Epoch 2 and 4 combined
    %-----------------------
    if epcomb==1
        plot(ripsizeedgesn,[0,Conmeanhist],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist+Conerrhist],...
            [Conmeanhist-Conerrhist],'b','b',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist+Experrhist],...
            [Expmeanhist-Experrhist],'r','r',1,0.3);
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize,p_ripsize] = ttest2([Conripsize_ep2_alldays, Conripsize_ep2_alldays],[Expripsize_ep2_alldays, Expripsize_ep4_alldays]);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist,h_ripsizehist] = ranksum([Conhistvec_ep2_alldays;Conhistvec_ep4_alldays],[Exphistvec_ep2_alldays;Exphistvec_ep4_alldays]);
        if h_ripsize == 1
            plot(7, 0.5, 'r*','MarkerSize',12);
        end
        text(7,0.5,['p = ',num2str(p_ripsize)],'FontSize',xfont);
        %---------------------------------
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Con vs Exp - Ripple Size during Run. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        if savefig1==1,
            figfile = [figdir,'RippleSizeEg_EpComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    else
        % Epoch 2 and 4 separate
        %-----------------------
        redimscreen_2horsubplots;
        
        subplot(1,2,1); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist2],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist2+Conerrhist2],...
            [Conmeanhist2-Conerrhist2],'b','b',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist2],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist2+Experrhist2],...
            [Expmeanhist2-Experrhist2],'r','r',1,0.3);
        % ----------- Stats -------------
        [h_ripsize_ep2,p_ripsize_ep2] = ttest2(Conripsize_ep2_alldays,Expripsize_ep2_alldays);
        [p_ripsizehist_ep2,h_ripsizehist_ep2] = ranksum(Conhistvec_ep2_alldays,Exphistvec_ep2_alldays);
        if h_ripsize_ep2 == 1
            plot(7, 0.5, 'r*','MarkerSize',12);
        end
        text(7,0.5,['p = ',num2str(p_ripsize_ep2)],'FontSize',xfont);
        %---------------------------------
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 11]);
        title(['Rip Size during Run. All Days: Ep2'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        
        subplot(1,2,2); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist4],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist4+Conerrhist4],...
            [Conmeanhist4-Conerrhist4],'b','b',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist4],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist4+Experrhist4],...
            [Expmeanhist4-Experrhist4],'r','r',1,0.3);
        % ----------- Stats -------------
        [h_ripsize_ep4,p_ripsize_ep4] = ttest2(Conripsize_ep4_alldays,Expripsize_ep4_alldays);
        [p_ripsizehist_ep4,h_ripsizehist_ep4] = ranksum(Conhistvec_ep4_alldays,Exphistvec_ep4_alldays);
        if h_ripsize_ep4 == 1
            plot(7, 0.5, 'r*','MarkerSize',12);
        end
        text(7,0.4,['p = ',num2str(p_ripsize_ep4)],'FontSize',xfont);
        %---------------------------------
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['All Days: Ep4'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
         if savefig1==1,
            figfile = [figdir,'RippleSizeEg_EpSep'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    end % end epcomb
    
end % end figopt2



if dosepdays==1
    
     ripsizeedgesn =[3,ripsizeedges];
   
    % Ripple Size - For each day separately
    %-------------------------
    days=1:size(Conripsizehist_day_ep2,1);
    %days=[1,2,3,5,8]; 
    for n = 1:length(days)
        currday = days(n);
        % Ep2
        Conmeanhist2 = Conripsizehist_day_ep2(currday,:); % Currday
        Expmeanhist2 = Expripsizehist_day_ep2(currday,:);
        % Ep4
        Conmeanhist4 = Conripsizehist_day_ep4(currday,:); 
        Expmeanhist4 = Expripsizehist_day_ep4(currday,:);
        % Combine across epochs
        Conmeanhist = mean([Conmeanhist2;Conmeanhist4]);
        Expmeanhist = mean([Expmeanhist2;Expmeanhist4]);
        % Vectors of histogram for current day - for stats
        Conhistvec_ep2(currday,:) = Conmeanhist2; 
        Exphistvec_ep2(currday,:) = Expmeanhist2; 
        Conhistvec_ep4(currday,:) = Conmeanhist4; 
        Exphistvec_ep4(currday,:) = Expmeanhist4; 
        
        figure; hold on;
        if forppr==1
            redimscreen_figforppr1;
        else
            redimscreen_figforppt1;
        end
        
        % Epoch 2 and 4 combined
        %-----------------------
        if epcomb==1
            plot(ripsizeedgesn,[0,Conmeanhist],'b.-','Linewidth',2,'MarkerSize',18);
            plot(ripsizeedgesn,[0,Expmeanhist],'r.-','Linewidth',2,'MarkerSize',18);
            % ----------- Stats -------------
            % Raw values of ripsize in std combined across animals for current day
            [h_dayripsize(n),p_dayripsize(n)] = ttest2([Conripsize_ep2_day{currday}, Conripsize_ep4_day{currday}],[Expripsize_ep2_day{currday}, Expripsize_ep4_day{currday}]);
            % ripsizehist in std combined across animals for current day
            [p_dayripsizehist(n),h_dayripsizehist(n)] = ranksum([Conhistvec_ep2(currday,:),Conhistvec_ep4(currday,:)],[Exphistvec_ep2(currday,:),Exphistvec_ep4(currday,:)]);
            if h_dayripsize(n) == 1
                plot(8.5, 0.5, 'r*','MarkerSize',12);
            end
            text(7,0.5,['p = ',num2str(p_dayripsize(n))],'FontSize',xfont);
            text(7,0.35,['NripCon = ',num2str(length([Conripsize_ep2_day{currday}, Conripsize_ep4_day{currday}]))],'FontSize',xfont);
            text(7,0.20,['NripExp = ',num2str(length([Expripsize_ep2_day{currday}, Expripsize_ep4_day{currday}]))],'FontSize',xfont);
            %---------------------------------
            set(gca,'YLim',[0 1]);
            set(gca,'XLim',[0 max(ripsizeedges)+1]);
            title(['Con vs Exp - Ripple Size during Run. Day ',num2str(n)],'FontSize',tfont,'Fontweight','normal');
            ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
            if savefig1==1,
                figfile = [figdir,'RippleSizeEg_EpComb_Day',num2str(n)];
                print('-dpdf', figfile);
                print('-djpeg', figfile);
                saveas(gcf,figfile,'fig');
            end
        else
            
            % Epoch 2 and 4 separate
            %-----------------------
            redimscreen_2horsubplots;
            
            subplot(1,2,1); hold on;
            plot(ripsizeedgesn,[0,Conmeanhist2],'b.-','Linewidth',2,'MarkerSize',18);
            plot(ripsizeedgesn,[0,Expmeanhist2],'r.-','Linewidth',2,'MarkerSize',18);
            % ----------- Stats -------------
            [h_dayripsize_ep2(n),p_dayripsize_ep2(n)] = ttest2(Conripsize_ep2_day{currday},Expripsize_ep2_day{currday});
            [p_dayripsizehist_ep2(n),h_dayripsizehist_ep2(n)] = ranksum(Conhistvec_ep2(currday,:),Exphistvec_ep2(currday,:));
            if h_dayripsize_ep2(n) == 1
                plot(7, 0.5, 'r*','MarkerSize',12);
            end
            text(7,0.5,['p = ',num2str(p_dayripsize_ep2(n))],'FontSize',xfont);
            text(7,0.4,['NripCon = ',num2str(length(Conripsize_ep2_day{currday}))],'FontSize',xfont);
            text(7,0.3,['NriprateCon = ',num2str(roundn(Conriprate(n,1),-2)),'Hz'],'FontSize',xfont);
            text(7,0.2,['NripExp = ',num2str(length(Expripsize_ep2_day{currday}))],'FontSize',xfont);
            text(7,0.1,['NriprateExp = ',num2str(roundn(Expriprate(n,1),-2)),'Hz'],'FontSize',xfont);
            %---------------------------------
            set(gca,'YLim',[0 1]);
            set(gca,'XLim',[0 11]);
            title(['Rip Size during Run. Day ',num2str(n),' Ep2'],'FontSize',tfont,'Fontweight','normal');
            ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
            
            subplot(1,2,2); hold on;
            plot(ripsizeedgesn,[0,Conmeanhist4],'b.-','Linewidth',2,'MarkerSize',18);
            plot(ripsizeedgesn,[0,Expmeanhist4],'r.-','Linewidth',2,'MarkerSize',18);
            % ----------- Stats -------------
            [h_dayripsize_ep4(n),p_dayripsize_ep4(n)] = ttest2(Conripsize_ep4_day{currday},Expripsize_ep4_day{currday});
            [p_dayripsizehist_ep4(n),h_dayripsizehist_ep4(n)] = ranksum(Conhistvec_ep4(currday,:),Exphistvec_ep4(currday,:));
            if h_dayripsize_ep4(n) == 1
                plot(7, 0.5, 'r*','MarkerSize',12);
            end
            text(7,0.5,['p = ',num2str(p_dayripsize_ep4(n))],'FontSize',xfont);
            text(7,0.4,['NripCon = ',num2str(length(Conripsize_ep4_day{currday}))],'FontSize',xfont);
            text(7,0.3,['NriprateCon = ',num2str(roundn(Conriprate(n,2),-1)),'Hz'],'FontSize',xfont);
            text(7,0.2,['NripExp = ',num2str(length(Expripsize_ep4_day{currday}))],'FontSize',xfont);
            text(7,0.1,['NriprateExp = ',num2str(roundn(Expriprate(n,2),-1)),'Hz'],'FontSize',xfont);
            %---------------------------------
            set(gca,'YLim',[0 1]);
            set(gca,'XLim',[0 max(ripsizeedges)+1]);
            title(['Day ',num2str(n),' Ep4'],'FontSize',tfont,'Fontweight','normal');
            ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
            if savefig1==1,
                figfile = [figdir,'RippleSizeEg_EpSep_Day',num2str(n)];
                print('-dpdf', figfile);
                print('-djpeg', figfile);
                saveas(gcf,figfile,'fig');
            end
        end % if epcomb
    end % end days
    
    
    keyboard;
    
    
end % end dosepdays





