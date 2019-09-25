
% Similar to DFSsj_getstimstats_summ
% Some plots are made by sj_stimstatssumm1: Ntrials, Mean speed, Mean speed vs day

% DFSsj_getstimstats_summ calls sj_stim_behstats. 
% This script will call sj_behstats which does the same things but for Normal group as well
% and skips stimulation stuff

% Run the analysis as DFAsj_behstats
% Analysis function would load "pos","linpos", and would use the iterator "epoch_behaveanal" 

clear; close all;
runscript = 0; % if is 0, just loads saved processed data and plots stuff 
savedata = 0; % save data option - only works if runscript is also on
figopt1=1; % Figure Options - main summ figs 
dodistr=0; % Figures for speed distributions
dosepdays=0; % Separate figures for speed distributions on each day
dowelltime=1;

timestep = 1/29.97;

savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/';
%savefile = [savedir 'Behstats_summ'];
%savefile = [savedir 'Behstats_summ_rev']; % After Rev
savefile = [savedir 'Behstats_summ_rev2']; % revision2 - more animals in Exp group 


% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    Expanimals = {'REc','REd','REe','REf','REg','REh'};
    Conanimals = {'RCa','RCb','RCc','RCd'};
    Noranimals = {'RNa','RNb','RNc','RNd'};
    
    %Filter creation
    %--------------------------------------------------------
    % epoch filter
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    epochfilter = 'isequal($type, ''run'')';
    
    % iterator
    iterator = 'epochbehaveanal';
    
    % Create filter    
    Expbehf = createfilter('animal',Expanimals,'days',dayfilter,'epochs',epochfilter,'iterator', iterator);
    Conbehf = createfilter('animal',Conanimals,'days',dayfilter,'epochs',epochfilter,'iterator', iterator);
    Norbehf = createfilter('animal',Noranimals,'days',dayfilter,'epochs',epochfilter,'iterator', iterator);
    
    % Set filter function
    % --------------------
    Expbehf = setfilterfunction(Expbehf, 'DFAsj_behstats',{'pos','linpos'});
    Conbehf = setfilterfunction(Conbehf, 'DFAsj_behstats',{'pos','linpos'});
    Norbehf = setfilterfunction(Norbehf, 'DFAsj_behstats',{'pos','linpos'});
    
    % run analysis
    % -------------
    Expbehf = runfilter(Expbehf); 
    Conbehf = runfilter(Conbehf);
    Norbehf = runfilter(Norbehf);
    
    %--------------------- Finished Filter Function Run -------------------
    
    disp('Finished running filter');
    
    if savedata == 1
        clear figopt1 runscript savedata dodistr dosepdays dowelltime
        save(savefile);
    end
    
else
    
    load(savefile);
    
end % end runscript

if ~exist('savedata') % Just ran and saved data. Dont go below
    return
end

%--------------------- End Run Script -------------------

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/Behavior/BehControls/';
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
savefig1=0;
% --------Parameters-------------------------------
lockout=0.25; % in sec
nov=1:3;
fam=5:8;
day_cutoff=8;
lastdays=2;

% ---------Extract Data and Calculate ---------

str=['Exp';'Con';'Nor'];
days = unique(Expbehf(1).epochs{1}(:,1));
allepochs = unique(Expbehf(1).epochs{1}(:,2)); % ep 2 and 4

for g = 1:size(str,1)   % Do Exp and Con groups separately

    eval(['totanim = length(',str(g,:),'behf);']);
    for an=1:totanim % Across animals
         for d=1:length(days)
            currday = days(d);
            for ep = 1:length(allepochs)
                currep = allepochs(ep);
                index = eval(['find( (',str(g,:),'behf(an).epochs{1}(:,1)==currday) & (',str(g,:),'behf(an).epochs{1}(:,2)==currep) );']);
               
                eval([str(g,:),'vel{an}{d}{ep} =',str(g,:),'behf(an).output{1}(index).vel;']);
                eval([str(g,:),'pos{an}{d}{ep} =',str(g,:),'behf(an).output{1}(index).pos;']);
                
                eval([str(g,:),'wellpos{an}{d}{ep} =',str(g,:),'behf(an).output{1}(index).wellpos;']);
                eval([str(g,:),'intpos{an}{d}{ep} =',str(g,:),'behf(an).output{1}(index).intpos;']);
               
                eval([str(g,:),'vel_moving{an}{d}{ep} =',str(g,:),'behf(an).output{1}(index).vel_moving;']);
                eval([str(g,:),'vel_still{an}{d}{ep} =',str(g,:),'behf(an).output{1}(index).vel_still;']);
                eval([str(g,:),'time_moving{an}{d}{ep} =',str(g,:),'behf(an).output{1}(index).time_moving;']);
                eval([str(g,:),'time_still{an}{d}{ep} =',str(g,:),'behf(an).output{1}(index).time_still;']);    
                
                eval([str(g,:),'totaltimewells(an,d,ep) =',str(g,:),'behf(an).output{1}(index).totaltimewells;']); % Single number
                eval([str(g,:),'Ntrials(an,d,ep) =',str(g,:),'behf(an).output{1}(index).Ntrials;']); % Single number
                eval([str(g,:),'avgwelltime(an,d,ep) =',str(g,:),'behf(an).output{1}(index).avgwelltime;']);
                eval([str(g,:),'avgwelltime1(an,d,ep) =',str(g,:),'behf(an).output{1}(index).avgwelltime1;']);
                eval([str(g,:),'avgwelltime2(an,d,ep) =',str(g,:),'behf(an).output{1}(index).avgwelltime2;']);
                eval([str(g,:),'avgwelltime3(an,d,ep) =',str(g,:),'behf(an).output{1}(index).avgwelltime3;']);
            end % end epoch
            % Combine vel vectors for day across epochs
            d;
            eval([str(g,:),'vel_day{an}{d} =[',str(g,:),'vel{an}{d}{1};',str(g,:),'vel{an}{d}{2}];']); % column vector above
            % Get mean vel for day
            eval([str(g,:),'meanvel_day(an,d)= nanmean(',str(g,:),'vel_day{an}{d});']);
            % Combine posn matr for day across epochs
            eval([str(g,:),'pos_day{an}{d} =[',str(g,:),'pos{an}{d}{1};',str(g,:),'pos{an}{d}{2}];']); % nX2 matrix of x-y positions
            
            %Combine well visits for day across epochs
            eval([str(g,:),'avgwelltime_day(an,d) = mean([',str(g,:),'avgwelltime(an,d,1);',str(g,:),'avgwelltime(an,d,2)]);']);
            eval([str(g,:),'Ntrials_day(an,d) = mean([',str(g,:),'Ntrials(an,d,1);',str(g,:),'Ntrials(an,d,2)]);']);
            eval([str(g,:),'totaltimewells_day(an,d) = mean([',str(g,:),'totaltimewells(an,d,1);',str(g,:),'totaltimewells(an,d,2)]);']);
            eval([str(g,:),'avgwelltime1_day(an,d) = mean([',str(g,:),'avgwelltime1(an,d,1);',str(g,:),'avgwelltime1(an,d,2)]);']);
            eval([str(g,:),'avgwelltime2_day(an,d) = mean([',str(g,:),'avgwelltime2(an,d,1);',str(g,:),'avgwelltime2(an,d,2)]);']);
            eval([str(g,:),'avgwelltime3_day(an,d) = mean([',str(g,:),'avgwelltime3(an,d,1);',str(g,:),'avgwelltime3(an,d,2)]);']);
            
         end % end day
    end % end anim
    
    % Speed Distributions. like ripple size distribution    
    % -------------------------------------------------
    % histogram edges
    xhist =[0:1:50];
    % Initialize to store for alldays for stats
    eval([str(g,:),'vel_alldays=[];']);
    % Loop over each day and get separately
    for d=1:length(days)
        currday=days(d);
        % Initialize to store raw values for this day across animals for stats
        eval([str(g,:),'allvel_day{',num2str(currday),'}=[];']);
        
        % Get from each animal 
         for ani=1:totanim % Loop over anim
             an=ani;
             % Get current data
             eval(['currvel =',str(g,:),'vel_day{an}{currday};']);
             % Make histogram for current day and epoch 
             h = histc(currvel,xhist); hnorm=h./max(h);
             eval([str(g,:),'velhist_dayanim(currday,an,:) = hnorm;']);
             % Store raw values for day stats 
             eval([str(g,:),'allvel_day{',num2str(currday),'}=[',str(g,:),'allvel_day{',num2str(currday),'};currvel];']);
             % Store raw values for alldays for grp stats
             eval([str(g,:),'vel_alldays=[',str(g,:),'vel_alldays;currvel];']);
         end % end animal
    end % end day
    
    % Position matrices. Can I Combine across animals for REd,REe,REf and RCb,RCc,RCd - same tracks 
    % --------------------
    % Initialize to store for alldays 
    eval([str(g,:),'pos_alldays=[];']);
    % Loop over each day and get separately
    for d=1:length(days)
        currday=days(d);
        % Initialize to store raw values for this day across animals for stats
        eval([str(g,:),'allpos_day{',num2str(currday),'}=[];']);
        
        % Get from each animal 
        % for ani=2:totanim % Loop over anim
        % an=ani-1;
         for ani=2:2 % Loop over anim
             an=ani;
             % Get current data
             eval(['currpos =',str(g,:),'pos_day{an}{currday};']);
             
             % Store raw values for day across animals
             eval([str(g,:),'allpos_day{',num2str(currday),'}=[',str(g,:),'allpos_day{',num2str(currday),'};currpos];']); % nX2 matrix of x-y positions
             % Store raw values for alldays for grp stats
             eval([str(g,:),'pos_alldays=[',str(g,:),'pos_alldays;currpos];']);
          end % end animal
    end % end day
       
end % end str=Exp or Con
    

%****************************************
% Figures -
%****************************************

if figopt1 == 1
    
 
    
    % ************* Mean speed ******************
    % 1) Mean Speed for each animal
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    Exp_animmeanvel = mean(Expmeanvel_day,2);
    Con_animmeanvel = mean(Conmeanvel_day,2);
    Nor_animmeanvel = mean(Normeanvel_day,2);
    plot(1*ones(size(Exp_animmeanvel)), Exp_animmeanvel,'ro','MarkerSize',12,'Linewidth',3);
    plot(2*ones(size(Con_animmeanvel)), Con_animmeanvel,'bd','MarkerSize',12,'Linewidth',3);
    plot(3*ones(size(Nor_animmeanvel)), Nor_animmeanvel,'ks','MarkerSize',12,'Linewidth',3);
    [p_expcon_animmeanvel,h_expcon_animmeanvel] = ranksum(Con_animmeanvel,Exp_animmeanvel); %p=0.1714
    [p_expnor_animmeanvel,h_expnor_animmeanvel] = ranksum(Nor_animmeanvel,Exp_animmeanvel); %p=0.1143
    [p_connor_animmeanvel,h_connor_animmeanvel] = ranksum(Nor_animmeanvel,Con_animmeanvel); %p=0.1143
    axis([0.5 3.5 0 15]);
    set(gca,'XTick',[1 2 3 ],'XTickLabel',{'Exp';'Con';'Nor'}','FontSize',xfont,'Fontweight','normal');
    title('Mean Speed','FontSize',tfont,'Fontweight','normal');
    ylabel('Mean Speed','FontSize',yfont,'Fontweight','normal');
    if savefig1==1,
        figfile = [figdir,'0ConExpNor_MeanSpeed'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % 2) Avg Speed vs day
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    Exp_daymeanvel = mean(Expmeanvel_day,1); Exp_daymeanvelerr = sem(Expmeanvel_day,1);
    Con_daymeanvel = mean(Conmeanvel_day,1); Con_daymeanvelerr = sem(Conmeanvel_day,1);
    Nor_daymeanvel = mean(Normeanvel_day,1); Nor_daymeanvelerr = sem(Normeanvel_day,1);
    errorbar(days-0.1,Exp_daymeanvel,Exp_daymeanvelerr,'ro-','MarkerSize',12,'Linewidth',3);
    errorbar(days,Con_daymeanvel,Con_daymeanvelerr,'bd-','MarkerSize',12,'Linewidth',3);
    %errorbar(days+0.1,Nor_daymeanvel,Nor_daymeanvelerr,'ks-','MarkerSize',12,'Linewidth',3);
    [panovarm_meanvel] = anova_rm({Expmeanvel_day Conmeanvel_day},'off');
    [panovarm_meanvel] = anova_rm({Expmeanvel_day Conmeanvel_day Normeanvel_day},'off');
    %Only Con vs Exp: pgrp=0.1467, ptime=0, pint=0.2087
    % Con, Exp Nor: pgrp=0.0526, ptime=0, pint=0.1272
    
    % Make Plot presentable
    title('Mean Speed vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Mean Speed','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 max(days)+0.5 0 15]); 
    if savefig1==1,
        figfile = [figdir,'0ConExp_MeanSpeedVsDays'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
 
    
end % end figopt1


% *************************************
% Distributions - Combined for all days
% *************************************
     
if dodistr==1
   
     Conmeanvelhist = mean(squeeze(mean(Convelhist_dayanim,2))); % Mean along animals, then days
     Expmeanvelhist = mean(squeeze(mean(Expvelhist_dayanim,2))); % Mean along animals, then days
     Conerrvelhist = sem(squeeze(sem(Convelhist_dayanim,2))); 
     Experrvelhist = sem(squeeze(sem(Expvelhist_dayanim,2))); 
     
     % Vectors of histogram - for stats
     Convelhistvec_alldays = Convelhist_dayanim(:);
     Expvelhistvec_alldays = Expvelhist_dayanim(:);
     
     % All speeds during behavior
     figure; hold on;
     if forppr==1
         redimscreen_figforppr1;
     else
         redimscreen_figforppt1;
     end
     plot(xhist,Conmeanvelhist,'b.-','Linewidth',3,'MarkerSize',18);
     plot(xhist,Expmeanvelhist,'r.-','Linewidth',3,'MarkerSize',18);
     jbfill(xhist, [Conmeanvelhist+Conerrvelhist],[Conmeanvelhist-Conerrvelhist],'b','b',1,0.3);
     jbfill(xhist, [Expmeanvelhist+Experrvelhist],[Expmeanvelhist-Experrvelhist],'r','r',1,0.3);
     set(gca,'YLim',[0 1.02]); set(gca,'XLim',[-1 max(xhist)+1]);
     title(['Con and Exp - Speed distr during beh. All Days'],'FontSize',16,'Fontweight','normal');
     ylabel('Proportion','FontSize',yfont,'Fontweight','normal');
     xlabel('Speed(cm/s)','FontSize',xfont,'Fontweight','normal');
     % ----------- Stats -------------
     [h_vel,p_vel] = ttest2(Convel_alldays,Expvel_alldays);
     [h_velhist,p_velhist] = kstest2(Convelhistvec_alldays,Expvelhistvec_alldays);
     
end % end dodistr


% ***************************************
% Distributions - Separately for all days
% ***************************************

if dosepdays==1
    
    for n = 1:length(days)
            currday = days(n);
            Conmeanvelhist = squeeze(mean(Convelhist_dayanim(currday,:,:),2))'; % Mean along 2nd dimension of animals
            Conerrvelhist = squeeze(sem(Convelhist_dayanim(currday,:,:),2))'; % sem along 2nd dimension of animals
            Expmeanvelhist = squeeze(mean(Expvelhist_dayanim(currday,:,:),2))'; % Mean along 2nd dimension of animals
            Experrvelhist = squeeze(sem(Expvelhist_dayanim(currday,:,:),2))'; % sem along 2nd dimension of animals
            
            % Vectors of histogram for current day - for stats
            Convelhistvec_currday = squeeze(Convelhist_dayanim(currday,:,:)); Convelhistvec(currday,:) = Convelhistvec_currday(:);
            Expvelhistvec_currday = squeeze(Expvelhist_dayanim(currday,:,:)); Expvelhistvec(currday,:) = Expvelhistvec_currday(:);
            
            % Plot for current day
            figure; hold on;
            if forppr==1
                redimscreen_figforppr1;
            else
                redimscreen_figforppt1;
            end
            subplot(1,2,1); hold on;
            plot(xhist,Conmeanvelhist,'b.-','Linewidth',3,'MarkerSize',18);
            plot(xhist,Expmeanvelhist,'r.-','Linewidth',3,'MarkerSize',18);
            jbfill(xhist, [Conmeanvelhist+Conerrvelhist],[Conmeanvelhist-Conerrvelhist],'b','b',1,0.3);
            jbfill(xhist, [Expmeanvelhist+Experrvelhist],[Expmeanvelhist-Experrvelhist],'r','r',1,0.3);
            set(gca,'YLim',[0 1.02]); set(gca,'XLim',[-1 max(xhist)+1]);
            title(['Con and Exp - Speed distr during beh. Day ',num2str(currday)],'FontSize',16,'Fontweight','normal');
            ylabel('Proportion','FontSize',yfont,'Fontweight','normal');
            xlabel('Speed(cm/s)','FontSize',xfont,'Fontweight','normal');
            % ----------- Stats -------------
            [h_velday(n),p_velday(n)] = ttest2(Conallvel_day{currday},Expallvel_day{currday});
            [h_velhistday(n),p_velhistday(n)] = kstest2(Convelhistvec(currday,:),Expvelhistvec(currday,:));
            
      
    end % end days
end % end dosepdays


% AFTER REVIEW
% -------------


if dowelltime==1
    
      % ************* AVG Well Time ******************
    % 1) Avg Well Time for each animal
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
   
    Exp_animavgwelltime = mean(Expavgwelltime_day,2);
    Con_animavgwelltime = mean(Conavgwelltime_day,2);
    Nor_animavgwelltime = mean(Noravgwelltime_day,2);
    
    jit=[0,0,0,0.05]; jit3=[0,0.05,0,0]; jite=[jit,0.2,0.2];
    plot(1+jite, Exp_animavgwelltime,'ro','MarkerSize',12,'Linewidth',3);
    plot(2*ones(size(Con_animavgwelltime)), Con_animavgwelltime,'bd','MarkerSize',12,'Linewidth',3);
    %plot(3+jit3, Nor_animavgwelltime,'ks','MarkerSize',12,'Linewidth',3);
    [p_expcon_animavgwelltime,h_expcon_animavgwelltime] = ranksum(Con_animavgwelltime,Exp_animavgwelltime); %p=0.6095
    [p_expnor_animavgwelltime,h_expnor_animavgwelltime] = ranksum(Nor_animavgwelltime,Exp_animavgwelltime); %p=0.2571
    [p_connor_animavgwelltime,h_connor_animavgwelltime] = ranksum(Nor_animavgwelltime,Con_animavgwelltime); %p=0.4857
    axis([0.5 2.5 0 20]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp';'Con'}','FontSize',xfont,'Fontweight','normal');
    title('Avg Time spent at reward wells','FontSize',tfont,'Fontweight','normal');
    ylabel('Avg Time spent at reward wells','FontSize',yfont,'Fontweight','normal');
    if savefig1==1,
        figfile = [figdir,'0ConExp_AvgWellTime'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % 2) Avg Well Time vs day
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    Exp_dayavgwelltime = mean(Expavgwelltime_day,1); Exp_dayavgwelltimeerr = sem(Expavgwelltime_day,1);
    Con_dayavgwelltime = mean(Conavgwelltime_day,1); Con_dayavgwelltimeerr = sem(Conavgwelltime_day,1);
    Nor_dayavgwelltime = mean(Noravgwelltime_day,1); Nor_dayavgwelltimeerr = sem(Noravgwelltime_day,1);
    errorbar(days-0.1,Exp_dayavgwelltime,Exp_dayavgwelltimeerr,'ro-','MarkerSize',12,'Linewidth',3);
    errorbar(days,Con_dayavgwelltime,Con_dayavgwelltimeerr,'bd-','MarkerSize',12,'Linewidth',3);
    errorbar(days+0.1,Nor_dayavgwelltime,Nor_dayavgwelltimeerr,'ks-','MarkerSize',12,'Linewidth',3);
    [panovarm0_avgwelltime] = anova_rm({Expavgwelltime_day Conavgwelltime_day Noravgwelltime_day},'off');
    [panovarm_avgwelltime] = anova_rm({Expavgwelltime_day Conavgwelltime_day},'off');
    %Con, Exp, Nor pgrp=0.5533, ptime=0, pint=0.4223
    % Only Con Exp: pgrp=0.8315, ptime=0.0009, pint=0.24
    
    % Make Plot presentable
    title('Avg Well Time vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Avg Well Time','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 max(days)+0.5 0 20]); 
    if savefig1==1,
        figfile = [figdir,'0ConExp_AvgWellTimeVsDays'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
    % 3) AVG Well Time for each well separately
    % -----------------------------------------
      figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    Exp_animavgwelltime1 = mean(Expavgwelltime1_day,2);
    Con_animavgwelltime1 = mean(Conavgwelltime1_day,2);
    Exp_animavgwelltime2 = mean(Expavgwelltime2_day,2);
    Con_animavgwelltime2 = mean(Conavgwelltime2_day,2);
    Exp_animavgwelltime3 = mean(Expavgwelltime3_day,2);
    Con_animavgwelltime3 = mean(Conavgwelltime3_day,2);
    
    jit=[0,0,0,0.05]; jit2=[0.05,0,0,0]; jite=[jit,0.2,0.2];jit2e=[jit2,0.2,0.2];
    plot(0.8+jite, Exp_animavgwelltime2,'ro','MarkerSize',12,'Linewidth',3);
    plot(1.2*ones(size(Con_animavgwelltime2)), Con_animavgwelltime1,'bd','MarkerSize',12,'Linewidth',3);
    plot(2.8+jit2e, Exp_animavgwelltime1,'ro','MarkerSize',12,'Linewidth',3);
    plot(3.2*ones(size(Con_animavgwelltime1)), Con_animavgwelltime2,'bd','MarkerSize',12,'Linewidth',3);
    plot(4.8+jite, Exp_animavgwelltime3,'ro','MarkerSize',12,'Linewidth',3);
    plot(5.2*ones(size(Con_animavgwelltime3)), Con_animavgwelltime3,'bd','MarkerSize',12,'Linewidth',3);
    
    [p_expcon_animavgwelltime1,h_expcon_animavgwelltime1] = ranksum(Con_animavgwelltime1,Exp_animavgwelltime1); %p=0.4857
    [p_expcon_animavgwelltime2,h_expcon_animavgwelltime2] = ranksum(Con_animavgwelltime2,Exp_animavgwelltime2); %p=0.3429
    [p_expcon_animavgwelltime3,h_expcon_animavgwelltime3] = ranksum(Con_animavgwelltime3,Exp_animavgwelltime3); %p=1
    
    [p_expexp_animavgwelltime12,h_expexp_animavgwelltime12] = ranksum(Exp_animavgwelltime1,Exp_animavgwelltime2);%0.8857
    [p_expexp_animavgwelltime13,h_expexp_animavgwelltime13] = ranksum(Exp_animavgwelltime1,Exp_animavgwelltime3);%0.1143
    [p_expexp_animavgwelltime23,h_expexp_animavgwelltime23] = ranksum(Exp_animavgwelltime2,Exp_animavgwelltime3);%0.1143
    
    [p_concon_animavgwelltime12,h_concon_animavgwelltime12] = ranksum(Con_animavgwelltime1,Con_animavgwelltime2);%0.4857
    [p_concon_animavgwelltime13,h_concon_animavgwelltime13] = ranksum(Con_animavgwelltime1,Con_animavgwelltime3);%0.6843
    [p_concon_animavgwelltime23,h_concon_animavgwelltime23] = ranksum(Con_animavgwelltime2,Con_animavgwelltime3);%0.8857
    
    mean(Exp_animavgwelltime1); mean(Con_animavgwelltime1); %8.30, 9.07
    mean(Exp_animavgwelltime2); mean(Con_animavgwelltime2); %8.43, 11.35
    mean(Exp_animavgwelltime3); mean(Con_animavgwelltime3);  %11.97, 11.35
    
    axis([0 6 0 20]);
    set(gca,'XTick',[1 3 5],'XTickLabel',{'SideWell';'CtrWell';'SideWell'}','FontSize',xfont,'Fontweight','normal');
    title('Avg Time spent at reward wells','FontSize',tfont,'Fontweight','normal');
    ylabel('Avg Time spent at reward wells','FontSize',yfont,'Fontweight','normal');
    if savefig1==1,
        figfile = [figdir,'ConExp_AvgWellTimeSeparate'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
end










