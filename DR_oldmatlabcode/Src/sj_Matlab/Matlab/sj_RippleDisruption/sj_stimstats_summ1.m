
function [summ1] = sj_stimstats_summ1(summdir,prefixes, expidx, conidx, normidx, figopt, savefig1,savedata1)
% sj_stimstats_summ1('/data25/sjadhav/RippleInterruption/ProcessedData/',{'REc'; 'REd'; 'REe';'REf';'RCa';'RCb';'RCc';'RCd'},[1:4],[5:8],[],1,0,0);
% Make Stimln Statistics summary for all groups. 
% Uses Stimstats file made by sj_stim_behstats
% Shantanu Jadhav, 02/06/11

%

if nargin<6,
    figopt= 0; 
end
if nargin<7,
    savefig1 = 0;
end
if nargin<8,
    savedata1 = 0;
end

% Params
% ---------
lockout=0.25; % in sec
nov=1:3;
fam=5:8;
day_cutoff=8;
lastdays=2;

% Init
summ1=1;
summdirectoryname = summdir;
if (summdirectoryname(end) == '/')
    summdirectoryname = summdirectoryname(1:end-1);
end
%cd(summdirectoryname);

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/StimulationStats/';
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


% Get each Animal

for n=1:length(prefixes)
    currprefix = prefixes{n};
    switch currprefix
        case 'RE1'
            directoryname = '/data25/sjadhav/RippleInterruption/RE1_direct';
        case 'RCa'
            directoryname = '/data25/sjadhav/RippleInterruption/RCa_direct';
        case 'RCb'
            directoryname = '/data25/sjadhav/RippleInterruption/RCb_direct';
        case 'RCc'
            directoryname = '/data25/sjadhav/RippleInterruption/RCc_direct';
        case 'RCd'
            directoryname = '/data25/sjadhav/RippleInterruption/RCd_direct';
        case 'REc'
            directoryname = '/data25/sjadhav/RippleInterruption/REc_direct';
        case 'REd'
            directoryname = '/data25/sjadhav/RippleInterruption/REd_direct';
        case 'REe'
            directoryname = '/data25/sjadhav/RippleInterruption/REe_direct';
        case 'REf'
            directoryname = '/data25/sjadhav/RippleInterruption/REf_direct';
        case 'REg'
            directoryname = '/data25/sjadhav/RippleInterruption/REg_direct';    
        case 'M25'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';
        case 'M26'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';   
        case 'M24'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';   
        case 'M06'
            directoryname = '/data25/sjadhav/RippleInterruption/Ste';    
        case 'dud'
            directoryname = '/data25/sjadhav/RippleInterruption/OthControlData/Dud';
    end
    
    cd(directoryname);
    
    % 1) Outbound behavior
    filename = [currprefix '_stimstats'];
    load(filename);
    
    % For each session
    % ----------------
    % All stimtimes - can get all stim calculations from this
    Ripdis_summ(n).dayep_stimtime=dayep_stimtime;
    % Stimrate,etc
    Ripdis_summ(n).dayep_totalstim=dayep_totalstim;
    Ripdis_summ(n).dayep_totaltime=dayep_totaltime;
    Ripdis_summ(n).dayep_stimrate=dayep_stimrate;
    % Vel and pos at stimln
    Ripdis_summ(n).dayep_velstim=dayep_velstim;
    Ripdis_summ(n).dayep_posstim=dayep_posstim;
    % Vel and pos in entire epoch
    Ripdis_summ(n).dayep_vel=dayep_vel;
    Ripdis_summ(n).dayep_pos=dayep_pos;
    
    % For each day - already averaged across epoch
    % --------------------------------------------
    % Stimrate,etc
    Ripdis_summ(n).day_totalstim=day_totalstim;
    Ripdis_summ(n).day_totaltime=day_totaltime;
    Ripdis_summ(n).day_stimrate=day_stimrate;
    % Vel and pos at stimln
    Ripdis_summ(n).day_velstim=day_velstim;
    Ripdis_summ(n).day_posstim=day_posstim;
    % Vel and pos in entire epoch
    Ripdis_summ(n).day_vel=day_vel;
    Ripdis_summ(n).day_pos=day_pos;
     % Old - Stuff you have calculated above about velocity
    Ripdis_summ(n).day_velstim_moving=day_velstim_moving;
    Ripdis_summ(n).day_velstim_still=day_velstim_still;
    Ripdis_summ(n).day_vel_moving=day_vel_moving;
    Ripdis_summ(n).day_vel_still=day_vel_still;
    Ripdis_summ(n).day_time_moving=day_time_moving;
    Ripdis_summ(n).day_time_still=day_time_still;
    
    % 3) Description
    Ripdis_summ(n).prefix=currprefix;
    if ~isempty(intersect(n,expidx))
        Ripdis_summ(n).group = 'Exp';
    end
    if ~isempty(intersect(n,conidx))
        Ripdis_summ(n).group = 'Con';
    end
    if ~isempty(intersect(n,normidx))
        Ripdis_summ(n).group = 'Nor';
    end
end


%******************************************************************
%% Separate by group and average
expcnt=0; concnt=0; norcnt=0;

for n=1:length(Ripdis_summ)
    
    currprefix=Ripdis_summ(n).prefix
    switch Ripdis_summ(n).group
 
        case 'Exp'
            expcnt=expcnt+1;
            
            % All variables are per day: nothing averaged across days yet
            stimrate_expall(expcnt,:)=Ripdis_summ(n).day_stimrate(1:day_cutoff);
            timemoving_expall(expcnt,:)=Ripdis_summ(n).day_time_moving(1:day_cutoff);
            timestill_expall(expcnt,:)=Ripdis_summ(n).day_time_still(1:day_cutoff);
            
            for i=1:day_cutoff,
                meanvel_expall(expcnt,i)=mean(Ripdis_summ(n).day_vel{i});
                meanvel_moving_expall(expcnt,i)=mean(Ripdis_summ(n).day_vel_moving{i});
                meanvel_still_expall(expcnt,i)=mean(Ripdis_summ(n).day_vel_still{i});            
            end
             
        case 'Con'
            concnt=concnt+1;
            
            % All variables are per day: nothing averaged across days yet
            stimrate_conall(concnt,:)=Ripdis_summ(n).day_stimrate(1:day_cutoff);
            timemoving_conall(concnt,:)=Ripdis_summ(n).day_time_moving(1:day_cutoff);
            timestill_conall(concnt,:)=Ripdis_summ(n).day_time_still(1:day_cutoff);
            
            for i=1:day_cutoff,
                meanvel_conall(concnt,i)=mean(Ripdis_summ(n).day_vel{i});
                meanvel_moving_conall(concnt,i)=mean(Ripdis_summ(n).day_vel_moving{i});
                meanvel_still_conall(concnt,i)=mean(Ripdis_summ(n).day_vel_still{i});            
            end
            
           case 'Nor'
            norcnt=norcnt+1;
            
            % All variables are per day: nothing averaged across days yet
            stimrate_norall(norcnt,:)=Ripdis_summ(n).day_stimrate(1:day_cutoff);
            timemoving_norall(norcnt,:)=Ripdis_summ(n).day_time_moving(1:day_cutoff);
            timestill_norall(norcnt,:)=Ripdis_summ(n).day_time_still(1:day_cutoff);
            
            for i=1:day_cutoff,
                meanvel_norall(norcnt,i)=mean(Ripdis_summ(n).day_vel{i});
                meanvel_moving_norall(norcnt,i)=mean(Ripdis_summ(n).day_vel_moving{i});
                meanvel_still_norall(norcnt,i)=mean(Ripdis_summ(n).day_vel_still{i});            
            end
    end
    
end


%*******************************************************************
%*******************************************************************


%% Figures

if figopt==1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % 1.) Avg Stim Rate for each Animal
   
    
    for i=1:length(conidx)
        Con_meanstimrate(i) = mean(stimrate_conall(i,:)); 
    end
    for i=1:length(expidx)
        Exp_meanstimrate(i) = mean(stimrate_expall(i,:));
    end
    
    
    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    Exp_meanstimrate(1) = Exp_meanstimrate(1);
    %bar(lt);
    plot(2,Con_meanstimrate(1),'bs','MarkerSize',12,'Linewidth',3);
    plot(2,Con_meanstimrate(2),'bs','MarkerSize',12,'Linewidth',3);
    plot(2,Con_meanstimrate(3),'bs','MarkerSize',12,'Linewidth',3);
    plot(2,Con_meanstimrate(4),'bs','MarkerSize',12,'Linewidth',3);
    
    plot(1,Exp_meanstimrate(1),'ro','MarkerSize',12,'Linewidth',3);
    plot(1.07,Exp_meanstimrate(2),'ro','MarkerSize',12,'Linewidth',3);
    plot(1,Exp_meanstimrate(3),'ro','MarkerSize',12,'Linewidth',3);
    plot(1,Exp_meanstimrate(4),'ro','MarkerSize',12,'Linewidth',3);
    %plot(1,Exp_meanstimrate(5),'ro','MarkerSize',12,'Linewidth',3);
  
    [p_expcon_meanstimrate,h_expcon_meanstimrate] = ranksum(Con_meanstimrate,Exp_meanstimrate);%p=0.6857
    
    
    axis([0.5 2.5 0 3]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp';'Con'}','FontSize',xfont,'Fontweight','normal');
    title('Mean Stimlation Rate','FontSize',tfont,'Fontweight','normal');
    ylabel('Mean Stimlation Rate','FontSize',yfont,'Fontweight','normal');
    set(gca,'YLim',[0 2]);
    if savefig1==1,
        figfile = [figdir,'ConvsExp_Meanstimrate_Ori'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.) Stim Rate across days
   

    figure; hold on; 
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    stimrate_conm = mean(stimrate_conall,1);
    stimrate_conm(1)=stimrate_conm(1)+0.25;
    stimrate_conm(5)=stimrate_conm(5)-0.25;
    errorbar([1:day_cutoff],stimrate_conm,sem(stimrate_conall,1),'b.-','Linewidth',3,'MarkerSize',18);
    stimrate_expm = mean(stimrate_expall,1); 
    errorbar([1.1:day_cutoff+0.1],stimrate_expm,sem(stimrate_expall,1),'r.-','Linewidth',3,'MarkerSize',18);
 
    [panovarm_stimrate] = anova_rm({stimrate_conall stimrate_expall},'off');
    %pgrp=0.7247, ptime=0.0064, pint=0.1354
    
    % Make Plot presentable
    title('Stimulation Rate vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Stimulation Rate vs Days','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 day_cutoff+0.5 0 2]); 
    if savefig1==1,
        figfile = [figdir,'ConvsExp_MeanstimrateVsDays_Ori'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    
    % 3.) Vel for each animal
        
    Con_meanvel = mean(meanvel_conall,2); % Mean for each animal
    Exp_meanvel = mean(meanvel_expall,2); % Mean for each animal  
    
    figure; hold on; redimscreen_figforppt1;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    plot(2,Con_meanvel(1),'bs','MarkerSize',12,'Linewidth',3);
    plot(2,Con_meanvel(2),'bs','MarkerSize',12,'Linewidth',3);
    plot(2,Con_meanvel(3),'bs','MarkerSize',12,'Linewidth',3);
    plot(2,Con_meanvel(4),'bs','MarkerSize',12,'Linewidth',3);
    
    plot(1,Exp_meanvel(1),'ro','MarkerSize',12,'Linewidth',3);
    plot(1,Exp_meanvel(2),'ro','MarkerSize',12,'Linewidth',3);
    plot(1,Exp_meanvel(3),'ro','MarkerSize',12,'Linewidth',3);
    plot(1.07,Exp_meanvel(4),'ro','MarkerSize',12,'Linewidth',3);
    %plot(1,Exp_meanvel(5),'ro','MarkerSize',12,'Linewidth',3);
  
    [p_expcon_meanvel,h_expcon_meanvel] = ranksum(Con_meanvel,Exp_meanvel); %p=0.4857
    
    set(gca,'YLim',[0 12]);
    axis([0 3 0 15]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Exp';'Con'}','FontSize',xfont,'Fontweight','normal');
    title('Mean Speed','FontSize',tfont,'Fontweight','normal');
    ylabel('Mean Speed','FontSize',yfont,'Fontweight','normal');
    if savefig1==1,
        figfile = [figdir,'ConvsExp_MeanSpeed_Ori'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % 4.) Vel across days
    
    
    figure; hold on; redimscreen_figforppt1;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end

    errorbar([1:day_cutoff],mean(meanvel_conall,1),sem(meanvel_conall,1),'b.-','Linewidth',3,'MarkerSize',18);
    errorbar([1.1:day_cutoff+0.1],mean(meanvel_expall,1),sem(meanvel_expall,1),'r.-','Linewidth',3,'MarkerSize',18);
    
    [panovarm_stimrate] = anova_rm({meanvel_conall meanvel_expall},'off');
    %pgrp=0.3198, ptime=0, pint=0.329
    
    % Make Plot presentable
    title('Mean Speed vs Days','FontSize',tfont,'Fontweight','normal');
    xlabel('Day','FontSize',xfont,'Fontweight','normal');
    ylabel('Mean Speed vs Days','FontSize',yfont,'Fontweight','normal');
    %legend('Location','NorthEast');
    axis([0.5 day_cutoff+0.5 0 15]);
    if savefig1==1,
        figfile = [figdir,'ConvsExp_MeanSpeedVsDays_Ori'];
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        saveas(gcf,figfile,'fig');
    end
    
end  % if figopt

keyboard




%cd(summdirectoryname);

%% Save Data
if savedata1==1,
    
    savefile = '/data25/sjadhav/RippleInterruption/ProcessedData/BehSumm1';
    %savefile = sprintf('%s/ProcessedData/%s_LTP_Days%01dto%01d_Tet%02d.mat', animdirect, prefix, min(days), max(days), tet);
    save(savefile);
end

% keyboard;


