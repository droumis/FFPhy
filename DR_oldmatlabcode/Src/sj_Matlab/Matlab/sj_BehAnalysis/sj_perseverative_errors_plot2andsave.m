
function [day1_perserr] = sj_perseverative_errors_plot2andsave(summdir,prefixes, figopt, saveg1,savedata1)

%% sj_perseverative_errors_plot ('/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals',{'sjc';'RE1'},0,0,0);
% sj_perseverative_errors_plot ('/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals',{'sjc';'RE1';'RCa';'Cor';'dud';'Eig';'Fiv';'Sev';'Six';'ten'},0,0,0);
% Calls sj_day_findinbound for EACH day, and returns vector of
% inbound_logic (Correct incorrect), wells exited and entered, and number
% of inbound trajectories per day, AND PERSEVERATIVE ERROR STATS
% Saves In File For Current Animal
% Shantanu Jadhav, 03/30/10

%%

if nargin<3,
    figopt= 0; % Still to update plotting to account for new data structure
end

if nargin<4,
    saveg1 = 0;
end

if nargin<5,
    savedata1 = 0;
end

%%
directoryname = summdir;
if (directoryname(end) == '/')
    directoryname = directoryname(1:end-1);
end
cd(directoryname);
clr = {'b','g','c','m','y','k','r'};

set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
set(0,'defaultaxeslinewidth',2);

%% Get each Animal

for n=1:length(prefixes)
    
    currprefix = prefixes{n};
    filename = [currprefix '_perserr'];
    load(filename);
    
    Rip_inbound_summary(n).proportion_correct = frac_corr;
    Rip_inbound_summary(n).proportion_perserr = frac_perserr;
    Rip_inbound_summary(n).proportion_turnerr = frac_turnerr;
    Rip_inbound_summary(n).number_trials = ntrajs_perday;
    
    if n<=2,
        Rip_inbound_summary(n).group = 'Ripple Disruption';
    end
    
    if n==3,
        Rip_inbound_summary(n).group = 'Control Stimulation';
    end
    
    if n>3,
        Rip_inbound_summary(n).group = 'Old Controls';
    end
    
    
end


%% Load Steve's data

load steve_Wtrack_inbound_summary;

Steve_inbound_summary = inbound_summary;


%% Plot

if figopt==1
    % Day 1
    figure(1); hold on;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    plot(Rip_inbound_summary(1:2).proportion_perserr(1),Rip_inbound_summary(1:2).proportion_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');   %sjc and RE1
    plot(Rip_inbound_summary(3).proportion_perserr(1),Rip_inbound_summary(3).proportion_turnerr(1),'bd','MarkerSize',10,'MarkerFaceColor','b');   %sjc and RE1
    plot(Rip_inbound_summary(4:end).proportion_perserr(1),Rip_inbound_summary(4:end).proportion_turnerr(1),'ko','MarkerSize',10,'MarkerFaceColor','g');   %Old controls
    plot(Steve_inbound_summary(1:4).proportion_perserr(1),Steve_inbound_summary(1:4).proportion_turnerr(1),'ks','MarkerSize',10,'MarkerFaceColor','k');   %Steves controls
    plot(Steve_inbound_summary(5:end).proportion_perserr(1),Steve_inbound_summary(5:end).proportion_turnerr(1),'ro','MarkerSize',10,'MarkerFaceColor','r');   %Steves lesion animals
    
    
    
    % Day 2
    figure(2); hold on;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    plot(Rip_inbound_summary(1:2).proportion_perserr(2),Rip_inbound_summary(1:2).proportion_turnerr(2),'gd','MarkerSize',10,'MarkerFaceColor','g');   %sjc and RE1
    plot(Rip_inbound_summary(3).proportion_perserr(2),Rip_inbound_summary(3).proportion_turnerr(2),'bd','MarkerSize',10,'MarkerFaceColor','b');   %sjc and RE1
    plot(Rip_inbound_summary(4:end).proportion_perserr(2),Rip_inbound_summary(4:end).proportion_turnerr(2),'ko','MarkerSize',10,'MarkerFaceColor','g');   %Old controls
    plot(Steve_inbound_summary(1:4).proportion_perserr(2),Steve_inbound_summary(1:4).proportion_turnerr(2),'ks','MarkerSize',10,'MarkerFaceColor','k');   %Steves controls
    plot(Steve_inbound_summary(5:end).proportion_perserr(2),Steve_inbound_summary(5:end).proportion_turnerr(2),'ro','MarkerSize',10,'MarkerFaceColor','r');   %Steves lesion animals
    
    % Day 3
    figure(3); hold on;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    plot(Rip_inbound_summary(1:2).proportion_perserr(3),Rip_inbound_summary(1:2).proportion_turnerr(3),'gd','MarkerSize',10,'MarkerFaceColor','g');   %sjc and RE1
    plot(Rip_inbound_summary(3).proportion_perserr(3),Rip_inbound_summary(3).proportion_turnerr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');   %sjc and RE1
    plot(Rip_inbound_summary(4:end).proportion_perserr(3),Rip_inbound_summary(4:end).proportion_turnerr(3),'ko','MarkerSize',10,'MarkerFaceColor','g');   %Old controls
    plot(Steve_inbound_summary(1:4).proportion_perserr(3),Steve_inbound_summary(1:4).proportion_turnerr(3),'ks','MarkerSize',10,'MarkerFaceColor','k');   %Steves controls
    plot(Steve_inbound_summary(5:end).proportion_perserr(3),Steve_inbound_summary(5:end).proportion_turnerr(3),'ro','MarkerSize',10,'MarkerFaceColor','r');   %Steves lesion animals
    
    
    
    %% Make Plots presentable
    
    figure(1); hold on;
    title('Day1: Inbound Errors','FontSize',24,'Fontweight','bold');
    xlabel('Proportion of perseverative errors','FontSize',16,'Fontweight','bold');
    ylabel('Proportion of turn-around errors','FontSize',16,'Fontweight','bold');
    %legend('Location','NorthEast');
    axis([0 1 0 1]); line([0 1], [1 0]);
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['PersErrorSummAllConDay1'],'jpg');
        saveas(gcf,['PersErrorSummAllConDay1'],'fig');
    end
    
    figure(2); hold on;
    title('Day2: Inbound Errors','FontSize',24,'Fontweight','bold');
    xlabel('Proportion of perseverative errors','FontSize',16,'Fontweight','bold');
    ylabel('Proportion of turn-around errors','FontSize',16,'Fontweight','bold');
    %legend('Location','NorthEast');
    axis([0 1 0 1]); line([0 1], [1 0]);
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['PersErrorSummAllConDay2'],'jpg');
        saveas(gcf,['PersErrorSummAllConDay2'],'fig');
    end
    
    
    % Day 3
    figure(3); hold on;
    title('Day3: Inbound Errors','FontSize',24,'Fontweight','bold');
    xlabel('Proportion of perseverative errors','FontSize',16,'Fontweight','bold');
    ylabel('Proportion of turn-around errors','FontSize',16,'Fontweight','bold');
    axis([0 1 0 1]); line([0 1], [1 0]);
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['PersErrorSummAllConDay3'],'jpg');
        saveas(gcf,['PersErrorSummAllConDay3'],'fig');
    end
    
    
    
    
    
    
    
    %% Make bar graphs
    
    % Day 1
    figure(11); hold on;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    bar(1,mean([day1_perserr(4:end)'; con_prop_perserr(:,1)]),'EdgeColor','k','FaceColor','none','Linewidth',2); % Steve's and other control animals
    plot(1*ones(length([day1_perserr(4:end)'; con_prop_perserr(:,1)])), [day1_perserr(4:end)'; con_prop_perserr(:,1)],'ko','MarkerFaceColor','k','MarkerSize',8);
    
    bar(2,mean(day1_perserr(3)),'EdgeColor','b','FaceColor','none','Linewidth',2); % RCa
    plot(2*ones(length(day1_perserr(3))), day1_perserr(3),'bo','MarkerFaceColor','b','MarkerSize',8);
    
    bar(3,mean([day1_perserr(1) day1_perserr(2)]),'EdgeColor','g','FaceColor','none','Linewidth',2); % RE1 and sjc
    plot(3*ones(length([day1_perserr(1) day1_perserr(2)])), [day1_perserr(1) day1_perserr(2)],'go','MarkerFaceColor','g','MarkerSize',8);
    
    bar(4,mean(les_prop_perserr(:,1)),'EdgeColor','r','FaceColor','none','Linewidth',2); % Steve's lesion animals
    plot(4*ones(length(les_prop_perserr(:,1))), les_prop_perserr(:,1),'ro','MarkerFaceColor','r','MarkerSize',8);
    
    title('Day1: Perseverative Errors','FontSize',24,'Fontweight','bold');
    ylabel('Proportion of perseverative errors','FontSize',16,'Fontweight','bold');
    xlabel('Animal Group','FontSize',16,'Fontweight','bold');
    %legend('Control','Control Stimulation','Ripple Disruption','Lesion','Location','NorthWest');
    
    % errorbar(1,mean([day1_perserr(4:end)'; con_prop_perserr(:,1)]),sem([day1_perserr(4:end)'; con_prop_perserr(:,1)]),'k');
    % errorbar(3,mean([day1_perserr(1) day1_perserr(2)]),sem([day1_perserr(1) day1_perserr(2)]),'g');
    % errorbar(4,mean(les_prop_perserr(:,1)),sem(les_prop_perserr(:,1)),'r');
    set(gca,'XTick',[1 2 3 4]);
    set(gca,'YLim',[0 1.1]);
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['PersErrorBarDay1'],'jpg');
        saveas(gcf,['PersErrorBarDay1'],'fig');
    end
    
    % Day 2
    figure(12); hold on;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    bar(1,mean([day2_perserr(4:end)'; con_prop_perserr(:,2)]),'EdgeColor','k','FaceColor','none','Linewidth',2); % Steve's and other control animals
    plot(1*ones(length([day2_perserr(4:end)'; con_prop_perserr(:,2)])), [day2_perserr(4:end)'; con_prop_perserr(:,2)],'ko','MarkerFaceColor','k','MarkerSize',8);
    
    bar(2,day2_perserr(3),'EdgeColor','b','FaceColor','none','Linewidth',2); % RCa
    plot(2*ones(length(day2_perserr(3))), day2_perserr(3),'bo','MarkerFaceColor','b','MarkerSize',8);
    
    bar(3,mean([day2_perserr(1) day2_perserr(2)]),'EdgeColor','g','FaceColor','none','Linewidth',2); % RE1 and sjc
    plot(3*ones(length([day2_perserr(1) day2_perserr(2)])), [day2_perserr(1) day2_perserr(2)],'go','MarkerFaceColor','g','MarkerSize',8);
    
    bar(4,mean(les_prop_perserr(:,2)),'EdgeColor','r','FaceColor','none','Linewidth',2); % Steve's lesion animals
    plot(4*ones(length(les_prop_perserr(:,2))), les_prop_perserr(:,2),'ro','MarkerFaceColor','r','MarkerSize',8);
    
    title('Day2: Perseverative Errors','FontSize',24,'Fontweight','bold');
    ylabel('Proportion of perseverative errors','FontSize',16,'Fontweight','bold');
    xlabel('Animal Group','FontSize',16,'Fontweight','bold');
    %legend('Control','Control Stimulation','Ripple Disruption','Lesion','Location','NorthWest');
    
    %errorbar(1,mean([day2_perserr(4:end)'; con_prop_perserr(:,2)]),sem([day2_perserr(4:end)'; con_prop_perserr(:,2)]),'k');
    %errorbar(3,mean([day2_perserr(1) day2_perserr(2)]),sem([day2_perserr(1) day2_perserr(2)]),'g');
    %errorbar(4,mean(les_prop_perserr(:,2)),sem(les_prop_perserr(:,2)),'r');
    set(gca,'XTick',[1 2 3 4]);
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['PersErrorBarDay2'],'jpg');
        saveas(gcf,['PersErrorBarDay2'],'fig');
    end
    
    % Day 3
    figure(13); hold on;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    bar(1,mean([day3_perserr(4:end)'; con_prop_perserr(:,3)]),'EdgeColor','k','FaceColor','none','Linewidth',2); % Steve's and other control animals
    plot(1*ones(length([day3_perserr(4:end)'; con_prop_perserr(:,3)])), [day3_perserr(4:end)'; con_prop_perserr(:,3)],'ko','MarkerFaceColor','k','MarkerSize',8);
    
    bar(2,day3_perserr(3),'EdgeColor','b','FaceColor','none','Linewidth',2); % RCa
    plot(2*ones(length(day3_perserr(3))), day3_perserr(3),'bo','MarkerFaceColor','b','MarkerSize',8);
    
    bar(3,mean([day3_perserr(1) day3_perserr(2)]),'EdgeColor','g','FaceColor','none','Linewidth',2); % RE1 and sjc
    plot(3*ones(length([day3_perserr(1) day3_perserr(2)])), [day3_perserr(1) day3_perserr(2)],'go','MarkerFaceColor','g','MarkerSize',8);
    
    bar(4,mean(les_prop_perserr(:,3)),'EdgeColor','r','FaceColor','none','Linewidth',2); % Steve's lesion animals
    plot(4*ones(length(les_prop_perserr(:,3))), les_prop_perserr(:,3),'ro','MarkerFaceColor','r','MarkerSize',8);
    
    title('Day3: Perseverative Errors','FontSize',24,'Fontweight','bold');
    ylabel('Proportion of perseverative errors','FontSize',16,'Fontweight','bold');
    xlabel('Animal Group','FontSize',16,'Fontweight','bold');
    %legend('Control','Control Stimulation','Ripple Disruption','Lesion','Location','NorthWest');
    
    % errorbar(1,mean([day3_perserr(4:end)'; con_prop_perserr(:,3)]),sem([day3_perserr(4:end)'; con_prop_perserr(:,3)]),'k');
    % errorbar(3,mean([day3_perserr(1) day3_perserr(2)]),sem([day3_perserr(1) day3_perserr(2)]),'g');
    % errorbar(4,mean(les_prop_perserr(:,3)),sem(les_prop_perserr(:,3)),'r');
    set(gca,'XTick',[1 2 3 4]);
    
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['PersErrorBarDay3'],'jpg');
        saveas(gcf,['PersErrorBarDay3'],'fig');
    end
    
    
    
end  % if plot




%% Save Data
if savedata1==1,
    
    savefile = '/data25/sjadhav/RippleInterruption/ProcessedData/PersErrSumm_9Jun2010';
    %savefile = sprintf('%s/ProcessedData/%s_LTP_Days%01dto%01d_Tet%02d.mat', animdirect, prefix, min(days), max(days), tet);
    save(savefile,'Rip_inbound_summary','Steve_inbound_summary');
end



