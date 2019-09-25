
function [day1_perserr] = sj_perseverative_errors_plot (summdir,prefixes, saveg1,saveg2,savedata1)

%% sj_perseverative_errors_plot ('/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals',{'sjc';'RE1'},0,0,0);
% sj_perseverative_errors_plot ('/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals',{'REb';'RE1';'RCa';'Cor';'dud';'Eig';'Fiv';'Sev';'Six';'ten'},0,0,0);
% Calls sj_day_findinbound for EACH day, and returns vector of
% inbound_logic (Correct incorrect), wells exited and entered, and number
% of inbound trajectories per day, AND PERSEVERATIVE ERROR STATS
% Saves In File For Current Animal
% Shantanu Jadhav, 03/30/10

%%
if nargin<3,
    saveg1 = 0;
end
if nargin<4
    saveg2 = 0; % Save new plots of Bar graph and new spreads
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
    day1_perserr(n) = frac_perserr(1); day1_turnerr(n) = frac_turnerr(1); 
    day1_nerrs(n) = nerr(1); day1_ntrajs(n) = ntrajs_perday(1);
    day1_fraccorr(n) = frac_corr(1);
    
    day2_perserr(n) = frac_perserr(2); day2_turnerr(n) = frac_turnerr(2);
    day2_nerrs(n) = nerr(2); day2_ntrajs(n) = ntrajs_perday(2);
    day2_fraccorr(n) = frac_corr(2);
    
    day3_perserr(n) = frac_perserr(3); day3_turnerr(n) = frac_turnerr(3);
    day3_nerrs(n) = nerr(3); day3_ntrajs(n) = ntrajs_perday(3);
    day3_fraccorr(n) = frac_corr(3);
    
%     if n==1,
%         for s=1:6,
%             cmd=sprintf('epoch%d_perserr(1) = frac_perserr_epoch(%d,1);',s,s); eval(cmd);
%             cmd=sprintf('epoch%d_turnerr(1) = frac_turnerr_epoch(%d,1);',s,s); eval(cmd);
%         end           
%     else
        epoch1_perserr(n) = frac_perserr_epoch(1,1); epoch1_turnerr(n) = frac_turnerr_epoch(1,1);
        epoch2_perserr(n) = frac_perserr_epoch(1,2); epoch2_turnerr(n) = frac_turnerr_epoch(1,2);
        epoch3_perserr(n) = frac_perserr_epoch(2,1); epoch3_turnerr(n) = frac_turnerr_epoch(2,1);
        epoch4_perserr(n) = frac_perserr_epoch(2,2); epoch4_turnerr(n) = frac_turnerr_epoch(2,2);
        epoch5_perserr(n) = frac_perserr_epoch(3,1); epoch5_turnerr(n) = frac_turnerr_epoch(3,1);
        epoch6_perserr(n) = frac_perserr_epoch(3,2); epoch6_turnerr(n) = frac_turnerr_epoch(3,2);
%     end
    
    %day4_perserr(n) = frac_perserr(4); day4_turnerr(n) = frac_turnerr(4);
    %day5_perserr(n) = frac_perserr(5); day5_turnerr(n) = frac_turnerr(5);
    
end


%% Load Steve's data

load steve_Wtrack_inbound_summary;


%% Plot 

% Day 1
figure(1); hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
plot(day1_perserr(1),day1_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');   %sjc
plot(day1_perserr(2),day1_turnerr(2),'cd','MarkerSize',10,'MarkerFaceColor','c');   %RE1
plot(day1_perserr(3),day1_turnerr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');   %RCa
plot(day1_perserr(4:end),day1_turnerr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');


% Day 2
figure(2); hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
plot(day2_perserr(1),day2_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');
plot(day2_perserr(2),day2_turnerr(2),'cd','MarkerSize',10,'MarkerFaceColor','c');
plot(day2_perserr(3),day2_turnerr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');   
plot(day2_perserr(4:end),day2_turnerr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');

% Day 3
figure(3); hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
plot(day3_perserr(1),day3_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');
plot(day3_perserr(2),day3_turnerr(2),'cd','MarkerSize',10,'MarkerFaceColor','c');
plot(day3_perserr(3),day3_turnerr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');   
plot(day3_perserr(4:end),day3_turnerr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');


% Plot Pers Error vs Prop Corr
figure(4); hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
plot(day1_perserr(1),day1_fraccorr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');   %sjc
plot(day1_perserr(2),day1_fraccorr(2),'gd','MarkerSize',10,'MarkerFaceColor','c');   %RE1
plot(day1_perserr(3),day1_fraccorr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');   %RCa
plot(day1_perserr(4:end),day1_fraccorr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');

figure(5); hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
plot(day2_perserr(1),day2_fraccorr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');   %sjc
plot(day2_perserr(2),day2_fraccorr(2),'gd','MarkerSize',10,'MarkerFaceColor','c');   %RE1
plot(day2_perserr(3),day2_fraccorr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');   %RCa
plot(day2_perserr(4:end),day2_fraccorr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');

figure(6); hold on;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
plot(day3_perserr(1),day3_fraccorr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');   %sjc
plot(day3_perserr(2),day3_fraccorr(2),'gd','MarkerSize',10,'MarkerFaceColor','c');   %RE1
plot(day3_perserr(3),day3_fraccorr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');   %RCa
plot(day3_perserr(4:end),day3_fraccorr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');



%% Overlay Steves Animals: Control = 1:4, Lesion = 5:10

cnt_con=0; cnt_les=0;
for i=1:length(inbound_summary)
    prop_perserr = inbound_summary(i).proportion_side_error;
    prop_turnerr = inbound_summary(i).proportion_revisit_error;
    prop_corr = inbound_summary(i).proportion_correct;
    
    if i<=4,
        cnt_con=cnt_con+1;
        con_prop_perserr(cnt_con,:)=prop_perserr;
        figure(1); hold on;
        plot(prop_perserr(1),prop_turnerr(1),'ks','MarkerSize',10,'MarkerFaceColor','k');
        figure(2); hold on;
        plot(prop_perserr(2),prop_turnerr(2),'ks','MarkerSize',10,'MarkerFaceColor','k');
        figure(3); hold on;
        plot(prop_perserr(3),prop_turnerr(3),'ks','MarkerSize',10,'MarkerFaceColor','k');
        
        figure(4); hold on;
        plot(prop_perserr(1),prop_corr(1),'ks','MarkerSize',10,'MarkerFaceColor','k');
        figure(5); hold on;
        plot(prop_perserr(2),prop_corr(2),'ks','MarkerSize',10,'MarkerFaceColor','k');
        figure(6); hold on;
        plot(prop_perserr(3),prop_corr(3),'ks','MarkerSize',10,'MarkerFaceColor','k');
        
%         figure(4); hold on;
%         plot(prop_perserr(4),prop_turnerr(4),'ks','MarkerSize',10,'MarkerFaceColor','k');
%         figure(5); hold on;
%         plot(prop_perserr(5),prop_turnerr(5),'ks','MarkerSize',10,'MarkerFaceColor','k');
    else
        cnt_les=cnt_les+1;
        les_prop_perserr(cnt_les,:)=prop_perserr;
        figure(1); hold on;
        plot(prop_perserr(1),prop_turnerr(1),'ro','MarkerSize',10,'MarkerFaceColor','r');
        figure(2); hold on;
        plot(prop_perserr(2),prop_turnerr(2),'ro','MarkerSize',10,'MarkerFaceColor','r');    
        figure(3); hold on;
        plot(prop_perserr(3),prop_turnerr(3),'ro','MarkerSize',10,'MarkerFaceColor','r');
        
        figure(4); hold on;
        plot(prop_perserr(1),prop_corr(1),'ro','MarkerSize',10,'MarkerFaceColor','r');
        figure(5); hold on;
        plot(prop_perserr(2),prop_corr(2),'ro','MarkerSize',10,'MarkerFaceColor','r');
        figure(6); hold on;
        plot(prop_perserr(3),prop_corr(3),'ro','MarkerSize',10,'MarkerFaceColor','r');
        
%         figure(4); hold on;
%         plot(prop_perserr(4),prop_turnerr(4),'ro','MarkerSize',10,'MarkerFaceColor','r');  
%         figure(5); hold on;
%         plot(prop_perserr(5),prop_turnerr(5),'ro','MarkerSize',10,'MarkerFaceColor','r');  
    end
    
end



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


% New plots with Pers Err vs Corr

figure(4); hold on;
title('Day1: Inbound Errors','FontSize',24,'Fontweight','bold');
xlabel('Proportion of perseverative errors','FontSize',16,'Fontweight','bold');
ylabel('Proportion of correct trajectories','FontSize',16,'Fontweight','bold');
%legend('Location','NorthEast');
axis([0 1 0 1]); line([0 1], [1 0]);

if saveg2==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['PersErrorvsCorrDay1'],'jpg');
    saveas(gcf,['PersErrorvsCorrDay1'],'fig');
end

figure(5); hold on;
title('Day2: Inbound Errors','FontSize',24,'Fontweight','bold');
xlabel('Proportion of perseverative errors','FontSize',16,'Fontweight','bold');
ylabel('Proportion of correct trajectories','FontSize',16,'Fontweight','bold');
%legend('Location','NorthEast');
axis([0 1 0 1]); line([0 1], [1 0]);

if saveg2==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['PersErrorvsCorrDay2'],'jpg');
    saveas(gcf,['PersErrorvsCorrDay2'],'fig');
end

figure(6); hold on;
title('Day3: Inbound Errors','FontSize',24,'Fontweight','bold');
xlabel('Proportion of perseverative errors','FontSize',16,'Fontweight','bold');
ylabel('Proportion of correct trajectories','FontSize',16,'Fontweight','bold');
%legend('Location','NorthEast');
axis([0 1 0 1]); line([0 1], [1 0]);

if saveg2==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['PersErrorvsCorrDay3'],'jpg');
    saveas(gcf,['PersErrorvsCorrDay3'],'fig');
end



% figure(4); hold on;
% title('Day4: Inbound Errors','FontSize',24,'Fontweight','bold');
% xlabel('Proportion of perseverative errors','FontSize',16,'Fontweight','bold');
% ylabel('Proportion of turn-around errors','FontSize',16,'Fontweight','bold');
% axis([0 1 0 1]); line([0 1], [1 0]);
% 
% figure(5); hold on;
% title('Day5: Inbound Errors','FontSize',24,'Fontweight','bold');
% xlabel('Proportion of perseverative errors','FontSize',16,'Fontweight','bold');
% ylabel('Proportion of turn-around errors','FontSize',16,'Fontweight','bold');
% axis([0 1 0 1]); line([0 1], [1 0]);






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




%% Save Data
if savedata1==1,    
    
    savefile = '/data25/sjadhav/RippleInterruption/ProcessedData/PersErrSumm_AllConAnim';
    %savefile = sprintf('%s/ProcessedData/%s_LTP_Days%01dto%01d_Tet%02d.mat', animdirect, prefix, min(days), max(days), tet);
    save(savefile);    
end






% % Day 4
% figure(4); hold on;
% orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
% plot(day4_perserr(1),day4_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');
% plot(day4_perserr(2),day4_turnerr(2),'cd','MarkerSize',10,'MarkerFaceColor','c');
% 
% % Day 5
% figure(5); hold on;
% orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
% plot(day5_perserr(1),day5_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');
% plot(day5_perserr(2),day5_turnerr(2),'cd','MarkerSize',10,'MarkerFaceColor','c');


%% EPOCH PLOTS

% figure; hold on;
% redimscreen_halfvert;
% orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
% 
% subplot(3,2,1); hold on;
% plot(epoch1_perserr(1),epoch1_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');
% plot(epoch1_perserr(2),epoch1_turnerr(2),'cd','MarkerSize',10,'MarkerFaceColor','c');
% plot(epoch1_perserr(3),epoch1_turnerr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');
% plot(epoch1_perserr(4:end),epoch1_turnerr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');
% title('Epoch1','FontSize',16,'Fontweight','bold');
% axis([0 1 0 1]); line([0 1], [1 0]);
% 
% subplot(3,2,2); hold on;
% plot(epoch2_perserr(1),epoch2_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');
% plot(epoch2_perserr(2),epoch2_turnerr(2),'cd','MarkerSize',10,'MarkerFaceColor','c');
% plot(epoch2_perserr(3),epoch2_turnerr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');
% plot(epoch2_perserr(4:end),epoch2_turnerr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');
% title('Epoch2','FontSize',16,'Fontweight','bold');
% axis([0 1 0 1]); line([0 1], [1 0]);
% 
% subplot(3,2,3); hold on;
% plot(epoch3_perserr(1),epoch3_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');
% plot(epoch3_perserr(2),epoch3_turnerr(2),'cd','MarkerSize',10,'MarkerFaceColor','c');
% plot(epoch3_perserr(3),epoch3_turnerr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');
% plot(epoch3_perserr(4:end),epoch3_turnerr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');
% title('Epoch3','FontSize',16,'Fontweight','bold');
% axis([0 1 0 1]); line([0 1], [1 0]);
% 
% subplot(3,2,4); hold on;
% plot(epoch4_perserr(1),epoch4_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');
% plot(epoch4_perserr(2),epoch4_turnerr(2),'cd','MarkerSize',10,'MarkerFaceColor','c');
% plot(epoch4_perserr(3),epoch4_turnerr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');
% plot(epoch4_perserr(4:end),epoch4_turnerr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');
% title('Epoch4','FontSize',16,'Fontweight','bold');
% axis([0 1 0 1]); line([0 1], [1 0]);
% 
% subplot(3,2,5); hold on;
% plot(epoch5_perserr(1),epoch5_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');
% plot(epoch5_perserr(2),epoch5_turnerr(2),'cd','MarkerSize',10,'MarkerFaceColor','c');
% plot(epoch5_perserr(3),epoch5_turnerr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');
% plot(epoch5_perserr(4:end),epoch5_turnerr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');
% title('Epoch5','FontSize',16,'Fontweight','bold');
% axis([0 1 0 1]); line([0 1], [1 0]);
% 
% subplot(3,2,6); hold on;
% plot(epoch6_perserr(1),epoch6_turnerr(1),'gd','MarkerSize',10,'MarkerFaceColor','g');
% plot(epoch6_perserr(2),epoch6_turnerr(2),'cd','MarkerSize',10,'MarkerFaceColor','c');
% plot(epoch6_perserr(3),epoch6_turnerr(3),'bd','MarkerSize',10,'MarkerFaceColor','b');
% plot(epoch6_perserr(4:end),epoch6_turnerr(4:end),'ko','MarkerSize',10,'MarkerFaceColor','k');
% title('Epoch6','FontSize',16,'Fontweight','bold');
% axis([0 1 0 1]); line([0 1], [1 0]);

