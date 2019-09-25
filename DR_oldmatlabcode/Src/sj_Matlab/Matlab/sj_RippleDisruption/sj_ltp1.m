
function [ep] = sj_ltp1 (processdir,prefix,pdays,alldays,allepochs,tets,figopt1,saveg1,figopt2, saveg2,savedata)

%%% This is the old, original, pre-May2011 summary EPSP-LTP plot where
%%% I am loading one file for the given animal for each given tetrode which
%%% has all the days for that tetrode and plotting summary normalized EPSP
%%% plots for given days - individually for days and normalized across days

%%%% PLOT and compare Response to Probe Stimulation in different Sleep
%%%% Epochs across days, controlling for velocity of animal
% eg
% sj_ltp1('/data25/sjadhav/RippleInterruption/RCb_direct','RCb',1:2,1:10,[1 3 5],[3 4 9 10 11 12],1,0,0,0,0);
% sj_ltp1('/data25/sjadhav/RippleInterruption/REd_direct','REd',[1:3],[1:9],[1 3 5],[3 4 5 6 10 11],1,0,0,0,0);

% sj_ltp1('/data25/sjadhav/RippleInterruption/ProcessedData/EPSP','REd',[1:8],[1:8],[1 3 5],[2 3 4 5 6 10],1,0,1,0,0);
% sj_ltp1('/data25/sjadhav/RippleInterruption/ProcessedData/EPSP','REe',[1:7],[1:7],[1 3 5],[3 4 6 11 12 13],1,0,1,0,0);

% pdays are days you want to generate plots for
% alldays are days in file

% figopt1: Make normalized plot with Z-scores across days and tetrodes
% saveg1: Save normalized figure
% figopt2: Make normalized plots per day across all tetrodes
% saveg2: Save individual day graphs
% savedata: Save mat file with data

%%
if nargin<7,
    figopt1 = 0;
end
if nargin<8,
    saveg1 = 0;
end
if nargin<9,
    figopt2 = 0;
end
if nargin<10,
    saveg2 = 0;
end
if nargin<11,
    savedata = 0;
end
%%

pre_frac=0.5;

Fs_e = 1500; %Hz

startidx = round(pre_frac*Fs_e)+7;
endidx = round(pre_frac*Fs_e)+75;

% OR, better way - go by time
stimtime=pre_frac*Fs_e*(1000/Fs_e);
resp_win = stimtime + [5, 50]; % ms
resp_win_idx = floor(resp_win*Fs_e/1000);

set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','normal');
set(0,'defaultaxeslinewidth',2);

%% Loop Over Tetrodes

%respamp=[]; respZ=[];

for uset=1:length(tets)
    
    tet=tets(uset);
    
    % Load file
    cd(processdir);
    file = sprintf('%s/%s_EPSP_Days%01dto%01d_Tet%02d.mat', processdir, prefix, min(alldays), max(alldays), tet);
    load(file);
    
    
    for d = 1:length(pdays)
        day = pdays(d);
        
        for ep = 1:length(allepochs)
            curr_e = all_estim{day}{ep};
            for n = 1:size(curr_e,1)
                respamp{d}{ep}{uset}(n) = max(abs(curr_e(n,resp_win_idx(1):resp_win_idx(2))));
            end
        end % end epoch
        
        % For current day and tet, get baselines for epoch 1 for Z-scores
        base_mean(d,uset) = mean(respamp{d}{1}{uset});
        base_std(d,uset) = std(respamp{d}{1}{uset});
        
        % For current day, get z-scores
        for ep = 1:length(allepochs)
            curramp = respamp{d}{ep}{uset};
            respZ{d}{ep}{uset} = (curramp - repmat(base_mean(d,uset),size(curramp)))./repmat(base_std(d,uset),size(curramp));
        end % end epoch
        
        %         if figopt1==1 && uset==2
        %             figure; hold on; redimscreen_figforppt1;
        %             plot(respamp{d}{1}{uset},'ko'); plot(respamp{d}{2}{uset},'ro'); plot(respamp{d}{3}{uset},'go');
        %         end
        
        
    end  % end days
    
end % end tets


% Across Tets and Across pdays
respZep1=[]; respZep2=[]; respZep3=[];

for uset=1:length(tets)
    
    for d=1:length(pdays),
        respZep1 = [respZep1, respZ{d}{1}{uset}];
        respZep2 = [respZep1, respZ{d}{2}{uset}];
        respZep3 = [respZep1, respZ{d}{3}{uset}];
    end
    
end

%[h1,p1,ci1] = ttest2(respZep1, respZep2, 0.05);
[hk12m, pk12m, cik12m] = kstest2(respZep1, respZep2, 0.05, 'unequal');
[hk13m, pk13m, cik13m] = kstest2(respZep1, respZep3, 0.05, 'unequal');

[p12m,h12m,ci12m] = ranksum(respZep1, respZep2, 'alpha', 0.05);
[p13m,h13m,ci13m] = ranksum(respZep1, respZep3, 'alpha', 0.05);


figure; hold on;  redimscreen_figforppt1;
bar(1,mean(respZep1),'k'); errorbar(1,mean(respZep1),sem(respZep1),'k');
bar(2,mean(respZep2),'r'); errorbar(2,mean(respZep2),sem(respZep2),'r');
bar(3,mean(respZep3),'g'); errorbar(3,mean(respZep3),sem(respZep3),'g');

set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre','Post1','Post2'},'FontSize',16,'Fontweight','normal');
title('EPSP Amplitude Z-score normalized to Sleep1','FontSize',20,'Fontweight','normal');
ylabel('EPSP Amplitude Z-score','FontSize',20,'Fontweight','normal');
axis([0 4 -0.05 0.05])

if hk12m==1,
    mul = sign(mean(respZep2));
    plot(2, mean(respZep2)+1.4*mul*sem(respZep2), 'r*','MarkerSize',8);
end

if hk13m==1,
    mul = sign(mean(respZep3));
    plot(3, mean(respZep3)+1.4*mul*sem(respZep3), 'r*','MarkerSize',8);
end


%set(gca, 'XTickLabel',[]);
%set(gca, 'YTickLabel',[]);

if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,[prefix '_LTPsumm'],'fig');
    saveas(gcf,[prefix '_LTPsumm'],'jpg');
end

%%% Plotting for individual days %%%

if figopt2==1
    
    %    figure; hold on; redimscreen_figforppt1;
    %     redimscreen_2x2subplots;
    %     orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    
    
    %% Collapse Across Tetrodes
    %drespZep1=[]; drespZep2=[]; drespZep3=[];
    drespZep1 = cell(length(pdays),1);
    drespZep2 = cell(length(pdays),1);
    drespZep3 = cell(length(pdays),1);
    
    for d=1:length(pdays)
        for uset=1:length(tets)
            drespZep1{d} = [drespZep1{d}, respZ{d}{1}{uset}];
            drespZep2{d} = [drespZep2{d}, respZ{d}{2}{uset}];
            drespZep3{d} = [drespZep3{d}, respZ{d}{3}{uset}];
        end
    end
    
    %% OR
    
    %% Pick a Tetrode
    %     for d=1:length(pdays)
    %         drespZep1{d} = respZ{d}{1}{uset};  % Whatever the last uset is
    %         drespZep2{d} = respZ{d}{2}{uset};
    %         drespZep3{d} = respZ{d}{3}{uset};
    %     end
    
    
    for d=1:length(pdays)
        
        %         subplot(ceil(length(pdays)/2),2,d); hold on;
        figure; hold on; redimscreen_figforppt1;
        currZ = [mean(drespZep1{d}), mean(drespZep2{d}), mean(drespZep3{d})];
        errZ = [sem(drespZep1{d}), sem(drespZep2{d}), sem(drespZep3{d})];
        %bar([1:3],currZ,'k'); errorbar([1:3],currZ,errZ,'k.');
        
        bar(1,mean(drespZep1{d}),'k'); errorbar(1,mean(drespZep1{d}),sem(drespZep1{d}),'k');
        bar(2,mean(drespZep2{d}),'r'); errorbar(2,mean(drespZep2{d}),sem(drespZep2{d}),'r');
        bar(3,mean(drespZep3{d}),'g'); errorbar(3,mean(drespZep3{d}),sem(drespZep3{d}),'g');
        
        
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre','Post1','Post2'},'FontSize',16,'Fontweight','normal');
        title(['Day' num2str(pdays(d)) ' : EPSP Amplitude Z-score norm to Sleep1'],'FontSize',20,'Fontweight','normal');
        ylabel('EPSP Amplitude Z-score','FontSize',20,'Fontweight','normal');
        axis([0 4 -0.5 0.5])
        
        %         if d==1,
        %             title([prefix ' Tet' num2str(tet) ' : Day ' num2str(pdays(d)) ],'FontSize',16,'Fontweight','normal');
        %             ylabel(['Z-scores']);
        %         else
        %             title(['Day ' num2str(pdays(d)) ],'FontSize',16,'Fontweight','normal');
        %         end
        
        [hk12(d), pk12(d), cik12{d}] = kstest2(drespZep1{d}, drespZep2{d}, 0.05, 'unequal');
        [hk13(d), pk13(d), cik13{d}] = kstest2(drespZep1{d}, drespZep3{d}, 0.05, 'unequal');
        
        [p12(d),h12(d),ci12{d}] = ranksum(drespZep1{d}, drespZep2{d}, 'alpha', 0.05);
        [p13(d),h13(d),ci13{d}] = ranksum(drespZep1{d}, drespZep3{d}, 'alpha', 0.05);
        
        if hk12(d)==1,
            mul = sign(currZ(2));
            plot(2, (currZ(2)+1.4*mul*errZ(2)), 'r*','MarkerSize',8);
        end
        
        if hk13(d)==1,
            mul = sign(currZ(3));
            plot(3, (currZ(3)+1.4*mul*errZ(3)), 'r*','MarkerSize',8);
        end
        
        if saveg1==1,
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
            saveas(gcf,['Day' num2str(pdays(d)) '_LTPsumm'],'fig');
            saveas(gcf,['Day' num2str(pdays(d)) '_LTPsumm'],'jpg');
        end
        
    end
    
end



if savedata==1,
    savefile = sprintf('/data25/sjadhav/RippleInterruption/ProcessedData/VelAtStimln_Compare.mat');
    savefile = sprintf('%s/%s_Tet%02d_LTPsumm.mat', processdir, prefix, tet);
    save(savefile);
end

keyboard;








