
function [ep] = sj_ltp1_summ1 (processdir,prefixes,pdays,allepochs,figopt1,saveg1,figopt2, saveg2,savedata)

%% Adding across animal summary to sj_ltp1, which uses files for each
%% tetrode across days (so indiv tet files must include all days used)

%%% This is the old, original, pre-May2011 summary EPSP-LTP plot where
%%% I am loading one file for the given animal for each given tetrode which
%%% has all the days for that tetrode and plotting summary normalized EPSP
%%% plots for given days - individually for days and normalized across days

%%%% PLOT and compare Response to Probe Stimulation in different Sleep
%%%% Epochs across days, controlling for velocity of animal
% eg

%  sj_ltp1_summ1('/data25/sjadhav/RippleInterruption/ProcessedData/EPSP',{'REd';'REe';'REf'},[1:8],[1 3 5],1,0,1,0,0);

% sj_ltp1('/data25/sjadhav/RippleInterruption/RCb_direct','RCb',1:2,1:10,[1 3 5],[3 4 9 10 11 12],1,0,0,0,0);
% sj_ltp1('/data25/sjadhav/RippleInterruption/REd_direct','REd',[1:3],[1:9]
% ,[1 3 5],[3 4 5 6 10 11],1,0,0,0,0);
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
if nargin<5,
    figopt1 = 0;
end
if nargin<6,
    saveg1 = 0;
end
if nargin<7,
    figopt2 = 0;
end
if nargin<8,
    saveg2 = 0;
end
if nargin<9,
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

novel =1:3;
fam = 4:8;

set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','normal');
set(0,'defaultaxeslinewidth',2);

%------------------------------------------------------

cd(processdir);

% Initialize arrays which will hold Z-scores across tetrodes across
% animals, for each Day individually, and also across days
respZep1=[]; respZep2=[]; respZep3=[]; % Across tetrodes, across days, across animals
drespZep1 = cell(length(pdays),1); % For each day, across tetrodes and animals
drespZep2 = cell(length(pdays),1);
drespZep3 = cell(length(pdays),1);

% For novel and fimiliar days
respZep1n=[]; respZep2n=[]; respZep3n=[];
respZep1f=[]; respZep2f=[]; respZep3f=[];

for anim=1:length(prefixes)
    
    % Initalize for each animal %
    uset=[]; days=[]; respZ=[];
    
    prefix=prefixes{anim}
    files=dir([prefix '*']);
    
    % Loop over files which have data for 1 tetrode across all days and get
    % terode list and day range for which data is available for current
    % animal
    
    for f=1:length(files),
        [name,ext]=strtok(files(f).name,'.');
        tets(f) = str2num(name(end-1:end));
        
        % Old method - get days from name
%          if f==1,
%             [part,part2]=strtok(name,'s');
%             day_a=str2num(part2(2));
%             day_b=str2num(part2(5));
%             days=[day_a:day_b];
%         end
        
        % New Method - get days from stored variable in file in loop below
        

    end
    
    for f=1:length(files), % Loop over each file aka each tetrode
        
        uset = tets(f);
        files(f).name;
        load(files(f).name);
        
        % New Method - get days from stored variable in file
        
        if isempty(days)
            if(strcmp(prefix,'REd'))
                days=1:8;
            end
            if(strcmp(prefix,'REe'))
                days=1:7;
            end
        end
        
        for d = 1:length(days)
            day = days(d);
            
            % Initialize
            respamp=[];
            
            for ep = 1:length(allepochs)
                curr_e = all_estim{day}{ep};
                for n = 1:size(curr_e,1)
                    %respamp{d}{ep}{uset}(n) = max(abs(curr_e(n,resp_win_idx(1):resp_win_idx(2))));
                    respamp{ep}(n) = max(abs(curr_e(n,resp_win_idx(1):resp_win_idx(2))));
                end
            end % end epoch
            
            % For current day and tet, get baselines for epoch 1 for Z-scores
            %base_mean(d,uset) = mean(respamp{d}{1}{uset});
            %base_std(d,uset) = std(respamp{d}{1}{uset});
            base_mean = mean(respamp{1});
            base_std = std(respamp{1});
            
            % For current day, get z-scores
            for ep = 1:length(allepochs)
                curramp = respamp{ep};
                respZ{ep} = (curramp - repmat(base_mean,size(curramp)))./repmat(base_std,size(curramp));
            end % end epoch
            
            %         if figopt1==1 && uset==2
            %             figure; hold on; redimscreen_figforppt1;
            %             plot(respamp{1}{uset},'ko'); plot(respamp{2}{uset},'ro'); plot(respamp{3}{uset},'go');
            %         end
            
            
            % Save for current day and current tetrode - already in loops
            respZep1 = [respZep1, respZ{1}];
            respZep2 = [respZep1, respZ{2}];
            respZep3 = [respZep1, respZ{3}];
            drespZep1{day} = [drespZep1{day}, respZ{1}];
            drespZep2{day} = [drespZep2{day}, respZ{2}];
            drespZep3{day} = [drespZep3{day}, respZ{3}];
            
            if ~isempty(intersect(day,novel))
                respZep1n = [respZep1n, respZ{1}];
                respZep2n = [respZep1n, respZ{2}];
                respZep3n = [respZep1n, respZ{3}];
            end
            
            if ~isempty(intersect(day,fam))
                respZep1f = [respZep1f, respZ{1}];
                respZep2f = [respZep1f, respZ{2}];
                respZep3f = [respZep1f, respZ{3}];
            end
            
        end  % end days
        
    end % end files - loop across tetrodes for current animal
    
    
end  % End animal


%%% Master plot

[ht12m,pt12m,cit12m] = ttest2(respZep1, respZep2, 0.05);
[ht13m,pt13m,cit13m] = ttest2(respZep1, respZep3, 0.05);

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
axis([0 4 -0.08 0.08])

if ht12m==1,
    mul = sign(mean(respZep2));
    plot(2, mean(respZep2)+1.4*mul*sem(respZep2), 'r*','MarkerSize',8);
end

if ht13m==1,
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


%%% Novel Days Plot

[ht12n,pt12n,cit12n] = ttest2(respZep1n, respZep2n, 0.05);
[ht13n,pt13n,cit13n] = ttest2(respZep1n, respZep3n, 0.05);


figure; hold on;  redimscreen_figforppt1;
bar(1,mean(respZep1n),'k'); errorbar(1,mean(respZep1n),sem(respZep1n),'k');
bar(2,mean(respZep2n),'r'); errorbar(2,mean(respZep2n),sem(respZep2n),'r');
bar(3,mean(respZep3n),'g'); errorbar(3,mean(respZep3n),sem(respZep3n),'g');

set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre','Post1','Post2'},'FontSize',16,'Fontweight','normal');
title('Novel Days: EPSP Amplitude Z-score normalized to Sleep1','FontSize',20,'Fontweight','normal');
ylabel('EPSP Amplitude Z-score','FontSize',20,'Fontweight','normal');
axis([0 4 -0.08 0.08])

if ht12n==1,
    mul = sign(mean(respZep2n));
    plot(2, mean(respZep2n)+1.4*mul*sem(respZep2n), 'r*','MarkerSize',8);
end

if ht13n==1,
    mul = sign(mean(respZep3n));
    plot(3, mean(respZep3n)+1.4*mul*sem(respZep3n), 'r*','MarkerSize',8);
end

%set(gca, 'XTickLabel',[]);
%set(gca, 'YTickLabel',[]);

if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,[prefix '_NovelLTPsumm'],'fig');
    saveas(gcf,[prefix '_NovelLTPsumm'],'jpg');
end


%%% Familiar Days Plot

[ht12f,pt12f,cit12f] = ttest2(respZep1f, respZep2f, 0.05);
[ht13f,pt13f,cit13f] = ttest2(respZep1f, respZep3f, 0.05);


figure; hold on;  redimscreen_figforppt1;
bar(1,mean(respZep1f),'k'); errorbar(1,mean(respZep1f),sem(respZep1f),'k');
bar(2,mean(respZep2f),'r'); errorbar(2,mean(respZep2f),sem(respZep2f),'r');
bar(3,mean(respZep3f),'g'); errorbar(3,mean(respZep3f),sem(respZep3f),'g');

set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre','Post1','Post2'},'FontSize',16,'Fontweight','normal');
title('Fam Days: EPSP Amplitude Z-score normalized to Sleep1','FontSize',20,'Fontweight','normal');
ylabel('EPSP Amplitude Z-score','FontSize',20,'Fontweight','normal');
axis([0 4 -0.08 0.08])

if ht12f==1,
    mul = sign(mean(respZep2f));
    plot(2, mean(respZep2f)+1.4*mul*sem(respZep2f), 'r*','MarkerSize',8);
end

if ht13f==1,
    mul = sign(mean(respZep3f));
    plot(3, mean(respZep3f)+1.4*mul*sem(respZep3f), 'r*','MarkerSize',8);
end

%set(gca, 'XTickLabel',[]);
%set(gca, 'YTickLabel',[]);

if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,[prefix '_FamLTPsumm'],'fig');
    saveas(gcf,[prefix '_FamlLTPsumm'],'jpg');
end



%%% Plotting for individual days %%%

if figopt2==1
    
    
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
        
        [ht12(d), pt12(d), cit12{d}] = ttest2(drespZep1{d}, drespZep2{d}, 0.05);
        [ht13(d), pt13(d), cit13{d}] = ttest2(drespZep1{d}, drespZep3{d}, 0.05);
        
        [hk12(d), pk12(d), cik12{d}] = kstest2(drespZep1{d}, drespZep2{d}, 0.05, 'unequal');
        [hk13(d), pk13(d), cik13{d}] = kstest2(drespZep1{d}, drespZep3{d}, 0.05, 'unequal');
        
        [p12(d),h12(d),ci12{d}] = ranksum(drespZep1{d}, drespZep2{d}, 'alpha', 0.05);
        [p13(d),h13(d),ci13{d}] = ranksum(drespZep1{d}, drespZep3{d}, 'alpha', 0.05);
        
        if ht12(d)==1,
            mul = sign(currZ(2));
            plot(2, (currZ(2)+1.4*mul*errZ(2)), 'r*','MarkerSize',8);
        end
        
        if ht13(d)==1,
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

%keyboard;








