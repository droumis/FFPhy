
function [ep] = sj_stim_behstats_velcompare (alldatadir,prefixes,figopt1,saveg1, savedata)

%% Compare Velocity at Stimulations across animals, and also Stimln Vel vs
%% All Vel
% eg
% sj_stim_behstats_velcompare('/data25/sjadhav/RippleInterruption/ProcessedData',{'sjc','RE1'},1,0,0);

%%
if nargin<3,
    figopt1 = 1;
end
if nargin<4,
    saveg1 = 0;
end
if nargin<5,
    savedata = 0;
end


% Variable parameters
thrsvel=5;  %<5cm per sec is still on track (3 for sleep session in sj_stimresp2_withvel.m)
dur = 15; % Standard epoch duration in mins

%% Initialize
directoryname = alldatadir;
if (directoryname(end) == '/')
    directoryname = directoryname(1:end-1);
end
cd(directoryname);
clr = {'b','r','g','c','m','y','k','b','r','g','c'};

%% Loop over animals
store_histv=[]; store_hv=[];
for a = 1:length(prefixes)
    
    currprefix = prefixes{a};
    if currprefix=='sjc',
        allepochs = 2;
    else
        allepochs = [2,4];
    end
    
    histv = []; hv = [];  norm_histv = []; norm_hv = [];
    directoryname = alldatadir;
    file = sprintf('%s/%s_StimlnBehStats.mat', directoryname, currprefix);
    load(file);
    
    %for d=1:length(day_velstim)
    for d=1:4
        
        for ep = 1:length(allepochs)
            
            currvelstim = day_velstim{d}{ep};
            histv = [histv; hist(currvelstim,[0:1:50])];
            norm_histv = [norm_histv; hist(currvelstim,[0:1:50])./max(hist(currvelstim,[0:1:50]))];
            
            velo = day_vel{d}{ep};
            hv = [hv; hist(velo,[0:1:50])];
            norm_hv = [norm_hv; hist(velo,[0:1:50])./max(hist(velo,[0:1:50]))];
            
        end
    end
    
    meanv(a,:) = mean(histv); norm_meanv(a,:) = mean(norm_histv);
    err1v(a,:) = std(histv);  norm_err1v(a,:) = std(norm_histv);
    err2v(a,:) = sem(histv);  norm_err2v(a,:) = sem(norm_histv);
    
    meanhv(a,:) = mean(hv); norm_meanhv(a,:) = mean(norm_hv);
    err1hv(a,:) = std(hv); norm_err1hv(a,:) = std(norm_hv);
    err2hv(a,:) = sem(hv); norm_err2hv(a,:) = sem(norm_hv);
    
    store_histv{a} = histv;
    store_hv{a} = hv;
    
end


%% Figure;
figopt1=1;
if figopt1==1
    % Plot across animal
    figure; hold on;
    redimscreen_halfvert
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
    set(0,'defaultaxeslinewidth',2);
    
    for a = 1:length(prefixes)
        plot([0:1:50],meanv(a,:),[clr{a} '.-'],'MarkerSize',8,'LineWidth',2);
    end
    
    xlabel('Velocity (cm/sec)');
    ylabel(['No of Stimulations']);
    title(['Vel Distributions at Stimuln: Days 1' ' to ' num2str(length(day_velstim)-1) ],'FontSize',20,'Fontweight','bold');
    legend(prefixes);
    
    for a = 1:length(prefixes)
        uperr = meanv(a,:) + err2v(a,:);
        lowerr = meanv(a,:) - err2v(a,:);
        jbfill([0:1:50],uperr,lowerr,clr{a},clr{a},1,0.2);
    end
    
    saveg1=1;
    if saveg1==1,
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['VelAtStimln_Compare'],'fig');
        saveas(gcf,['VelAtStimln_Compare'],'jpg');
    end
    
    %%%% Plot for within animal
    
    for a = 1:length(prefixes)
        
        currprefix = prefixes{a};
        figure; hold on;
        redimscreen_halfvert
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
        set(0,'defaultaxeslinewidth',2);
                
        plot([0:1:50],norm_meanv(a,:),['k.-'],'MarkerSize',8,'LineWidth',2);
        plot([0:1:50],norm_meanhv(a,:),['g.-'],'MarkerSize',8,'LineWidth',2);
        
        xlabel('Velocity (cm/sec)');
        ylabel(['Normalized no. of occurences']);
        title([currprefix ': Vel Distributions at Stimuln: Days 1 to 4'],'FontSize',20,'Fontweight','bold');
        legend('Stimln Vel','All Vel');
        
        uperr = norm_meanv(a,:) + norm_err2v(a,:); lowerr = norm_meanv(a,:) - norm_err2v(a,:);
        jbfill([0:1:50],uperr,lowerr,'k','k',1,0.2);
        
        uperr = norm_meanhv(a,:) + norm_err2hv(a,:); lowerr = norm_meanhv(a,:) - norm_err2hv(a,:);
        jbfill([0:1:50],uperr,lowerr,'g','g',1,0.2);
        
        if saveg1==1,
            orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
            saveas(gcf,[currprefix 'VelAtStimln'],'fig');
            saveas(gcf,[currprefix 'VelAtStimln'],'jpg');
        end
        
        
    end
    
    
end



savedata=1;
if savedata==1,
    savefile = sprintf('/data25/sjadhav/RippleInterruption/ProcessedData/VelAtStimln_Compare.mat');
    save(savefile);
end


