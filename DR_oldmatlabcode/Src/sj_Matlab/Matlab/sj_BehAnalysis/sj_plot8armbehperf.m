
saveg=0;

load /data25/sjadhav/SJStimC_direct/ProcessedData/StimC_8arm_track_perf.mat

figure(1); hold on;
%redimscreen;
orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
set(0,'defaultaxeslinewidth',2);

t=1:length(p);
plot(t,p,'r.-', 'Linewidth', 2, 'Markersize',16);

 %% Title
    title(['SJStimC 8arm track'],'FontSize',24,'Fontweight','bold');
    
    %% Axes Names
    xlabel('Day Number')
    ylabel('Performance Index')
    
axis([0.5 7.5 0.2 0.8]);
    
if saveg==1,
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['/data25/sjadhav/SJStimC_direct/Figures/SJStimC 8arm track'],'jpg');
    saveas(gcf,['/data25/sjadhav/SJStimC_direct/Figures/SJStimC 8arm track'],'fig');
    
end
