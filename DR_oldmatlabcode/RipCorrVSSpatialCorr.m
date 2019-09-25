
cyclemaps = 1;
savefigs = 1;
fonttype = 'Arial';
titlesize = 16;


load('/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/unmodRs')
load('/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/modRs')


unmodRs_signnonsig = load('/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/unmodRs_signonsig')
unmodRs_signnonsig = unmodRs_signnonsig.unmodRs;

modRs_signnonsig = load('/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/modRs_signonsig')
modRs_signnonsig = modRs_signnonsig.modRs;

allRs = [modRs_signnonsig; unmodRs_signnonsig];
allmodRs = [modRs; unmodRs];

for i=1:length(allmodRs(:,1));
    clear tempexcl
    tempexcl = rowfind(allmodRs(i,:),allRs(:,:));
    if tempexcl ~= 0;
        allRs(tempexcl,:) = [];
    end
end
    

figure
hold on;

scatter(unmodRs_signnonsig(:,1), unmodRs_signnonsig(:,2),10,'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7]);
scatter(modRs_signnonsig(:,1), modRs_signnonsig(:,2),10,'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7]);

scatter(modRs(modRs(:,1)>0,1), modRs(modRs(:,1)>0,2), 50, 'MarkerEdgeColor','k','MarkerFaceColor','b');
scatter(modRs(modRs(:,1)<0,1), modRs(modRs(:,1)<0,2), 50, 'MarkerEdgeColor','k','MarkerFaceColor',[.8 .7 0]);


    set(gca,'XLim',[min(allRs(:,1)) max(allRs(:,1))]);
  set(gca,'YLim',[min(allRs(:,2)) max(allRs(:,2))]);
  xlabel('Ripple CorrCoef')
  ylabel('Spatial CorrCoef')


hold on;
[b, bint, resi, rint, stats] = regress(allmodRs(:,2), [allmodRs(:,1), ones(size(allmodRs(:,1)))]);
r = polyval(b,allmodRs(:,1));
plot(allmodRs(:,1), r, '-g')

[b2, bint, resi, rint, stats2] = regress(allRs(:,2), [allRs(:,1), ones(size(allRs(:,1)))]);
r2 = polyval(b2,allRs(:,1));
plot(allRs(:,1), r2, '-k')


  
title({['mPFC-CA1 Pairs']; [sprintf('Green = ripmodPFC/sigripcorrPairs R2(%0.5f) p(%0.5f)', stats(1), stats(3))]; [sprintf('Black = ALL R2(%0.5f) p(%0.5f)', stats2(1), stats2(3))] });
set(gca, 'Fontsize', titlesize,'FontName',fonttype)
          set(gcf,'PaperPositionMode','auto');
        set(gcf, 'Position', [100 100 500 600])

  mkdir('/mnt/data25/sjadhav/HPExpt/Figures_DR/Ripcorr_SpatCorr');
  figfile = '/mnt/data25/sjadhav/HPExpt/Figures_DR/Ripcorr_SpatCorr/ripcorr_spatcorr_excluded';
      if ~cyclemaps == 0
        keyboard;
    end
    if savefigs==1
        print('-djpeg', figfile, '-r300');
        print('-dpdf', figfile, '-r300');
    end
    close