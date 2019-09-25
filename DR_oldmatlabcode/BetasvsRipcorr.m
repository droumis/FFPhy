
savefigs = 1;
cyclemaps =1;


% load('/home/droumis/MATLAB/allPFCCA1sigidxs_converged');
% load('/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/HP_allPFCCA1sigcorridxs_feb14_DR')
% ripcorrstruct = allPFCCA1sigcorridxs;
% betastruct = allPFCCA1sigidxs;

load('/home/droumis/MATLAB/allPFCCA1sigidxs_withcorr');
allripcorrbetas = []; sigripcorr_matchglm= []; matchripcorr_sigglm = []; sigripcorr_sigglm = [];
paircnt = 0;
for i = 1:length(allPFCCA1sigidxs);
    if ~isempty(allPFCCA1sigidxs(i).pvals);
    allripcorrbetas = [allripcorrbetas; allPFCCA1sigidxs(i).corrrvals allPFCCA1sigidxs(i).betas];
    sigripcorr_matchglm = [sigripcorr_matchglm; allPFCCA1sigidxs(i).corrrvals(allPFCCA1sigidxs(i).corrpvals<0.05) allPFCCA1sigidxs(i).betas(allPFCCA1sigidxs(i).corrpvals<0.05)];
    matchripcorr_sigglm = [matchripcorr_sigglm; allPFCCA1sigidxs(i).corrrvals(allPFCCA1sigidxs(i).pvals<0.05) allPFCCA1sigidxs(i).betas(allPFCCA1sigidxs(i).pvals<0.05)];
    sigripcorr_sigglm = [sigripcorr_sigglm; allPFCCA1sigidxs(i).corrrvals(allPFCCA1sigidxs(i).pvals <0.05 & allPFCCA1sigidxs(i).corrpvals<0.05) allPFCCA1sigidxs(i).betas(allPFCCA1sigidxs(i).pvals <0.05 & allPFCCA1sigidxs(i).corrpvals<0.05)];
    end
end



%scatter
    figure
    hold on
    scatter(allripcorrbetas(:,1), allripcorrbetas(:,2),'MarkerEdgeColor','none','MarkerFaceColor',[.8 .8 .8]);
    scatter(matchripcorr_sigglm(:,1), matchripcorr_sigglm(:,2),'MarkerEdgeColor','k','MarkerFaceColor','b');
    scatter(sigripcorr_matchglm(:,1), sigripcorr_matchglm(:,2),'MarkerEdgeColor','k','MarkerFaceColor','r');
    scatter(sigripcorr_sigglm(:,1), sigripcorr_sigglm(:,2),'MarkerEdgeColor','k','MarkerFaceColor','g');
    set(gca,'XLim',[min(allripcorrbetas(:,1)) max(allripcorrbetas(:,1))]);
  set(gca,'YLim',[min(allripcorrbetas(:,2)) max(allripcorrbetas(:,2))]);
  xlabel('RipCorr R')
  ylabel('GLMbeta')
    
  %regression
      [b,bint,r,rint,stats] = regress(allripcorrbetas(:,1), allripcorrbetas(:,2));
      [r1 p1]= corrcoef(allripcorrbetas(:,1), allripcorrbetas(:,2));
    xpts = min(allripcorrbetas(:,1)):0.01:max(allripcorrbetas(:,1));
    bfit = b(1)*xpts;
    plot(xpts,bfit,'-','Color', [.8 .8 .8],'LineWidth',2); 
    
        [b2,bint,r,rint,stats] = regress(matchripcorr_sigglm(:,1), matchripcorr_sigglm(:,2));
      [r2 p2]= corrcoef(matchripcorr_sigglm(:,1), matchripcorr_sigglm(:,2));
    bfit = b2(1)*xpts;
    plot(xpts,bfit,'-','Color', 'b','LineWidth',2);   
    
            [b3,bint,r,rint,stats] = regress(sigripcorr_matchglm(:,1), sigripcorr_matchglm(:,2));
      [r3 p3]= corrcoef(sigripcorr_matchglm(:,1), sigripcorr_matchglm(:,2));
    bfit = b3(1)*xpts;
    plot(xpts,bfit,'-','Color', 'r','LineWidth',2);
    
                [b4,bint,r,rint,stats] = regress(sigripcorr_sigglm(:,1), sigripcorr_sigglm(:,2));
      [r4 p4]= corrcoef(sigripcorr_sigglm(:,1), sigripcorr_sigglm(:,2));
    bfit = b4(1)*xpts;
    plot(xpts,bfit,'-','Color', 'g','LineWidth',2);
    
      set(gca,'XLim',[min(allripcorrbetas(:,1)) max(allripcorrbetas(:,1))]);
  set(gca,'YLim',[min(allripcorrbetas(:,2)) max(allripcorrbetas(:,2))]);
  xlabel('RipCorr R')
 ylabel('GLMbeta')
    title('(Green = SigR & Sigbeta)  (Blue = Sigbeta)   (Red = SigR)');
  
  mkdir('/mnt/data25/sjadhav/HPExpt/Figures_DR/Ripcorr_GLMbetas');
  figfile = '/mnt/data25/sjadhav/HPExpt/Figures_DR/Ripcorr_GLMbetas/ripcorr_GLMbetas_scatter';
      if ~cyclemaps == 0
        keyboard;
    end
    if savefigs==1
        print('-djpeg', figfile, '-r300');
    end
    close

% 
% 
% 
% for i = 1:length(ripcorrstruct);
%     ripcorrPFCList(i,:) = str2num(sprintf('%d',ripcorrstruct(i).PFCidx(1,:)));
% end
% 
% for ii = 1:length(betastruct);
%     betaPFCList(ii,:) = str2num(sprintf('%d',betastruct(ii).PFCidx(1,:)));
% end
% 
% paircnt = 0;
% for iii = 1:length(ripcorrstruct);
%     if ~isempty(ripcorrstruct(iii).CA1sigidxs); %if any ca1 cells
%         clear ripCA1cellsList pfcmatch
%         for j = 1:length(ripcorrstruct(iii).CA1sigidxs(:,1)); %num ca1 cells
%             ripCA1cellsList(j,:) = str2num(sprintf('%d',ripcorrstruct(iii).CA1sigidxs(j,:)));
%         end
%         pfcmatch = find(ripcorrPFCList(iii) == betaPFCList(:));
%         if ~isempty(pfcmatch) && ~isempty(betastruct(pfcmatch).pvals); %if theres a corresponding pfc cell set and with any ca1 cells
%             clear betaCA1cellsList
%             for jj = 1:length(betastruct(pfcmatch).pvals);
%                 betaCA1cellsList(jj,:) = str2num(sprintf('%d',betastruct(pfcmatch).CA1idx(jj,:)));
%             end
%             
%             for jjj = 1:length(ripCA1cellsList); %for each rip corr ca1 cell to be paired
%                 clear pairmatch
%                 pairmatch = find(ripCA1cellsList(jjj) == betaCA1cellsList(:));
%                 if ~isempty(pairmatch);
%                     paircnt  = paircnt +1;
%                     RipBetas(paircnt, [1 2]) = [ripcorrstruct(iii).rsig(jjj) betastruct(pfcmatch).betas(pairmatch)];
%                 end
%             end %for each rip corr ca1 cell
%             end %if any ca1 cells in betastruct
%         end %if any ca1 cells in ripcorr struct
%     end
%     
%     %scatter
%     figure
%     hold on
%     scatter(RipBetas(:,1), RipBetas(:,2))%,'MarkerEdgeColor','k','MarkerFaceColor',[.5 .5 .5],'LineWidth',1)
%   %regression
%       [b,bint,r,rint,stats] = regress(RipBetas(:,1), RipBetas(:,2));
%       [r2 p2]= corrcoef(RipBetas(:,1), RipBetas(:,2));
%     xpts = min(RipBetas(:,1)):0.01:max(RipBetas(:,1));
%     bfit = b(1)*xpts;
%     plot(xpts,bfit,'k-','LineWidth',1);  % Theta vs Rip
%       set(gca,'XLim',[min(RipBetas(:,1)) max(RipBetas(:,1))]);
%   set(gca,'YLim',[min(RipBetas(:,2)) max(RipBetas(:,2))]);
%   ylabel('RipCorr R')
%  xlabel('GLMbeta')
%   title({['GLMbetas x SpatCorrCoefs']; sprintf('R^2(%0.5f) p(%0.5f)', r2(2,1), p2(2,1))});
%   
%   mkdir('/mnt/data25/sjadhav/HPExpt/Figures_DR/Ripcorr_GLMbetas');
%   figfile = '/mnt/data25/sjadhav/HPExpt/Figures_DR/Ripcorr_GLMbetas/ripcorr_GLMbetas_scatter';
%     if savefigs==1
%         print('-djpeg', figfile, '-r300');
%     end
%     if ~cyclemaps == 0
%         keyboard;
%     end
%     close