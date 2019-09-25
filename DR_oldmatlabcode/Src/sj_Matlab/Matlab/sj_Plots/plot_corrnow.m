%-------------------------------------------------------
%


% ------------------------------
% Figure and Font Sizes

figdir = '/data25/sjadhav/HPExpt/HPa_direct/Figures/';
forppr = 0;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

summdir = figdir;
%set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',12);
    tfont = 12; % title font
    xfont = 12;
    yfont = 12;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end

clr = {'b','r','g','c','m','y','k','r'};

% ---------------------------------------

prefix='HPa';

% Ripple Corr
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end

plot(sleepcorrtime, plotcorr_sl_ep1,'k','LineWidth',3);
%plot(sleepcorrtime, plotcorr_sl_ep3,'b','LineWidth',3);
plot(sleepcorrtime, plotcorr_sl_ep5,'r','LineWidth',3);
plot(sleepcorrtime, plotcorr_runsl,'r--','LineWidth',3);
line([0 0], [0 max([plotcorr_sl_ep1, plotcorr_runsl, plotcorr_sl_ep5])],'Color',[0.5 0.5 0.5],'LineWidth',2);
xlabel('Time (sec)');
ylabel('Prob-Sleep Corrln in ripples')
set(gca,'XLim',[-0.1 0.1]);
set(gca,'YLim',[0 max([plotcorr_sl_ep1, plotcorr_runsl, plotcorr_sl_ep5])+0.001]);

title(['Nev1: ',num2str(Neventscorr_sl_ep1),'; Nevrr: ',num2str(Neventscorr_runsleps),'; Nev3: ',num2str(Neventscorr_sl_ep5)]...
    ,'FontSize',20,'Fontweight','normal');

saveg=1;
if saveg==1
    figfile = [figdir,prefix,'CA1t4c2_PFCt18c2_RipCorrln'];
    print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end


% Run Corr
% ----------

% Ripple Corr
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end

figdir = '/data25/sjadhav/HPExpt/HPa_direct/Figures/RunCorr/'
plot(runcorrtime, plotcorr_runpl,'b','LineWidth',3);
line([0 0], [0 max(plotcorr_runpl)],'Color',[0.5 0.5 0.5],'LineWidth',2);
xlabel('Time (sec)');
ylabel('Prob-Run Corrln')
set(gca,'XLim',[-0.31 0.31]);
set(gca,'YLim',[0 max(plotcorr_runpl)+0.001]);
title(['Nevspl: ',num2str(Neventscorr_runpleps)],'FontSize',20,'Fontweight','normal');

if saveg==1
    figfile = [figdir,prefix,'CA1t4c2_PFCt18c2_RunCorrln'];
    print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end














% cmap = jet(1024);
% cmap(1,:) = 1;
% colormap(cmap);
% % set up the bounds to make it look good
% bounds = [0 0];
% if (isempty(peakrate))
%     peakrate = max(rmap(:));
%     bounds(2) = peakrate * 0.65; %0.65;
% else
%     bounds(2) = peakrate;
% end
% if (peakrate > 20)
%     minbound = -1;
% elseif (peakrate > 10)
%     minbound = -.5;
%     rmap(find(rmap == -1)) = -.5;
% elseif (peakrate > 3)
%     minbound = -.1;
%     rmap(find(rmap == -1)) = -.1;
% else
%     minbound = -.01;
%     rmap(find(rmap == -1)) = -.01;
% end
% 
% bounds(1) = minbound;
% if (peakrate < .1)
%     bounds(2) = 1;
% end
% % set up the colormap to be white at some negative values
% if (peakrate >= minrate)
%     h = imagesc(flipud(rmap), bounds);
%     ch = colorbar;
%     set(ch, 'FontSize', fontsize);
%     if (showmax)
%         set(ch, 'YTick', floor(peakrate * 0.65), 'YTickLabel', num2str(floor(peakrate * 1))); 
%     end
%     %title(['Peakrate: ',num2str(floor(peakrate * 1))]);
% end
% 
% axis off;
% axis equal

