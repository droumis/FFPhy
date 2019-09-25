figdir = '/data25/sjadhav/HPExpt/Figures/Ripplemod/IndividualCells';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',25);
tfont = 25;
xfont = 25;
yfont = 25;


if doindivfigs
figcnt = 0; totalplots = 0; %Count total plots across figures
nsubplots=5; % per figure
xaxis = -pret:binsize:postt; currhist = [];

for i=1:cntcells
    
if allsigshuf_sleep(i)==1 || allsigshuf_run(i)==1  % Either sleep or run
        
        curridx = allmod(i).idx;
        switch curridx(1)
            case 1
                prefix = 'HPa';
            case 2
                prefix = 'HPb';
            case 3
                prefix = 'HPc';
        end
        
        day = curridx(2); tet = curridx(3); cell = curridx(4);
        str = ''; if allsigshuf_run(i)==1, str = '*'; end
        str_sleep = ''; if allsigshuf_sleep(i)==1, str_sleep = '*'; end
        
        if mod(totalplots,5)==0 % 5 subplots finished
            figcnt=figcnt+1;
        end
        
        figure(figcnt); redimscreen; hold on;
        % Run plot
        subplot(2,5,mod(totalplots,5)+1); hold on;
        currhist = allmod(i).runhist;
        meanrate = mean(mean(currhist)); meanvec = mean(currhist); stdrate = std(meanvec);
        plot(xaxis,mean(currhist),'r','Linewidth',3);
        set(gca,'XLim',[-pret postt]);
        title(sprintf('%s D%d t%d c%d Run%s',prefix, day, tet, cell, str),'FontSize',tfont);
        % Plot vertical lines
        set(gca,'XTick',[-pret:500:postt],'XTickLabel',num2str([-pret:500:postt]'));
        ylow = min(mean(currhist)); yhigh = max(mean(currhist));
        y1 =  ylow - 0.2*meanrate;  y2 =  yhigh + 0.2*meanrate;
        set(gca,'YLim',[y1 y2]);
        ypts = y1:0.1:y2;
        % Plot Line at 0 ms - Onset of SWR
        xpts = 0*ones(size(ypts)); plot(xpts , ypts, 'b--','Linewidth',2);
        % Plot lines at rwin and bckwin
        xpts = rwin(1)*ones(size(ypts)); plot(xpts , ypts, 'b--','Linewidth',1);
        xpts = rwin(2)*ones(size(ypts)); plot(xpts , ypts, 'b--','Linewidth',1);
        xpts = bckwin(1)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
        xpts = bckwin(2)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
        % Plot horizontal Line at meanrate and n sds above/below mean
        xpts = -pret:postt; 
        ypts = meanrate*ones(size(xpts)); plot(xpts , ypts, 'k--','Linewidth',2);
        ypts = (meanrate+1*stdrate)*ones(size(xpts)); plot(xpts , ypts, 'k--','Linewidth',1);
        ypts = (meanrate-1*stdrate)*ones(size(xpts)); plot(xpts , ypts, 'k--','Linewidth',1);
        if mod(totalplots,5)==0
           ylabel('Fir rate (Hz)','FontSize',yfont,'Fontweight','normal');
           xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal'); 
        end
        
        % Sleep plot
        subplot(2,5,mod(totalplots,5)+1+5); hold on;
        currhist = allmod(i).sleephist;
        meanrate = mean(mean(currhist)); meanvec = mean(currhist); stdrate = std(meanvec);
        plot(xaxis,mean(currhist),'r','Linewidth',3);
        set(gca,'XLim',[-pret postt]);
        title(sprintf('%s D%d t%d c%d Slp%s',prefix, day, tet, cell, str_sleep),'FontSize',tfont);
        % Plot vertical lines
        set(gca,'XTick',[-pret:500:postt],'XTickLabel',num2str([-pret:500:postt]'));
        ylow = min(mean(currhist)); yhigh = max(mean(currhist));
        y1 =  ylow - 0.2*meanrate;  y2 =  yhigh + 0.2*meanrate;
        set(gca,'YLim',[y1 y2]);
        ypts = y1:0.1:y2;
        % Plot Line at 0 ms - Onset of SWR
        xpts = 0*ones(size(ypts)); plot(xpts , ypts, 'b--','Linewidth',2);
        % Plot lines at rwin and bckwin
        xpts = rwin(1)*ones(size(ypts)); plot(xpts , ypts, 'b--','Linewidth',1);
        xpts = rwin(2)*ones(size(ypts)); plot(xpts , ypts, 'b--','Linewidth',1);
        xpts = bckwin(1)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
        xpts = bckwin(2)*ones(size(ypts)); plot(xpts , ypts, 'g--','Linewidth',1);
        % Plot horizontal Line at meanrate
        xpts = -pret:postt; 
        ypts = meanrate*ones(size(xpts)); plot(xpts , ypts, 'k--','Linewidth',2);
        ypts = (meanrate+1*stdrate)*ones(size(xpts)); plot(xpts , ypts, 'k--','Linewidth',1);
        ypts = (meanrate-1*stdrate)*ones(size(xpts)); plot(xpts , ypts, 'k--','Linewidth',1);
        if mod(totalplots,5)==0
           ylabel('Fir rate (Hz)','FontSize',yfont,'Fontweight','normal');
           xlabel('Time(ms)','FontSize',xfont,'Fontweight','normal'); 
        end
            
        
        % Update plot number
        totalplots=totalplots+1;
        
        % Saving fig
        if mod(totalplots,5)==0
            figfile = [figdir,area,'_RippleAlign_IndivCells',num2str(figcnt)];
            keyboard;
            %print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
            
        end
        
        
    end % if sigshuf
end

end  % do indivfigs-- 