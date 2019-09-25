


% ***Before running this, Run extractandTrackperformance_1day for each day with saveperformance == 1
%make this agnostic to which pairs instead of needing it ordered

close all; clear all;
% animals = 'T05'; %days  11.. 12
animals = 'T08'; %days  10.. 11

% animals = 'T09T10'; %days 1:10
% animals = 'T05T08'; %days 1:10
% animals = 'T01'
% animnum = 1; %dun use this nemore
days = [11];
plotcurves = 0;
plotTTC = 1; %trials to criterion -- uses the sliding window right now.
criterion = .75 ; %
savefigs = 1;
saveperformance = 1;
cyclefigs = 1;
dir = sprintf('/data19/droumis/TransInf/%s/',lower(animals));
figdir = sprintf('/data19/droumis/TransInf/%s/',lower(animals));
win = 20; %sliding window
% pairs = [{'AB'} {'BC'} {'DE'} {'EF'} {'CD'} {'BE'} {'AF'} {'BD'} {'CE'}]; %dun touch
pairs = [{'AB'} {'BC'} {'DE'} {'EF'} {'CD'} {'BE'} {'AF'} ]; %dun touch
set(0,'defaultlinelinewidth',1);

daysstr = sprintf('%d-%d', min(days), max(days));
%% load in 'performance' mats
multdayperf = [];
for i = days;
    if strcmp('T01',animals) | strcmp('T03', animals);
        tr = 1;
        while tr
            try
                load(sprintf('%sDay%d/Performance%s_d%d_t%d.mat',dir,i,lower(animals),i,tr));
                multdayperf = [multdayperf; performance];
                tr = tr+1;
            catch
                tr = 0;
            end
        end
    else
        load(sprintf('%s/performance_d%d.mat',dir,i));
        multdayperf = [multdayperf; performance];
    end
end

% %% use 'performance' to get sliding window averages
% g = 0;
% for f = unique(multdayperf(:,7))'; %for each an
%     g = g+1;
%     for i = 1:length(unique(multdayperf(multdayperf(:,7) == f,4))); %for each pair type
%         anpairperf = multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == i,6); %get performance for each animal and each pair
%         k = 0;
%         for j = win/2:(length(anpairperf)-win/2);
%             k = k+1;
%             slidepairs{g}{i}(k) = mean(anpairperf(j-(win/2-1):j+win/2));
%         end
%     end
% end

%%
% use 'performance' to get sliding window averages
g = 0;
for f = unique(multdayperf(:,7))'; %for each an
    g = g+1;
    for pa = unique(multdayperf(multdayperf(:,7) == f,4))'; %for each pair type
        anpairperf = multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa,6); %get performance for each animal and each pair
        k = 0;
        for j = win/2:(length(anpairperf)-win/2);
            k = k+1;
            slidepairs{g}{pa}(k) = mean(anpairperf(j-(win/2-1):j+win/2));
        end
    end
end

%% TTC plot
if plotTTC == 1;
    ff = 0;
    for f = unique(multdayperf(:,7))'; %for each an
        ff = ff +1;
        figure;
        i = 0;
        for pa = unique(multdayperf(multdayperf(:,7) == f,4))'; %for each pair
            i = i+1;
            %             subplot(1,length(unique(multdayperf(multdayperf(:,7) == f,4))),i);
            hold on;
            try
            TTC = find( slidepairs{ff}{pa} >= criterion, 1);
            if isempty
            bar(pa, TTC, 'FaceColor', [i/10 i/10 i/10]);
%             text(pa-.15,-1,pairs{pa},'Color',[.5 0 0])
            end
        end
        set(gca, 'XLim', [0 8]);
        set(gca,'XTickLabel', [' ', pairs,' '])
        title(sprintf('T%d d%s ',f,daysstr));
        ylabel(sprintf('#Trials To Criterion (%d%%; Twin: %d)',criterion*100,win))
        ax1 = gca;
        hold(ax1,'all');
         ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top', 'YAxisLocation','right','Color','none','XColor','k','YColor','k');
         set(gca,'XTickLabel','')
         ylabel('%correct')
         hold(ax2, 'all');
        i = 0;
        for pa = unique(multdayperf(multdayperf(:,7) == f,4))'; %for each pair
            i = i+1;
            d = 0;
            for da = unique(multdayperf(multdayperf(:,7) == f,8))'; %for each day
                d = d +1;
            numcorr(i,d) = length(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa & multdayperf(:,8) == da &  multdayperf(:,5) == 1));
            numerrlit(i,d) = length(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa &  multdayperf(:,8) == da & multdayperf(:,5) == 2));
            numerrunlit(i,d) = length(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa &  multdayperf(:,8) == da & multdayperf(:,5) == 3));
            pairperf(i,d) = numcorr(i,d)/(numcorr(i,d)+numerrlit(i,d)+numerrunlit(i,d));
            if da == max(unique(multdayperf(multdayperf(:,7) == f,8))')
                plot(pa, pairperf(i,d)*100,'db','MarkerSize',20,'MarkerFaceColor','b');
            line([pa pa],[0 pairperf(i,d)*100], 'Color', [0 .4 .8])
            axis([0 8 0 100]);
            else
            plot(pa, pairperf(i,d)*100,'ob','MarkerSize',10,'MarkerFaceColor',[.4 1-d/10 d/10], 'MarkerEdgeColor', 'none'); %[1-sqrt(d/10) d/10 d/10],
            line([pa pa],[0 pairperf(i,d)*100], 'Color', [0 .4 .8])
            axis([0 8 0 100]);
            end
            end
        end
                    line([0 8], [75  75], 'LineStyle', '--', 'Color', [0 .6 .9])
                    line([0 8], [50 50], 'Color', 'b')
                    
                    figfile = [figdir 'MultdayTTC_' sprintf('T%d_d%s',f,daysstr)];
                    if savefigs==1
                        %             print('-dpng', figfile,'-r300');
                        print('-dpdf', figfile,'-r300');
                    end
                    if cyclefigs == 1;
                        keyboard;
                    end
                    close
    end
end

%% PLOT STUFFS
if plotcurves == 1;
    ff = 0;
    for f = unique(multdayperf(:,7))'; %for each an, create a new figure
        ff = ff +1;
        figure;
        i = 0;
        for pa = unique(multdayperf(multdayperf(:,7) == f,4))'; %for each pair
            i = i+1;
            %     for i = 1:length(unique(multdayperf(multdayperf(:,7) == f,4)));
            
            subplot(4,length(unique(multdayperf(multdayperf(:,7) == f,4))),i);
            hold on;
            %         numcorr(i) = length(multdayperf(multdayperf(multdayperf(:,7) == f,4) == i & multdayperf(multdayperf(:,7) == f,5) == 1));
            %         numerrlit(i) = length(multdayperf(multdayperf(multdayperf(:,7) == f,4) == i & multdayperf(multdayperf(:,7) == f,5) == 2));
            %         numerrunlit(i) = length(multdayperf(multdayperf(multdayperf(:,7) == f,4) == i & multdayperf(multdayperf(:,7) == f,5) == 3));
            
            numcorr(i) = length(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa &  multdayperf(:,5) == 1));
            numerrlit(i) = length(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa &  multdayperf(:,5) == 2));
            numerrunlit(i) = length(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa &  multdayperf(:,5) == 3));
            
            bar(1, numcorr(i), 'g');
            bar(2, numerrlit(i), 'k');
            bar(3, numerrunlit(i), 'FaceColor', [.5 .5 .5]);
            pairperf(i) = numcorr(i)/(numcorr(i)+numerrlit(i)+numerrunlit(i));
            title(sprintf('%s %0.0f%%', pairs{pa},100*pairperf(i)));
            set(gca,'XTick', [0:4]);
            set(gca,'XTickLabel','')
            if i == 1;
                set(gca,'XTickLabel',{'', 'Cor','Elit','E~lit', ''})
                ylabel('#Trials ')
            end
            
            subplot(4, length(unique(multdayperf(multdayperf(:,7) == f,4))),length(unique(multdayperf(multdayperf(:,7) == f,4)))+i);
            [t, p, bmode, b05, b95, bp, background_prob, sige, pmatrix, pc, lt] = DR_getestprobcorrect(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa,6),.5,0);
            line([1 t(end)], [background_prob  background_prob],  'Color', [0 0 1]);
            line([1 t(end)], [.75  .75], 'LineStyle', '--', 'Color', [0 .6 .9])
            for u = [days]; %add lines to days
                try %need this if one of the animals doesn't have data from one of the days
                    dayind = find(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa,8) == u,1);
                    line([dayind dayind], [0 1],'Color', [.7 .7 .7]);
                catch
                end
            end
            hold on;
            plot(t, bmode(2:end),'r-');
            hold on;
            plot(t, b05(2:end),'k', t, b95(2:end), 'k');
            hold on; [y, x] = find(bp > 0);
            h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor','k');
            set(h, 'MarkerEdgeColor', 'k');
            hold on; [y, x] = find(bp == 0);
            h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor', [0.75 0.75 0.75]);
            set(h, 'MarkerEdgeColor', 'k');
            axis([1 t(end)  0 1]);
            
            set(gca,'XTick',[0: floor(t(end)/4) :t(end)]);
            %         line([1 t(end)], [background_prob  background_prob ]);
            title(['LearnT:' num2str(lt)]);
            if i == 1;
                xlabel('Trial #')
                ylabel('Prob Correct')
            end
            %         for u = [days]; %add lines to days
            %             try
            %             dayind = find(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa,8) == u,1);
            %             line([dayind dayind], [0 1],'Color', [.7 .7 .7]);
            %             catch
            %             end
            %         end
            
            try
                subplot(4,length(unique(multdayperf(multdayperf(:,7) == f,4))), length(unique(multdayperf(multdayperf(:,7) == f,4)))*2+i);
                %             plot(slidepairs{ff}{i}, 'r');
                %             axis([1 length(slidepairs{ff}{i})  0 1]);
                %             set(gca, 'YTick', [0:.5:1])
                %             line([1 length(slidepairs{ff}{i})], [.75 .75], 'Color', [0 1 0])
                %             line([1 length(slidepairs{ff}{i})], [.5 .5], 'Color', [0 0 1])
                %             set(gca,'XTick', [0:floor(length(slidepairs{ff}{i})/3):length(slidepairs{ff}{i})]);
                line([1 length(slidepairs{ff}{pa})], [.75 .75], 'LineStyle', '--', 'Color', [0 .6 .9])
                line([1 length(slidepairs{ff}{pa})], [.5 .5], 'Color', [0 0 1])
                for u = [days]; %add lines to days
                    try
                        dayind = find(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa,8) == u,1);
                        line([dayind dayind], [0 1],'Color', [.7 .7 .7]);
                    catch
                    end
                end
                hold on;
                plot(slidepairs{ff}{pa}, 'r');
                axis([1 length(slidepairs{ff}{pa})  0 1]);
                set(gca, 'YTick', [0:.5:1])
                
                set(gca,'XTick', [0:floor(length(slidepairs{ff}{pa})/3):length(slidepairs{ff}{pa})]);
                
                if i == 1;
                    ylabel(sprintf('%%cor TWin:%dt',win))
                    %         ylabel('%correct')
                end
                
            catch
                plot([1:1],[1:1]);
                text(0, .5, sprintf('~enuf trials '));
            end
            
            subplot(4,length(unique(multdayperf(multdayperf(:,7) == f,4))), length(unique(multdayperf(multdayperf(:,7) == f,4)))*3+i);
            % subplot(length(unique(performance(:,4)))+2,length(unique(performance(:,4))), length(unique(performance(:,4)))+i);
            hold on;
            %         respcorr(i,1) = mean(multdayperf(multdayperf(multdayperf(:,7) == f,4) == i & multdayperf(multdayperf(:,7) == f,5) == 1,3))/10000; %put into seconds
            %         resperrlit(i,1) = mean(multdayperf(multdayperf(multdayperf(:,7) == f,4) == i & multdayperf(multdayperf(:,7) == f,5) == 2,3))/10000;
            %         resperrunlit(i,1) = mean(multdayperf(multdayperf(multdayperf(:,7) == f,4) == i & multdayperf(multdayperf(:,7) == f,5) == 3,3))/10000;
            respcorr(i,1) = mean(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa & multdayperf(:,5) == 1,3))/10000; %put into seconds
            resperrlit(i,1) = mean(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa & multdayperf(:,5) ==2,3))/10000;
            resperrunlit(i,1) = mean(multdayperf(multdayperf(:,7) == f & multdayperf(:,4) == pa & multdayperf(:,5) == 3,3))/10000;
            if isnan(respcorr(i)); respcorr(i) = 0; end
            if isnan(resperrlit(i)); resperrlit(i) = 0; end
            if isnan(resperrunlit(i)); resperrunlit(i) = 0; end
            bar(1, respcorr(i), 'g');
            bar(2, resperrlit(i), 'k');
            bar(3, resperrunlit(i), 'FaceColor', [.5 .5 .5]);
            pairperf(i) = mean([respcorr(i) resperrlit(i) resperrunlit(i)]);
            title(sprintf('mean:%0.1fs', pairperf(i)));
            set(gca,'XTick', [0:4]);
            set(gca,'XTickLabel','')
            if i == 1;
                set(gca,'XTickLabel',{ '', 'Cor','Elit','E~lit', '', ''})
                ylabel('Latency(s)')
            end
        end
        totaltrials = (sum(numcorr)+sum(numerrlit)+sum(numerrunlit));
        allperf = sum(numcorr)/totaltrials;
        allresp = mean([respcorr; resperrlit; resperrunlit]);
        %     supertitle(sprintf('Performance and ResponseTime %s d%d t%d %0.0f%% %0.2fs #%d',animals,day,trial, 100*allperf, allresp, totaltrials));
        supertitle(sprintf('T%d d%s ___ %0.0f%% ___%0.2fs ___ #%d',f,daysstr, 100*allperf, allresp, totaltrials));
        
        
        figfile = [figdir 'Multday_' sprintf('T%d_d%s',f,daysstr)];
        if savefigs==1
            print('-dpng', figfile,'-r300');
            print('-dpdf', figfile,'-r300');
        end
        if cyclefigs == 1;
            keyboard;
        end
        close
        
    end
end

if saveperformance == 1;
    save(sprintf('%smultdayperformance_d%s', figdir, daysstr), 'multdayperf');
end