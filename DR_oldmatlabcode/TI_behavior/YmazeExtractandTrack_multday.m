

%this is weird for T1 and 3 because i split up the epochs into multiple saved files... i won't be doing this for the next cohort
%right now if i run it for these animals it only uses 1 of the epochs per day..

%Before running this, Run extractandTrackperformance_1day for each day with saveperformance == 1
close all; clear all;
set(0, 'DefaultAxesFontSize',10)
 set(0,'defaultlinelinewidth',1)
animals = 'Ycoh2';
animnum = 1;
days = [9:15];
daysstr = '9-15';
win = 20;
savefigs = 1;
saveperformance = 1;
cyclefigs = 1;
dir = sprintf('/data19/droumis/TransInf/%s/',animals);
figdir = sprintf('/data19/droumis/TransInf/%s/',animals);

% multdayslide = {};
multdayperf = [];
for i = 1:length(days);
    load(sprintf('%sperformance_d%d.mat',dir,days(i)));
    multdayperf = [multdayperf; performance];
    % if i == 1;
    %     multdayslide = slide;
    % else
    %     for ii = 1:length(slide);
    %         multdayslide{ii} = [multdayslide{ii} slide{ii}];
    %     end
    % end
end
anims = [unique(multdayperf(:,6))];

% clear slide
%%
% use 'performance' to get sliding window averages
for i = 1: length(anims); %for each an
    anperf = multdayperf(multdayperf(:,6) == anims(i,1),5); %get performance for each animal
    k = 0;
    for j = win/2:(length(anperf)-win/2);
        k = k+1;
        multdayslide{i}(k) = mean(anperf(j-(win/2-1):j+win/2));
    end
end

%%
anims = [unique(multdayperf(:,6))];
for ii = 1: length(anims); %for each an
    subplot(2,length(unique(multdayperf(:,6))),ii);
    [t, p, bmode, b05, b95, bp, background_prob, sige, pmatrix, pc, lt] = DR_getestprobcorrect(multdayperf(multdayperf(:,6) == anims(ii),5),.5,0);
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
    line([1 t(end)], [background_prob  background_prob], 'Color', [0 0 1]);
    line([1 t(end)], [.7  .7], 'Color', [0 1 0])
    %         title([{'IO(0.95)TrialLearn =' num2str(lt)}; {'RW var =' num2str(sige^2)}]);
    title([{sprintf('T%d',anims(ii))}; {['LearnT:' num2str(lt)]}]);
    set(gca,'XTick', [0:floor(length(t)/3):length(t)]);
    for i = [days]; %add lines to days
        dayind = find(multdayperf(multdayperf(:,6) == anims(ii),7) == i,1);
        line([dayind dayind], [0 1],'Color', [.7 .7 .7]);
    end
    if ii == 1;
        xlabel('Trial#')
        ylabel('Prob Correct')
    end
    
    subplot(2,length(unique(multdayperf(:,6))),length(unique(multdayperf(:,6)))+ii);
    plot(multdayslide{ii}, 'r');
    axis([1 length(multdayslide{ii})  0 1]);
    set(gca, 'YTick', [0:.5:1])
    line([1 length(multdayslide{ii})], [.75 .75], 'Color', [0 1 0])
    line([1 length(multdayslide{ii})], [.5 .5], 'Color', [0 0 1])
    set(gca,'XTick', [0:floor(length(multdayslide{ii})/3):length(multdayslide{ii})]);
    if ii == 1;
        xlabel(sprintf('TrialWin (%dt)',win))
        ylabel('% correct')
    end
    for i = [days]; %add lines to days
        dayind = find(multdayperf(multdayperf(:,6) == anims(ii),7) == i,1) - win/2;
        line([dayind dayind], [0 1],'Color', [.7 .7 .7]);
    end
end
supertitle(sprintf('%s d%s',animals,daysstr));

figfile = [figdir 'MultdayPerformance' sprintf('%s_d%s',animals,daysstr)];
if savefigs==1
    print('-djpeg', figfile,'-r300');
    print('-dpdf', figfile,'-r300');
end
if cyclefigs == 1;
    keyboard;
end
close

if saveperformance == 1;
    save(sprintf('%sMultdayPerformance_d%s', figdir, daysstr), 'multdayperf', 'multdayslide');
end