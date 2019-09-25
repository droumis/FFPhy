


%extract = 1 first then rerun for analysis with extract = 0

%to do: add state space model to this and multi day
%state space aint working yet
set(0, 'DefaultAxesFontSize',10)
%  set(0,'defaultlinelinewidth',2)
clear all; close all;
warning('off')
extract = 0;
day =15;
cohort = 'Ycoh2';
cyclefigs = 1;
savefigs = 1;
win = 20; %window for sliding average.. only use even numbers for the window... bc im lazy
saveperformance = 1;

loadtimes = sprintf('/data19/droumis/TransInf/demetrissnow19/%s_d%d.txt', cohort,day);
dir = sprintf('/data19/droumis/TransInf/%s/Day%d/%s_d%d_raw',cohort,day,cohort,day); 
outdir  = sprintf('/data19/droumis/TransInf/%s/Day%d/%s_d%d',cohort,day,cohort,day); % '/data19/droumis/TransInf/t03_28_1_out';
figdir = sprintf('/data19/droumis/TransInf/%s/',cohort);


%%
if extract == 1;
    mkdir(dir);
    mkdir(outdir);
    copyfile(sprintf('/data19/droumis/TransInf/demetrissnow19/%s_d%d.snow.dat',cohort,day), dir);
    cd(dir)
    system('nspike_extract -all ./*');
    generateTimesFileFromCont;
    diodayprocess(dir,outdir,cohort,day);
    return
end

%%
% load in the text file of [ animal start-hour start-min start-sec end-hour end-min end-sec drinking?]
%and convert it to timestamp rate [animal startstamp endstamp drinking?]
fid = fopen(loadtimes,'r+'); %open the file and give it an int handle
i = 0;
numrows = 0;
tline = fgetl(fid);
while tline ~=-1
    i = i+1;
    diostr = textscan(tline, '%d'); %this will separate the time stamp from the dio bits
    Timestamps(i,1) = diostr{1}(1);
    Timestamps(i,2) = [(diostr{1}(2)*3600 + diostr{1}(3)*60 + diostr{1}(4))*10000];
    Timestamps(i,3) = [(diostr{1}(5)*3600 + diostr{1}(6)*60 + diostr{1}(7))*10000];
    Timestamps(i,4) = diostr{1}(8);
    tline = fgetl(fid); %this will get the next line in the txt file..
end
fclose(fid)
clear fid

% for i = epochs;
if day >= 10;
load(sprintf('%s/%sDIO%d.mat',outdir,cohort,day));
else
load(sprintf('%s/%sDIO0%d.mat',outdir,cohort,day));
end
% end

%%
%Rew IN [ 1   2    3    4   5   6  H ] wells
%              [16 15 14 13 12 11 10] bits

%Rew Out [  1    2   3   4   5   6    H ] wells
%                  [ 24 23 22 21 20 19 18 ] bits

%LEDs [   1   2   3   4   5    6  H7 H8 ] wells
%           [ 32 31 30 29 28 27 26  25 ] bits

%don't change: % this assumes im using home reward out with arm1 well in and arm1 led... and then arms 2 and 3 for the outer Y
HL = 32; HW = 16; HR = 18;
AL = 31; AW = 15; AR = 23;
BL = 30; BW = 14; BR = 22;

Yperf = [];
for ii = unique(Timestamps(:,1))';
    numcorrect = 0; numhomes = 0;
    anind = find(Timestamps(:,1) == ii)';%
    for k = [anind];
        for i = [AR BR];
            numcorrect = numcorrect + length(diopulses{i}.pulsetimes(diopulses{i}.pulsetimes(:,1) > Timestamps(k,2) & diopulses{i}.pulsetimes(:,1)<Timestamps(k,3) ,1));
        end
        numhomes = numhomes + length(diopulses{HR}.pulsetimes(diopulses{HR}.pulsetimes(:,1)>Timestamps(k,2) & diopulses{HR}.pulsetimes(:,1)<Timestamps(k,3) ,1));
    end
    percentCorrect = numcorrect / numhomes;
    Yperf = [Yperf; ii numcorrect numhomes-numcorrect percentCorrect*100];
end

%% get sequential performance to feed into state space model
alldio = [];

for ii = double(unique(Timestamps(:,1))');
    %         try %this will spit out diopulses with a bit tag in the 2nd col and animal tag in the 3rd col.
    anind = find(Timestamps(:,1) == ii)';%
    for k = [anind]; % loop over each epoch for each animal. create alldio with each dio pulse, bit tag, anim tag.
        for i = 10:32;
            try
                alldio = [alldio; diopulses{i}.pulsetimes(diopulses{i}.pulsetimes(:,1) > Timestamps(k,2) & diopulses{i}.pulsetimes(:,1) < Timestamps(k,3),1) ones(length(diopulses{i}.pulsetimes(diopulses{i}.pulsetimes(:,1) > Timestamps(k,2) & diopulses{i}.pulsetimes(:,1) < Timestamps(k,3),1)),1).*i ones(length(diopulses{i}.pulsetimes(diopulses{i}.pulsetimes(:,1) > Timestamps(k,2) & diopulses{i}.pulsetimes(:,1) < Timestamps(k,3),1)),1).*ii]; %collect all stamps and tag with #bit
            end
        end
    end
end
alldiosort = sortrows(alldio,1);

performance = [];
for  j = unique(alldiosort(:,3))'; %for each animal..
    alldiosortAN = alldiosort(alldiosort(:,3) == j,:); 
    for i = 1:length(alldiosortAN(alldiosortAN(:,2) == HR & alldiosortAN(:,3) == j));%length(HomR(:,1));
        currhometime = alldiosortAN(max(find(alldiosortAN(:,2) == HR, i)),1);
        nextlightstime  = min(alldiosortAN(find(alldiosortAN(:,1)>currhometime & alldiosortAN(:,2) ~=HW,1),1));
        nextlightsbits  = alldiosortAN(find(alldiosortAN(:,1)>currhometime & alldiosortAN(:,2) ~=HW,1),2);
        nextwelltime = alldiosortAN(find(alldiosortAN(:,1)>currhometime & alldiosortAN(:,2) ==BW | alldiosortAN(:,2) ==AW,1),1);
        nextwellbit =  alldiosortAN(find(alldiosortAN(:,1)>currhometime & alldiosortAN(:,2) ==BW | alldiosortAN(:,2) ==AW,1),2);
        
        if isempty(nextwelltime);
            break
        end

        if nextlightsbits == AL;
            if nextwellbit == AW;
                performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 1 1 j day]; % pair corr anim day
            else
                performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 1 0 j day]; % pair incorr anim day
            end
        end
        
        if nextlightsbits == BL;
            if nextwellbit == BW;
                performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 2 1 j day]; % pair corr anim day
            else
                performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 2 0 j day]; % pair incorr anim day
            end
        end
        
    end
end


%%
% use 'performance' to get sliding window averages
for i = 1:length(Yperf(:,1)); %for each animal
    anperf = performance(performance(:,6) == Yperf(i,1),5); %get performance for each animal
    k = 0;
    for j = win/2:(length(anperf)-win/2);
        k = k+1;
        slide{i}(k) = mean(anperf(j-(win/2-1):j+win/2));                
    end
end

%%
% Yperf = double(Yperf);
figure; hold on;
for i = 1:length(Yperf(:,1));
subplot(3,length(Yperf(:,1)),i);
bar(1, Yperf(i,2), 'g'); hold on
bar(2, Yperf(i,3), 'b');
title(sprintf('T%d %d%%',Yperf(i,1), Yperf(i,4)));
set(gca,'XLim',[0 3]);
%  set(gca, 'Fontsize', 12);
set(gca,'XTickLabel',{'Cor','Err'})
 if i == 1;
         ylabel('#trials ')
end

subplot(3,length(Yperf(:,1)),length(Yperf(:,1))+i);
[t, p, bmode, b05, b95, bp, background_prob, sige, pmatrix, pc, lt] = DR_getestprobcorrect(performance(performance(:,6) == Yperf(i,1),5),.5,0);
        plot(t, bmode(2:end),'r-');
        hold on;
        plot(t, b05(2:end),'k', t, b95(2:end), 'k');
        hold on; [y, x] = find(bp > 0);
        h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor','k');
        set(h, 'MarkerEdgeColor', 'k');
        hold on; [y, x] = find(bp == 0);
        h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor', [0.75 0.75 0.75]);
        set(h, 'MarkerEdgeColor', 'k');
        axis([1 t(end)  0 1.05]);
        set(gca,'XTick',[0: floor(t(end)/4) :t(end)]);
        line([1 t(end)], [background_prob  background_prob ]);
%         title([{'IO(0.95)TrialLearn =' num2str(lt)}; {'RW var =' num2str(sige^2)}]);
        title(['LearnT:' num2str(lt)]);
%         set(gca, 'Fontsize', 12);
        if i == 1;
        xlabel('Trial#')
        ylabel('Prob Correct')
%         set(gca, 'Fontsize', 12);
        end
        
        %plot sliding window avg
        subplot(3, length(Yperf(:,1)),length(Yperf(:,1))*2+i);
        plot(slide{i})
        axis([1 length(slide{i})  0 1.05]);
        set(gca, 'YTick', [0:.5:1])
        line([1 length(slide{i})], [.75 .75], 'Color', [0 1 0])
        line([1 length(slide{i})], [.5 .5], 'Color', [1 0 0])
        set(gca,'XTick', [0:floor(length(slide{i})/4):length(slide{i})]);
        if i == 1;
        xlabel(sprintf('TrialWin (%dt)',win))
        ylabel('% correct')
%         set(gca, 'Fontsize', 12);
        end
%         set(gca, 'Fontsize', 12);

end
supertitle(sprintf('%s d%d',cohort,day));

figfile = [figdir 'Performance' sprintf('%s_d%d',cohort,day)];
if savefigs==1
    print('-djpeg', figfile,'-r300');
end
if cyclefigs == 1;
    keyboard;
end
close

if saveperformance == 1;
    save(sprintf('%sperformance_d%d', figdir, day), 'performance', 'slide');
end

