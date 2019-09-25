
%add varg for last day plotting 
%add lines to demarcate 1pair epochs and weighted epochs. 

%6/23/14 I messed up with the delays for big day t05 and set them to 1500 instead of 15000... so i had to edit the code to backup a second to look for the next lights bits... this should still be compatible with other animals/days..

%haven't finished getting t01/3 to work.. looks weird. 
%have i got to colect all t01/t03 trials yet?? i think maybe...?
%produce upated figures for these guys

warning('off')
clear all; close all;
extract = 0;
% extract = 0;
day = 17;

% animals = 'T05T08'; %this is the cohort of 2 animals
% animals = 'T08'; 
% animals = 'T05';
% spatconfig = 3;
animals = 'T09T10';
% animals = 'T01';
spatconfig = 4;

% phase = 2; %NOT USING ANYMORE
savefigs =1;
cyclefigs = 1;
win = 10; %sliding window
% statespace = 0;


usetimestxt = 1; %load in a txt file of the times and animal(s). it's this OR (exclusive) using timewise/useallpulses

timewise = 0; %to specify a certain timerange to look in.. it's either this OR (exclusive) loading a times file.. I probably only want to use this for the last day where i have multiple different kinds of epochs in one saved files and i want to look at a specific one.
useallpulses = 0; %if i don't know the start and end time and i just want to look at all recorded dio pulses from a specific animal.. use this OR (ex) timewise
animnum = 5; %just specify this for now if using timewise... T01 ~ 1. T08 ~8. etc
trial = 1; %use a single trial for extraction NOT USING THIS ANTMORE
starttime = [1 57 54]; %only if timewise = 1
endtime = [2 16 20]; %only if timewise = 1

% epochs = [1 2 3 4]; % use multiple epochs for analysis... don't use this
loadtimes = sprintf('/data19/droumis/TransInf/demetrissnow19/%s/%s_d%d.txt',animals,animals,day);
% dir = sprintf('/data19/droumis/TransInf/%s/Day%d/%s_d%d_t%d_raw',lower(animals),day,lower(animals),day,trial);
% outdir  = sprintf('/data19/droumis/TransInf/%s/Day%d/%s_d%d_t%d',lower(animals),day,lower(animals),day,trial); ;% '/data19/droumis/TransInf/t03_28_1_out';
dir = sprintf('/data19/droumis/TransInf/%s/Day%d/%s_d%d_raw',lower(animals),day,lower(animals),day);
outdir  = sprintf('/data19/droumis/TransInf/%s/Day%d/%s_d%d',lower(animals),day,lower(animals),day);% '/data19/droumis/TransInf/t03_28_1_out';
figdir = sprintf('/data19/droumis/TransInf/%s/',lower(animals));
saveperformance = 1;
set(0, 'DefaultAxesFontSize',8)



%which track is the corresponding letter on.. this is all you need to change between different animal track configurations
%%%%%%%%%%A  B  C  D  E  F
if spatconfig == 2;
    configsetup =  [  6  2   4  5   3  1 ] ; %ANIMALS T01 and T03
elseif spatconfig == 3;
    configsetup =  [  4  2   6  1   5  3 ] ; %ANIMALS T05 and T08
elseif spatconfig == 4;
    configsetup =  [  4  6   1  3   2  5 ] ; %ANIMALS T09 and T10
end

%% extract

if extract == 1;
    mkdir(dir);
    mkdir(outdir);
    %     copyfile(sprintf('/data19/droumis/TransInf/demetrissnow19/%s_ph%d_d%d_t%d.snow.dat',animals,phase, day, trial), dir);
    copyfile(sprintf('/data19/droumis/TransInf/demetrissnow19/%s/%s_d%d.snow.dat',animals, animals, day), dir);
    cd(dir)
    system('nspike_extract -all ./*');
    generateTimesFileFromCont;
    diodayprocess(dir,outdir,animals,day);
    return
end


%% use specified time range for a single animal

%%convert to time stamps
if timewise == 1;
    startstamp = [(starttime(1)*3600+starttime(2)*60+starttime(3))*10000];
    endstamp = [(endtime(1)*3600+endtime(2)*60+endtime(3))*10000];
    load(sprintf('%s/%sDIO%d.mat',outdir,animals,day));
elseif useallpulses == 1;
    if animals == 'T01' | animals == 'T03'; %for t01 and t03, loop load the trials within a day... eeeugghh why did i do this to myself.
        tr = 1;
        while tr
            try
                load(sprintf('%s_t%d/%sDIO%d.mat',outdir,tr,animals,day));
                if tr == 1;
                    diopulsesTMP = diopulses;
                else
                    for dp = 1: length(diopulses);
                        try
                            diopulsesTMP{dp}.pulsetimes = [diopulsesTMP{dp}.pulsetimes; diopulses{dp}.pulsetimes];
                            diopulsesTMP{dp}.timesincelast = [diopulsesTMP{dp}.timesincelast; diopulses{dp}.timesincelast];
                            diopulsesTMP{dp}.pulselength = [diopulsesTMP{dp}.pulselength; diopulses{dp}.pulselength];
                            diopulsesTMP{dp}.pulseind = [diopulsesTMP{dp}.pulseind; diopulses{dp}.pulseind];
                        catch
                        end
                    end
                end
                tr = tr +1;
            catch
                tr = 0;
            end
        end
        clear diopulses;
        diopulses = diopulsesTMP;
    else
    load(sprintf('%s/%sDIO%d.mat',outdir,animals,day));
    end
end

%% use txt times for mult an or single an

% load in the text file of [ animal start-hour start-min start-sec end-hour end-min end-sec drinking?]
%and convert it to timestamp rate [animal startstamp endstamp drinking?]
if usetimestxt == 1;
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
        %         Timestamps(i,4) = diostr{1}(8);
        tline = fgetl(fid); %this will get the next line in the txt file..
    end
    fclose(fid)
    clear fid
    
    % for i = epochs;
    if day < 10;
    load(sprintf('%s/%sDIO0%d.mat',outdir,animals,day));
    else
        load(sprintf('%s/%sDIO%d.mat',outdir,animals,day));
    end
    % end
end

%% bit setup

%don't you dare change this order
pairs = [{'AB'} {'BC'} {'DE'} {'EF'} {'CD'} {'BE'} {'AF'} {'BD'} {'CE'}];

%Rew IN [ 1   2    3    4   5   6  H ] wells
%              [16 15 14 13 12 11 10] bits

%Rew Out [  1    2   3   4   5   6    H ] wells
%                  [ 24 23 22 21 20 19 18 ] bits

%LEDs [   1   2   3   4   5    6  H7 H8 ] wells
%           [ 32 31 30 29 28 27 26  25 ] bits

%don't change:
W = [16 15 14 13 12 11 10];
R = [ 24 23 22 21 20 19 18 ];
L = [ 32 31 30 29 28 27 26  25 ];
AL = L(configsetup(1)); AW = W(configsetup(1)); AR = R(configsetup(1));
BL = L(configsetup(2)); BW = W(configsetup(2)); BR = R(configsetup(2));
CL = L(configsetup(3)); CW = W(configsetup(3)); CR = R(configsetup(3));
DL = L(configsetup(4)); DW = W(configsetup(4)); DR = R(configsetup(4));
EL = L(configsetup(5)); EW = W(configsetup(5)); ER = R(configsetup(5));
FL = L(configsetup(6)); FW = W(configsetup(6)); FR = R(configsetup(6));
HL1 = 25; HL2 = 26; HW = 10; HR = 18; %home doesn't change

pairbits = [{AL;BL; AW; BW; AR}...
    {BL;CL; BW; CW; BR}...
    {DL;EL; DW; EW; DR}...
    {EL;FL; EW; FW; ER}...
    {CL;DL; CW;DW;CR}...
    {BL;EL; BW; EW; BR}...
    {AL;FL; AW;FW; AR}...
    {BL;DL; BW;DW;BR}...
    {CL;EL; CW; EW;CR}];

%% Collect DIO pulses according to params
%create var for each bit list that exists.. . i.e.: A1W is arm 1 Well (sensors).. R ~ reward... L ~ light... H is home

alldio = [];

if timewise == 1 | useallpulses == 1;
    for i = 10:32;
        try %this will spit out diopulses with a bit tag in the 2nd col and animal tag in the 3rd col.
            if timewise == 1;
                alldio = [alldio; diopulses{i}.pulsetimes(diopulses{i}.pulsetimes(:,1) > startstamp & diopulses{i}.pulsetimes(:,1) < endstamp,1) ones(length(diopulses{i}.pulsetimes(diopulses{i}.pulsetimes(:,1) > startstamp & diopulses{i}.pulsetimes(:,1) < endstamp,1)),1).*i ones(length(diopulses{i}.pulsetimes(diopulses{i}.pulsetimes(:,1) > startstamp & diopulses{i}.pulsetimes(:,1) < endstamp,1)),1).*animnum]; %collect all stamps and tag with #bit
            elseif useallpulses == 1;
                alldio = [alldio; diopulses{i}.pulsetimes(:,1) ones(length(diopulses{i}.pulsetimes(:,1)),1).*i ones(length(diopulses{i}.pulsetimes(:,1)),1).*animnum]; %collect all stamps and tag with #bit
            end
        end
    end
elseif usetimestxt == 1;
    for ii = double(unique(Timestamps(:,1))'); %for each animal, collect the dio pulses within that time range
        anind = find(Timestamps(:,1) == ii)';% figure out how to separate animals to use later
        for k = [anind]; %cycle over the epochs of this animal
            for i = 10:32; %cycle through all the bits and collect any that exist for this animal for this epoch
                try %this will spit out diopulses in the time range of each epoch with a bit tag in the 2nd col and animal tag in the 3rd col.
                    alldio = [alldio; diopulses{i}.pulsetimes(diopulses{i}.pulsetimes(:,1) > Timestamps(k,2) & diopulses{i}.pulsetimes(:,1) < Timestamps(k,3),1) ones(length(diopulses{i}.pulsetimes(diopulses{i}.pulsetimes(:,1) > Timestamps(k,2) & diopulses{i}.pulsetimes(:,1) < Timestamps(k,3),1)),1).*i ones(length(diopulses{i}.pulsetimes(diopulses{i}.pulsetimes(:,1) > Timestamps(k,2) & diopulses{i}.pulsetimes(:,1) < Timestamps(k,3),1)),1).*ii]; %collect all stamps and tag with #bit and animal
                end
            end
        end
    end
end

alldiosort = sortrows(alldio,1);

%% to get percent correctoutboudnd for each individual pair types

%will spit out performance struct where the first column is the time stamp that the lights of the next pair went on.
%the second col is the time of the next well sensor, regardless of whether it was correct
%the third col is the pair type of that trial... AB = 1, BC = 2, DE = 3, EF = 4, CD = 5, BD = 6, BE = 7, CE = 8, AF = 9.
%fourth col is 1 = correct, 2 = went to the wrong lighted arm, 3 = went to the wrong unlighted arm.
%fifth col is 1 if correct, 0 if incorrect (regardless of incorrect type).. this get's fed into the state space model.

performance = [];
for  j = unique(alldiosort(:,3))'; %for each animal..
    alldiosortAN = alldiosort(alldiosort(:,3) == j,:);
    for i = 1:length(alldiosortAN(alldiosortAN(:,2) == HR));%length(HomR(:,1));
        currhometime = alldiosortAN(max(find(alldiosortAN(:,2) == HR, i)),1);
        currhometime = currhometime - 10000; %if the HR is delayed for .5 seconds after the HW, and the lights for T05 wrongly came on .15 sec after the HW, then HR-5000 should precede lights
%         nextlightstime  = min(alldiosortAN(find(alldiosortAN(:,1)>currhometime & alldiosortAN(:,2) ~=HW,2),1));
        nextlightstime  = min(alldiosortAN(find(alldiosortAN(:,1)>currhometime & alldiosortAN(:,2) >=27 & alldiosortAN(:,2) <=32,2),1));
        nextlightsbits  = alldiosortAN(find(alldiosortAN(:,1)>currhometime &  alldiosortAN(:,2) >=27 & alldiosortAN(:,2) <=32,2),2);
        nextwelltime = alldiosortAN(find(alldiosortAN(:,1)>currhometime & alldiosortAN(:,2) >=11 & alldiosortAN(:,2) <=16,1),1);
        nextwellbit =  alldiosortAN(find(alldiosortAN(:,1)>currhometime & alldiosortAN(:,2) >=11 & alldiosortAN(:,2) <=16,1),2);
        
        endepoch = 0;
        if isempty(nextwelltime); %if it's the end of the day and the file ended on a light
            break
        elseif nextwelltime-nextlightstime > 10000*60*1; %if the end of the epoch... change the last multiplier to extend or reduce epoch latency.. right now it's 1 minute
            endepoch = 1;
        end
        
        
        if endepoch == 0;
            for kk = 1:length(pairbits(1,:));
                try
                    if nextlightsbits == [pairbits{1,kk};pairbits{2,kk}] | nextlightsbits == [pairbits{2,kk};pairbits{1,kk}]
                        if nextwellbit == pairbits{3,kk};
                            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime kk 1 1 j day]; % pair corr corr anim day
                        elseif nextwellbit == pairbits{4,kk};
                            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime kk 2 0 j day]; % pair incorrLit incorr anim day
                        else
                            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime kk 3 0 j day]; % pair incorrUnlit incorr anim day
                        end
                    end
                end
            end
        end
    end
end

%%
% use 'performance' to get sliding window averages
g = 0;
for f = unique(performance(:,7))'; %for each an
    g = g+1;
        for pa = unique(performance(performance(:,7) == f,4))'; %for each pair type
            anpairperf = performance(performance(:,7) == f & performance(:,4) == pa,6); %get performance for each animal and each pair
        k = 0;
        for j = win/2:(length(anpairperf)-win/2);
            k = k+1;
            slidepairs{g}{pa}(k) = mean(anpairperf(j-(win/2-1):j+win/2));
        end
    end
end

%% plotstuffs
ff = 0;
for f = unique(performance(:,7))'; %for each an, create a new figure
    ff = ff +1;
    figure;
    i = 0;
    for pa = unique(performance(performance(:,7) == f,4))';
        i = i+1;
        %first row plot performance pairwise
        %         subplot(length(unique(performance(:,4)))+2,length(unique(performance(:,4))),i);
        subplot(4,length(unique(performance(performance(:,7) == f,4))),i);
        hold on;
        numcorr(i) = length(performance(performance(:,7) == f & performance(:,4) == pa &  performance(:,5) == 1));
        numerrlit(i) = length(performance(performance(:,7) == f & performance(:,4) == pa &  performance(:,5) == 2));
        numerrunlit(i) = length(performance(performance(:,7) == f & performance(:,4) == pa &  performance(:,5) == 3));
%         
%         numcorr(i) = length(performance(performance(performance(:,7) == f,4) == pa & performance(performance(:,7) == f,5) == 1));
%         numerrlit(i) = length(performance(performance(performance(:,7) == f,4) == pa & performance(performance(:,7) == f,5) == 2));
%         numerrunlit(i) = length(performance(performance(performance(:,7) == f,4) == pa & performance(performance(:,7) == f,5) == 3));
%         
        
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
        
        %         subplot(length(unique(performance(:,4)))+2, length(unique(performance(:,4))), ((i+1)*length(unique(performance(:,4))))+1:(((i+1)*length(unique(performance(:,4))))+ length(unique(performance(:,4)))));
        subplot(4, length(unique(performance(performance(:,7) == f,4))),length(unique(performance(performance(:,7) == f,4)))+i);
        [t, p, bmode, b05, b95, bp, background_prob, sige, pmatrix, pc, lt] = DR_getestprobcorrect(performance(performance(:,7) == f & performance(:,4) == pa,6),.5,0);
        plot(t, bmode(2:end),'r-');
        hold on;
        plot(t, b05(2:end),'k', t, b95(2:end), 'k');
        hold on; [y, x] = find(bp > 0);
        h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor','k');
        set(h, 'MarkerEdgeColor', 'k');
        hold on; [y, x] = find(bp == 0);
        h = plot(x,y+0.05,'s'); set(h, 'MarkerFaceColor', [0.75 0.75 0.75]);
        set(h, 'MarkerEdgeColor', 'k');
        %         axis([1 t(end)  0 1.05]);
        %         title([{'IO(0.95)TrialLearn =' num2str(lt)}; {'RW var =' num2str(sige^2)}]);
        axis([1 t(end)  0 1]);
        line([1 t(end)], [background_prob  background_prob], 'Color', [0 0 1]);
        line([1 t(end)], [.75  .75], 'Color', [0 1 0])
        set(gca,'XTick',[0: floor(t(end)/4) :t(end)]);
        %         line([1 t(end)], [background_prob  background_prob ]);
        title(['LearnT:' num2str(lt)]);
        if i == 1;
            xlabel('Trial #')
            ylabel('Prob Correct')
        end
        
        try
            subplot(4,length(unique(performance(performance(:,7) == f,4))), length(unique(performance(performance(:,7) == f,4)))*2+i);
            plot(slidepairs{ff}{pa}, 'r');
            axis([1 length(slidepairs{ff}{pa})  0 1]);
            set(gca, 'YTick', [0:.5:1])
            line([1 length(slidepairs{ff}{pa})], [.75 .75], 'Color', [0 1 0])
            line([1 length(slidepairs{ff}{pa})], [.5 .5], 'Color', [0 0 1])
            set(gca,'XTick', [0:floor(length(slidepairs{ff}{pa})/3):length(slidepairs{ff}{pa})]);
            if i == 1;
                ylabel(sprintf('%%cor TWin:%dt',win))
                %         ylabel('%correct')
            end
        catch
            plot([1:1],[1:1]);
            text(0, .5, sprintf('~enuf trials '));
        end
        
        subplot(4,length(unique(performance(performance(:,7) == f,4))), length(unique(performance(performance(:,7) == f,4)))*3+i);
        % subplot(length(unique(performance(:,4)))+2,length(unique(performance(:,4))), length(unique(performance(:,4)))+i);
        hold on;
        respcorr(i,1) = mean(performance(performance(:,7) == f & performance(:,4) == pa & performance(:,5) == 1,3))/10000; %put into seconds
        resperrlit(i,1) = mean(performance(performance(:,7) == f & performance(:,4) == pa & performance(:,5) ==2,3))/10000;
        resperrunlit(i,1) = mean(performance(performance(:,7) == f & performance(:,4) == pa & performance(:,5) == 3,3))/10000;
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
    supertitle(sprintf('T%d d%d ___ %0.0f%% ___%0.2fs ___ #%d',f,day, 100*allperf, allresp, totaltrials));
    
    %     figfile = [figdir 'Performance' sprintf('%s_d%d_t%d',lower(animals),day,trial)];
    figfile = [figdir sprintf('T%d_d%d',f,day)];
    if savefigs==1
        print('-dpng', figfile,'-r300');
    end
    if cyclefigs == 1;
        keyboard;
    end
    close
end

if saveperformance == 1;
    save(sprintf('%sperformance_d%d', figdir, day), 'performance');
end