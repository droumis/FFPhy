

warning('off')
clear all;
extract = 0;
trial = 1; %use a single trial for extraction
animal = 'T01';
phase = 2;
day = 17;
savefigs = 1;
cyclefigs = 0;
epochs = [1 2 3 4]; % use multiple epochs for analysis
dir = sprintf('/data19/droumis/TransInf/%s/Day%d/%s_d%d_t%d_raw',lower(animal),day,lower(animal),day,trial);
outdir  = sprintf('/data19/droumis/TransInf/%s/Day%d/%s_d%d_t%d',lower(animal),day,lower(animal),day,trial); ;% '/data19/droumis/TransInf/t03_28_1_out';
figdir = sprintf('/data19/droumis/TransInf/%s/Day%d/',lower(animal),day);


%%%%%%%%%A  B  C  D  E  F
configsetup =  [  6  2   4  5   3  1 ] ;%which track is the corresponding letter on

pairs = [{'AB'} {'BC'} {'DE'} {'EF'}];

%%

if extract == 1;
    mkdir(dir);
    mkdir(outdir);
    copyfile(sprintf('/data19/droumis/TransInf/demetrissnow19/%s_ph%d_d%d_t%d.snow.dat',animal,phase, day, trial), dir);
    cd(dir)
    system('nspike_extract -all ./*');
    generateTimesFileFromCont;
    diodayprocess(dir,outdir,animal,day);
    return
end

% for i = epochs;
load(sprintf('%s/%sDIO%d.mat',outdir,animal,day));
% end

%%
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


%% to get percent correct outbound for 1 day... can get rid of most of his except for alldiosort and HomR
%first pass at % correct.. not accounting for inbound errors.
%# outer rew out / #home rew out

%create var for each bit list that exists.. . i.e.: A1W is arm 1 Well (sensors).. R ~ reward... L ~ light... H is home
%also gather all types of bits separately
% load(sprintf('%s/%sDIO0%d.mat',outdir,animal,day));

bitexist = [];
alldio = [];
allrew = [];
allwell = [];
allled = [];
for i = 10:32;
    try
        alldio = [alldio; diopulses{i}.pulsetimes(:,1) ones(length(diopulses{i}.pulsetimes(:,1)),1).*i]; %collect all stamps and tag with #bit
        if i >=11 & i <=16;
            allwell = [allwell; diopulses{i}.pulsetimes(:,1) ones(length(diopulses{i}.pulsetimes(:,1)),1).*i]; %collect all stamps and tag with #bit
            eval(['T' num2str(i-10) 'W' '=diopulses{i}.pulsetimes(:,1);'])
            eval(['bitexist(i)' '= 1;']) %this should create a list of bits that exist
        elseif i == 10;
            allwell = [allwell; diopulses{i}.pulsetimes(:,1) ones(length(diopulses{i}.pulsetimes(:,1)),1).*i]; %collect all stamps and tag with #bit
            eval(['HomW' '=diopulses{i}.pulsetimes(:,1);'])
            eval(['bitexist(i)' '= 1;']) %this should create a list of bits that exist
        elseif i == 18;
            allrew = [allrew; diopulses{i}.pulsetimes(:,1) ones(length(diopulses{i}.pulsetimes(:,1)),1).*i]; %collect all stamps and tag with #bit
            eval(['HomR' '=diopulses{i}.pulsetimes(:,1);'])
            eval(['bitexist(i)' '= 1;']) %this should create a list of bits that exist
        elseif i >=19 & i <=24;
            allrew = [allrew; diopulses{i}.pulsetimes(:,1) ones(length(diopulses{i}.pulsetimes(:,1)),1).*i]; %collect all stamps and tag with #bit
            eval(['T' num2str(i-18) 'R' '=diopulses{i}.pulsetimes(:,1);'])
            eval(['bitexist(i)' '= 1;']) %this should create a list of bits that exist
        elseif i >=27 & i <=32;
            allled = [allled; diopulses{i}.pulsetimes(:,1) ones(length(diopulses{i}.pulsetimes(:,1)),1).*i]; %collect all stamps and tag with #bit
            eval(['T' num2str(i-26) 'L' '=diopulses{i}.pulsetimes(:,1);'])
            eval(['bitexist(i)' '= 1;']) %this should create a list of bits that exist
        elseif i == 25 | i == 26;
            allled = [allled; diopulses{i}.pulsetimes(:,1) ones(length(diopulses{i}.pulsetimes(:,1)),1).*i]; %collect all stamps and tag with #bit
            eval(['Hom' num2str(i-24) 'L' '=diopulses{i}.pulsetimes(:,1);'])
            eval(['bitexist(i)' '= 1;']) %this should create a list of bits that exist
        end
    catch
        eval(['bitexist(i)' '= 0;'])
    end
end

alldiosort = sortrows(alldio,1);
allrewsort = sortrows(allrew,1);
allledsort = sortrows(allled,1);
allwellsort = sortrows(allwell,1);

%get num outer arms rewarded
[o indexist] = find(bitexist == 1);
rewOutexist = indexist(indexist>= 19 & indexist <= 24);
numcorrect = 0;
for i = rewOutexist;
    numcorrect = numcorrect + length(diopulses{i}.pulsetimes(:,1));
end
%num outer arms rewards / num home rewards
percentCorrect = numcorrect/length(diopulses{18}.pulsetimes(:,1));


%% to get percent correctoutboudnd for each individual pair types

%will spit out performance struct where the first column is the time stamp that the lights of the next pair went on. 
%the second col is the time of the next well sensor, regardless of whether it was correct
%the third col is the pair type of that trial... AB = 1, BC = 2, DE = 3, EF = 4, CD = 5, BD = 6, BE = 7, CE = 8, AF = 9.
%fourth col is 1 = correct, 2 = went to the wrong lighted arm, 3 = went to the wrong unlighted arm. 

performance = [];
for i = 1:length(HomR(:,1));
    currhometime = alldiosort(max(find(alldiosort(:,2) == HR, i)),1);
    nextlightstime  = min(alldiosort(find(alldiosort(:,1)>currhometime & alldiosort(:,2) ~=10,2),1));
    nextlightsbits  = alldiosort(find(alldiosort(:,1)>currhometime & alldiosort(:,2) ~=10,2),2);
    nextwelltime = alldiosort(find(alldiosort(:,1)>currhometime & alldiosort(:,2) >=11 & alldiosort(:,2) <=16,1),1);
    nextwellbit =  alldiosort(find(alldiosort(:,1)>currhometime & alldiosort(:,2) >=11 & alldiosort(:,2) <=16,1),2);
    
    if isempty(nextwelltime);
        break
    end
        
    if nextlightsbits == [AL;BL] | nextlightsbits == [BL;AL]
        if nextwellbit == AW;
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 1 1]; % 1 ~ pair AB. 1 ~ correct
        elseif nextwellbit == BW;
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 1 2]; % 1 ~ pair AB. 2 ~ error wrong light
        else
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 1 3]; % 1 ~ pair AB. 3 ~ error no light
        end
    elseif nextlightsbits == [BL;CL] | nextlightsbits == [CL;BL];
        if nextwellbit == BW;
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 2 1]; % 2 ~ pair BC. 1 ~ correct
        elseif nextwellbit == CW;
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 2 2]; % 2 ~ pair BC. 2 ~ error wrong light
        else
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 2 3]; % 2 ~ pair BC. 3 ~ error no light
        end
    elseif nextlightsbits == [DL;EL] | nextlightsbits == [DL;EL];
        if nextwellbit == DW;
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 3 1]; % 3 ~ pair DE. 1 ~ correct
        elseif nextwellbit == EW;
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 3 2]; % 3 ~ pair DE. 2 ~ error wrong light
        else
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 3 3]; % 3 ~ pair DE. 3 ~ error no light
        end
    elseif nextlightsbits == [EL;FL] | nextlightsbits == [FL;EL];
        if nextwellbit == EW;
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 4 1]; % 4 ~ pair EF. 1 ~ correct
        elseif nextwellbit == FW;
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 4 2]; % 4 ~ pair EF. 2 ~ error wrong light
        else
            performance = [performance; nextlightstime nextwelltime nextwelltime-nextlightstime 4 3]; % 4 ~ pair EF. 3 ~ error no light
        end
    end
end

figure;
for i = 1:length(unique(performance(:,4)));
    subplot(1,length(unique(performance(:,4))),i);
    hold on;
    numcorr(i) = length(performance(performance(:,4) == i & performance(:,5) == 1));
    numerrlit(i) = length(performance(performance(:,4) == i & performance(:,5) == 2));
    numerrunlit(i) = length(performance(performance(:,4) == i & performance(:,5) == 3));
    bar(1, numcorr(i), 'g'); 
    bar(2, numerrlit(i), 'b'); 
    bar(3, numerrunlit(i), 'r');
    pairperf(i) = numcorr(i)/(numcorr(i)+numerrlit(i)+numerrunlit(i));
    title(sprintf('%s(%0.2f)', pairs{i},pairperf(i)));
end
totaltrials = (sum(numcorr)+sum(numerrlit)+sum(numerrunlit));
allperf = sum(numcorr)/totaltrials;
supertitle(sprintf('Correct %s d%d t%d (%0.2f) #%d',animal,day,trial, allperf,totaltrials));


figfile = [figdir 'Performance' sprintf('%s_d%d_t%d',lower(animal),day,trial)];
if savefigs==1
    print('-djpeg', figfile,'-r300');
end
if cyclefigs == 1;
    keyboard;
end
close

figure;
for i = 1:length(unique(performance(:,4)));
    subplot(1,length(unique(performance(:,4))),i);
    hold on;
    respcorr(i,1) = mean(performance(performance(:,4) == i & performance(:,5) == 1,3))/10000; %put into seconds
    resperrlit(i,1) = mean(performance(performance(:,4) == i & performance(:,5) == 2,3))/10000;
    resperrunlit(i,1) = mean(performance(performance(:,4) == i & performance(:,5) == 3,3))/10000;
    if isnan(respcorr(i)); respcorr(i) = 0; end
    if isnan(resperrlit(i)); resperrlit(i) = 0; end
    if isnan(resperrunlit(i)); resperrunlit(i) = 0; end
    bar(1, respcorr(i), 'g');
    bar(2, resperrlit(i), 'b'); 
    bar(3, resperrunlit(i), 'r');
    pairperf(i) = mean([respcorr(i) resperrlit(i) resperrunlit(i)]);
    title(sprintf('%s(%0.2f)', pairs{i},pairperf(i)));
end
allperf = mean([respcorr; resperrlit; resperrunlit]);
supertitle(sprintf('RespTime %s d%d t%d (%0.2f)',animal,day,trial, allperf));

figfile = [figdir 'ResponseRate' sprintf('%s_d%d_t%d',lower(animal),day,trial)];
if savefigs==1
    print('-djpeg', figfile,'-r300');
end
if cyclefigs == 1;
    keyboard;
end
close


%% get performance from entire day


%% get performance across days



