
% Similar to combine theta and ripple. Compare run and sleep ripplemodln

% USe the saved data files - HP_ripplemodsleep and HP_ripplemod to plot
% correlations between theta modulation and ripple modulation'

clear;

doindivfigs = 1; % Figures for individual cells. Multiple plots across figures for entire significant population

savefig1=0;
savedir = '/data25/sjadhav/HPExpt/ProcessedData/';

area = 'PFC';
%PFC or CA1 - specify area above
sleepfile = [savedir 'HP_ripplemodsleep_',area,'_gather'];
ripplefile = [savedir 'HP_ripplemod_',area,'_gather'];

load(sleepfile, 'allripplemod','allripplemod_idx');
allsleepmod = allripplemod; allsleepmod_idx = allripplemod_idx;
load(ripplefile, 'allripplemod','allripplemod_idx','pret','postt','binsize','rwin','bckwin');


% Match idxs as in xcorrmesaures2

cntcells=0; cnt_mismatch=0;

% for i=1:length(allripplemod)
%
%     rippleidx = allripplemod_idx(i,:);
%     match = rowfind(rippleidx, allthetamod_idx);

allsigshuf_sleep=[];
allsigshuf_run=[];

for i=1:length(allripplemod)
    
    curridx = allripplemod_idx(i,:);
    match = rowfind(curridx, allsleepmod_idx);
    
    if match~=0,
        cntcells = cntcells+1;
        allmod(cntcells).idx = curridx;
        % Sleep
        allmod(cntcells).sig_shufs = allsleepmod(match).sig_shuf; % Use this to determine significance
        %allmod(cntcells).pshuf = allsleepmod(i).pshuf;
        %allmod(cntcells).D = allsleepmod(i).Dm; % Distance metric
        allmod(cntcells).ripmodln_peaks = allsleepmod(match).modln_peak; % % peak change above/below baseline
        allmod(cntcells).ripmodlns = allsleepmod(match).modln; % Mean change over baseline
        allmod(cntcells).ripmodln_shufs = allsleepmod(match).modln_shuf; % %value of significance
        allmod(cntcells).sig_ttests = allsleepmod(match).sig_ttest;
        allmod(cntcells).sleephist = allsleepmod(match).hist;
        % Run Ripple
        allmod(cntcells).sig_shuf = allripplemod(i).sig_shuf; % Use this to determine significance
        %allmod(cntcells).pshuf = allripplemod(i).pshuf;
        %allmod(cntcells).D = allripplemod(i).Dm; % Distance metric
        allmod(cntcells).ripmodln_peak = allripplemod(i).modln_peak; % % peak change above baseline
        allmod(cntcells).ripmodln = allripplemod(i).modln; % Mean change over baseline
        allmod(cntcells).ripmodln_shuf = allripplemod(i).modln_shuf; % %value of significance
        allmod(cntcells).sig_ttest = allripplemod(i).sig_ttest;
        allmod(cntcells).runhist = allripplemod(i).hist;
        
        allsigshuf_sleep(cntcells) = allsleepmod(i).sig_shuf;
        allsigshuf_run(cntcells) = allripplemod(i).sig_shuf;
        
    end
    
end

% Figures for individual cells

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
        % Plot Line at 0 ms - Onset of stimulation
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
        % Plot Line at 0 ms - Onset of stimulation
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

end  % do indivfigs





% ------------------
% Population Figures
% ------------------
forppr = 0;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Figures/03Nov/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',40);
    tfont = 40;
    xfont = 40;
    yfont = 40;
end

% Get data
for i=1:length(allmod)
    % Sleep
    allripmodln_peaks(i) = allmod(i).ripmodln_peaks;
    allripmodlns(i) = allmod(i).ripmodlns;
    allripmodln_shufs(i) = allmod(i).ripmodln_shufs;
    allsigshufs(i) = allmod(i).sig_shufs;
    allsigttests(i) = allmod(i).sig_ttests;
    % Ripple
    allripmodln_peak(i) = allmod(i).ripmodln_peak;
    allripmodln(i) = allmod(i).ripmodln;
    allripmodln_shuf(i) = allmod(i).ripmodln_shuf;
    allsigshuf(i) = allmod(i).sig_shuf;
    allsigttest(i) = allmod(i).sig_ttest;
end

sigsleep = find(allsigshuf_sleep==1); % Same as, sigsleep = find(allsigshufs==1);  
sigrip = find(allsigshuf_run==1); % Same as, sigrip = find(allsigshuf==1); 

sigboth = find(allsigshufs==1 & allsigshuf==1);

%sigrip = find(allpshuf<0.05);
allsig = union(sigrip, sigsleep);


% Which ripple modln to use
% ------------------------
allripmodln = abs(allripmodln_peak); % Peak change over baseline
allripmodlns = abs(allripmodln_peaks);
%allripmodln = allripmodln_peak; % Peak change over baseline
%allripmodlns = allripmodln_peaks;

%allripmodln = abs(allripmodln); % Mean change over baseline
%allripmodlns = abs(allripmodlns);
%allripmodln = allripmodln; % Mean change over baseline
%allripmodlns = allripmodlns; 


%Sleep vs. Awake Ripple Modulation
%-------------------------------------
figure; hold on;
if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end

% Plot all cells if you want to
plot(allripmodlns, allripmodln, 'k.','MarkerSize',24); % All cells
%plot(allripmodlns(sigsleep), allripmodln(sigsleep), 'k.','MarkerSize',24);
%plot(allripmodlns(sigrip), allripmodln(sigrip), 'k.','MarkerSize',24);

plot(allripmodlns(sigsleep), allripmodln(sigsleep), 'cs','MarkerSize',22, 'LineWidth',2);
plot(allripmodlns(sigrip), allripmodln(sigrip), 'ro','MarkerSize',20,'LineWidth',2);
title(sprintf('%s: Sleep vs Run Ripple Modulation', area),'FontSize',tfont,'Fontweight','normal');
xlabel(['Sleep Ripple Modln'],'FontSize',xfont,'Fontweight','normal');
ylabel(['Awake Ripple Modln'],'FontSize',yfont,'Fontweight','normal');
legend('All Units','Sleep Signf','Run Signf');

xaxis = -100:1:400;
plot(xaxis,zeros(size(xaxis)),'k--','LineWidth',2);

[r1,p1] = corrcoef(abs(allripmodlns),abs(allripmodln))
[r1rt,p1rt] = corrcoef(abs(allripmodlns(sigboth)),abs(allripmodln(sigboth)))
[r1t,p1t] = corrcoef(abs(allripmodlns(sigsleep)),abs(allripmodln(sigsleep)))
[r1r,p1r] = corrcoef(abs(allripmodlns(sigrip)),abs(allripmodln(sigrip))) 

str = ''; rstr = '';
if p1(1,2)<0.05, str = '*'; end
if p1(1,2)<0.01, str = '**'; end
if p1(1,2)<0.001, str = '***'; end
if p1r(1,2)<0.05, rstr = '*'; end
if p1r(1,2)<0.01, rstr = '**'; end
if p1r(1,2)<0.001, rstr = '***'; end

if strcmp(area,'CA1')
    set(gca,'XLim',[-300 6800]); %set(gca,'YLim',[0 1550]);
    text(-200,2500,sprintf('Cells defined in both Run and Sleep: %d', length(allsigshuf)),'FontSize',30,'Fontweight','normal');
    text(1000,1050,sprintf('Signf in both: %d',length(sigboth)),'FontSize',30,'Fontweight','normal');
    text(1000,950,sprintf('Signf in Run: %d',length(sigrip)),'FontSize',30,'Fontweight','normal');
    text(1000,850,sprintf('Signf in Sleep: %d',length(sigsleep)),'FontSize',30,'Fontweight','normal');
else
    set(gca,'XLim',[-100 350]); % Peak
    %set(gca,'XLim',[-80 250]); % Mean
    %set(gca,'YLim',[-80 100]);
    text(0,190,sprintf('Cells defined in both Run and Sleep: %d', length(allsigshuf)),'FontSize',30,'Fontweight','normal');
    text(150,100,sprintf('Signf in both: %d',length(sigboth)),'FontSize',30,'Fontweight','normal');
    text(150,80,sprintf('Signf in Run: %d',length(sigrip)),'FontSize',30,'Fontweight','normal');
    text(150,60,sprintf('Signf in Sleep: %d',length(sigsleep)),'FontSize',30,'Fontweight','normal');
    text(150,40,sprintf('R = %0.2f%s',r1(1,2), str),'FontSize',30,'Fontweight','normal');
    text(150,20,sprintf('Rr = %0.2f%s',r1r(1,2), rstr),'FontSize',30,'Fontweight','normal');
end


 
% Do a shuffle test - both a normal shuffle test, and a regression shuffle
% ------------------------------------------------------------------------

rshuf=[]; pshuf=[];
for i = 1:1000
    rorder = randperm(length(allripmodlns));
    randmodln = allripmodlns(rorder);
    [rtmp, ptmp] = corrcoef(randmodln,allripmodln);
    rshuf(i) = rtmp(1,2); pshuf(i) = ptmp(1,2);
end

length(find(pshuf<0.05));
prctile95 = prctile(rshuf,95); prctile99 = prctile(rshuf,99);
if r1(1,2)>prctile95,
    Sig95 = 1,
else
    Sig95 = 0,
end
if r1(1,2)>prctile99,
    Sig99 = 1,
else
    Sig99 = 0,
end



% Regression
% -----------
[b00,bint00,r00,rint00,stats00] = regress(allripmodln', [ones(size(allripmodlns')) allripmodlns']);
xpts = min(allripmodlns):0.01:max(allripmodlns);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
rsquare = stats00(1);

% Regression for only awake SWR modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allripmodln(sigrip)', [ones(size(allripmodlns(sigrip)')) allripmodlns(sigrip)']);
xpts = min(allripmodlns):0.01:max(allripmodlns);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'r-','LineWidth',4);  % Theta vs Rip - Only SWR significant
rsquare = stats00(1);

% Regression for only sleep SWR modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allripmodln(sigsleep)', [ones(size(allripmodlns(sigsleep)')) allripmodlns(sigsleep)']);
xpts = min(allripmodlns):0.01:max(allripmodlns);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'c-','LineWidth',4);  % Theta vs Rip - Only SWR significant
rsquare = stats00(1);

% Regression for both awake and sleep SWR modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allripmodln(sigboth)', [ones(size(allripmodlns(sigboth)')) allripmodlns(sigboth)']);
xpts = min(allripmodlns):0.01:max(allripmodlns);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'b-','LineWidth',4);  % Theta vs Rip - Only SWR significant
rsquare = stats00(1);

% Do regression after shifting data to make intercept 0?
% ------------------------------------------------------
allripmodln_0 = allripmodln-mean(allripmodln);
allripmodlns_0 = allripmodlns-mean(allripmodlns);
[b0,bint0,r0,rint0,stats0] = regress(allripmodln_0',[ones(size(allripmodlns_0')) allripmodlns_0']);
bfit0 = b0(1)+b0(2)*xpts;
 
rshuffle = []; pshuffle = []; rsquare_shuffle = []; psquare_shuffle = []; b_shuffle = [];


% % Shuffling
% % ---------
for n=1:1000
    rorder = randperm(length(allripmodlns));
    randmodln = allripmodlns(rorder);
    % Get corrcoeff of shuffle
    [rsh,psh] = corrcoef(randmodln, allripmodln);
    r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
    % Get regression of shuffle after making intercept 0 / Or Not
    %shuffle_0 = shuffle - mean(shuffle);
    %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(allripmodlns_0')) allripmodlns_0']);
    [bsh,bintsh,rsh,rintsh,statssh] = regress(allripmodln', [ones(size(randmodln')) randmodln']);
    rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
    b_shuffle(n,:) = bsh;
end

% Significance from shuffle
prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
pshuf2 = length(find(r1(1,2)<r_shuffle))/n

% Get regression corresponding to 99 percentile
idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
plot(xpts,bfitsh,'k--','LineWidth',2);  % Theta vs Rip - 99% shuffle line
% 
% % and 95 %tile
% idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,95));
% idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
% bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
% plot(xpts,bfitsh,'g--','LineWidth',2);  % Theta vs Rip - 95% shuffle line



figfile = [figdir,area,'_SleepVsAwakeRipModln']
if savefig1==1,
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end




% --------------------------------------------------------------------------


