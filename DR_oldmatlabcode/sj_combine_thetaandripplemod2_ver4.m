
% Ver4 : Starting 10Feb2014 - Sync codes with everyone
% Currently, this is same as sj_combine_thetaandripplemod2_ver4var. Change change this to plot ripplemod using
% avg firng rate criterion instead of var criterion

% Version 2 - USing filter framework ripplemod 

% USe the saved data files - HP_thetamod and HP_ripplemod to plot
% correlations between theta modulation and ripple modulation'

clear;

savefig1=0;
%savedir = '/data15/gideon/ProcessedData/';
savedir = '/data25/sjadhav/HPExpt/ProcessedData/';

area = 'PFC'; 
%PFC or CA1 - specify area above
thetafile = [savedir 'HP_thetamod_',area,'_alldata_Nspk50_gather_2-12-2014']; 
%thetafile = [savedir 'HP_thetamod_',area,'_alldata_gather']; 
state = ''; %state = 'sleep'; %or state = '';
ripplefile = [savedir 'HP_ripplemod',state,'_',area,'_alldata_std3_speed4_ntet2_Nspk50_gather_2-12-2014']; % ripple mod in awake or sleep

load(ripplefile, 'allripplemod','allripplemod_idx'); % load allripplemod and allripplemod_idx. 
load(thetafile,'allthetamod','allthetamod_idx'); 

if strcmp(state,'sleep'),
    statename = 'Sleep';
else
    statename = 'Run';
end

% Match idxs as in xcorrmesaures2

cntcells=0; cnt_mismatch=0;

% for i=1:length(allripplemod)
%     
%     rippleidx = allripplemod_idx(i,:);
%     match = rowfind(rippleidx, allthetamod_idx);
    
for i=1:length(allthetamod)
    i;
    curridx = allthetamod_idx(i,:);
    
    % DONT SHIFT ANY DAYS
    % Temp - till Gideon's animal Nadal is fixed in DFSsj_plotthetamod
%     if curridx(1)==4
%         curridx(2)=curridx(2)-7; %Adjust day
%     end
%     
    match = rowfind(curridx, allripplemod_idx);
    
    if match~=0,    
        cntcells = cntcells+1;
        allmod(cntcells).idx = curridx;
        % Theta
        allmod(cntcells).sph = allthetamod(i).sph;
        allmod(cntcells).Nspk = allthetamod(i).Nspk;
        allmod(cntcells).kappa = allthetamod(i).kappa;
        allmod(cntcells).modln = allthetamod(i).modln;
        allmod(cntcells).meanphase = allthetamod(i).meanphase;
        allmod(cntcells).prayl = allthetamod(i).prayl;
        % Ripple
        allmod(cntcells).sig_shuf = allripplemod(match).sig_shuf; % Use this to determine significance
        %allmod(cntcells).pshuf = allripplemod(match).pshuf;
        %allmod(cntcells).D = allripplemod(match).D; % Distance metric
        allmod(cntcells).ripmodln_peak = allripplemod(match).modln_peak; % % peak change above baseline
        allmod(cntcells).ripmodln = allripplemod(match).modln; % Mean change over baseline
        allmod(cntcells).ripmodln_shuf = allripplemod(match).modln_shuf; % %value of significance
        allmod(cntcells).sig_ttest = allripplemod(match).sig_ttest;
         
        % New
        %allmod(cntcells).modln_raw = allripplemod(match).modln_raw;
        allmod(cntcells).ripmodln_div = allripplemod(match).modln_div; % % peak change above baseline
        allmod(cntcells).pvar = allripplemod(match).rasterShufP;
        allmod(cntcells).mvar = allripplemod(match).varRespAmp;
        allmod(cntcells).mvar2 = allripplemod(match).varRespAmp2;
        allmod(cntcells).pvar2 = allripplemod(match).rasterShufP2;

        %allmod(cntcells).mvar3 = allripplemod(match).var_changerespbck;
        
        % Prop
        match
        allmod(cntcells).cellfr = allripplemod(match).cellfr;
        
        % Rdm resp modln
        allmod(cntcells).rdmmodln_peak = allripplemod(match).modln_peak_rdm; % % peak change above baseline
        allmod(cntcells).rdmmodln = allripplemod(match).modln_rdm; % Mean change over baseline
        
    end
    
end

% ------------------
% Population Figures
% ------------------
forppr = 0; 
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Figures/Jan2014/';
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
    % Theta
    allidxs(i,:) = allmod(i).idx;
    allkappas(i) = allmod(i).kappa;
    allmodln(i) = allmod(i).modln;
    allmeanphase(i) = allmod(i).meanphase;
    allprayl(i) = allmod(i).prayl;      
    % Ripple
    allripmodln_peak(i) = allmod(i).ripmodln_peak;
    allripmodln_div(i) = allmod(i).ripmodln_div;
    allripmodln(i) = allmod(i).ripmodln;
    allripmodln_shuf(i) = allmod(i).ripmodln_shuf; %[]
    %allD(i) = allmod(i).D;
    %allpshuf(i) = allmod(i).pshuf;      
    allsigshuf(i) = allmod(i).sig_shuf;
    allsigttest(i) = allmod(i).sig_ttest;
    
    % New
    allpvar(i) = allmod(i).pvar;
    allripmodln_var(i) = allmod(i).mvar;
    allripmodln_var2(i) = allmod(i).mvar2;
    %allripmodln_var3(i) = allmod(i).mvar3;
    allpvar2(i) = allmod(i).pvar2;
    %allmodln_raw(i) = allmod(i).modln_raw;
    
    % Prop
    allcellfr(i) = allmod(i).cellfr;
    
    % Rdm resp
    allrdmmodln_peak(i) = allmod(i).rdmmodln_peak;
    allrdmmodln(i) = allmod(i).rdmmodln;
    
    
    
end



% WhichSWR modln to use
% ------------------------
%allripmodln = allripmodln_peak; % Peak change over baseline 
%allripmodln = allripmodln; % Mean change over baseline
%allripmodln = abs(allripmodln); % Mean change over baseline

%allripmodln = abs(allripmodln_peak); sigrip = find(allsigshuf==1); sigboth = find(allprayl<0.05 & allsigshuf==1); allrdmmodln = abs(allrdmmodln_peak); % 
allripmodln = abs(allripmodln_var2); 
% %Remove outlier in allripmodln_var
% if ~strcmp(state,'sleep')
%     rem = find(allripmodln>12);
%     allripmodln(rem)=[]; allmodln(rem)=[]; allkappas(rem)=[]; allprayl(rem)=[]; allpvar(rem)=[]; allpvar2(rem)=[]; allsigshuf(rem)=[];
%     allcellfr(rem) = [];
% end

sigrip = find(allpvar2<0.05); sigtheta = find(allprayl<0.05); sigboth = find(allprayl<0.05 & allpvar2<0.05==1);
allsig = union(sigrip, sigtheta); 




% Vector of 0s and 1s - sig ripple vs sig theta
ripvec = zeros(size(allripmodln)); ripvec(sigrip)=1;
thetavec = zeros(size(allripmodln)); thetavec(sigtheta)=1;
[rvec,pvec] = corrcoef(ripvec,thetavec);

%x = find(allmodln > 0.65 & allmodln < 0.7 & allripmodln > 80);
%allidxs(x,:);



 

% -------------------------------% -------------------------------% -------------------------------
% -------------------------------% -------------------------------% -------------------------------



% -------------------------------
%Theta Modln vs. Ripple Modulation
%-------------------------------------
figure; hold on;
if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end

plot(allmodln, allripmodln, 'k.','MarkerSize',24); % Plot all if you want to

%plot(allmodln(sigtheta), allripmodln(sigtheta), 'k.','MarkerSize',24);
%plot(allmodln(sigrip), allripmodln(sigrip), 'k.','MarkerSize',24);
plot(allmodln(sigtheta), allripmodln(sigtheta), 'cs','MarkerSize',22,'LineWidth',2);
plot(allmodln(sigrip), allripmodln(sigrip), 'ro','MarkerSize',20,'LineWidth',2);
title(sprintf('%s: Theta vs %s Ripple Modulation', area, statename),'FontSize',tfont,'Fontweight','normal');
xlabel(['Theta Modln'],'FontSize',xfont,'Fontweight','normal');
ylabel(sprintf('%s Ripple Modln',statename),'FontSize',yfont,'Fontweight','normal');
legend('All Units','Sig Theta Phlock',sprintf('Sig %s Ripple',statename));

[r2,p2] = corrcoef(allmodln,allripmodln)  
[r2rt,p2rt] = corrcoef(allmodln(allsig),allripmodln(allsig)) 
[r2t,p2t] = corrcoef(allmodln(sigtheta),allripmodln(sigtheta))
[r2r,p2r] = corrcoef(allmodln(sigrip),allripmodln(sigrip)) 

xaxis = 0.5:0.1:1;
plot(xaxis,zeros(size(xaxis)),'k--','LineWidth',2);

str = ''; rstr = ''; tstr = ''; rtstr = '';
if p2(1,2)<0.05, str = '*'; end
if p2(1,2)<0.01, str = '**'; end
if p2(1,2)<0.001, str = '***'; end
if p2r(1,2)<0.05, rstr = '*'; end
if p2r(1,2)<0.01, rstr = '**'; end
if p2r(1,2)<0.001, rstr = '***'; end
if p2t(1,2)<0.05, tstr = '*'; end
if p2t(1,2)<0.01, tstr = '**'; end
if p2t(1,2)<0.001, tstr = '***'; end
if p2rt(1,2)<0.05, rtstr = '*'; end
if p2rt(1,2)<0.01, rtstr = '**'; end
if p2rt(1,2)<0.001, rtstr = '***'; end


if strcmp(area,'CA1')
    if ~strcmp(state,'sleep')
        set(gca,'XLim',[0.55 1])
        text(0.56,1000,sprintf('Signf in both: %d',length(sigboth)),'FontSize',30,'Fontweight','normal');
        text(0.56,900,sprintf('Signf in Rip: %d',length(sigrip)),'FontSize',30,'Fontweight','normal');
        text(0.56,800,sprintf('Signf in Theta: %d',length(sigtheta)),'FontSize',30,'Fontweight','normal');
        text(0.56,700,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
    else
        set(gca,'XLim',[0.55 1])
        text(0.8,1450,sprintf('Signf in both: %d',length(sigboth)),'FontSize',30,'Fontweight','normal');
        text(0.8,1300,sprintf('Signf in Rip: %d',length(sigrip)),'FontSize',30,'Fontweight','normal');
        text(0.8,1200,sprintf('Signf in Theta: %d',length(sigtheta)),'FontSize',30,'Fontweight','normal');
        text(0.8,1100,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
    end
else
    if ~strcmp(state,'sleep'); 
        set(gca,'XLim',[0.5 1]); 
        %set(gca,'YLim',[-120 200]); % Modln peak
        %set(gca,'YLim',[-80 100]); % Modln
        text(0.55,40,sprintf('Total Valid Cells: %d', length(allsigshuf)),'FontSize',30,'Fontweight','normal');
        text(0.82,35,sprintf('Signf in both: %d',length(sigboth)),'FontSize',30,'Fontweight','normal');
        text(0.82,30,sprintf('Signf in Rip: %d',length(sigrip)),'FontSize',30,'Fontweight','normal');
        text(0.82,25,sprintf('Signf in Theta: %d',length(sigtheta)),'FontSize',30,'Fontweight','normal');
        text(0.82,20,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
        text(0.82,5,sprintf('Rr = %0.2f%s',r2r(1,2), rstr),'FontSize',30,'Fontweight','normal');
        text(0.72,5,sprintf('Rt = %0.2f%s',r2t(1,2), tstr),'FontSize',30,'Fontweight','normal');
        text(0.92,5,sprintf('Rrt = %0.2f%s',r2rt(1,2), rtstr),'FontSize',30,'Fontweight','normal');


    else
        set(gca,'XLim',[0.5 1.]);  
        %set(gca,'YLim',[-120 500]); % Modln peak
        % set(gca,'YLim',[-80 250]); % Modln
        text(0.6,40,sprintf('Total Valid Cells: %d', length(allsigshuf)),'FontSize',30,'Fontweight','normal');
        text(0.75,35,sprintf('Signf in both: %d',length(sigboth)),'FontSize',30,'Fontweight','normal');
        text(0.75,30,sprintf('Signf in Sleep SWR: %d',length(sigrip)),'FontSize',30,'Fontweight','normal');
        text(0.75,25,sprintf('Signf in Theta: %d',length(sigtheta)),'FontSize',30,'Fontweight','normal');
        text(0.75,20,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
        text(0.75,5,sprintf('Rr = %0.2f%s',r2r(1,2), rstr),'FontSize',30,'Fontweight','normal');
        text(0.65,5,sprintf('Rt = %0.2f%s',r2t(1,2), tstr),'FontSize',30,'Fontweight','normal');
        text(0.85,5,sprintf('Rrt = %0.2f%s',r2rt(1,2), rtstr),'FontSize',30,'Fontweight','normal');


    end
end

figfile = [figdir,area,sprintf('_ThetaModlnVs%sRipModln',statename)]
if savefig1==1,
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end



 
% Do a shuffle test - both a normal shuffle test, and a regression shuffle
% ------------------------------------------------------------------------

rshuf=[]; pshuf=[];
for i = 1:1000
    rorder = randperm(length(allmodln));
    randmodln = allmodln(rorder);
    [rtmp, ptmp] = corrcoef(randmodln,allripmodln);
    rshuf(i) = rtmp(1,2); pshuf(i) = ptmp(1,2);
end

length(find(pshuf<0.05));
prctile95 = prctile(rshuf,95); prctile99 = prctile(rshuf,99);
if r2(1,2)>prctile95,
    Sig95 = 1,
else
    Sig95 = 0,
end
if r2(1,2)>prctile99,
    Sig99 = 1,
else
    Sig99 = 0,
end



% Regression
% -----------
[b00,bint00,r00,rint00,stats00] = regress(allripmodln', [ones(size(allmodln')) allmodln']);
xpts = min(allmodln):0.01:max(allmodln);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
rsquare = stats00(1);

% Regression for only SWR modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allripmodln(sigrip)', [ones(size(allmodln(sigrip)')) allmodln(sigrip)']);
xpts = min(allmodln):0.01:max(allmodln);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'r-','LineWidth',4);  % Theta vs Rip - Only SWR significant
rsquare = stats00(1);

% Regression for only theta modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allripmodln(sigtheta)', [ones(size(allmodln(sigtheta)')) allmodln(sigtheta)']);
xpts = min(allmodln):0.01:max(allmodln);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'c-','LineWidth',4);  % Theta vs Rip - Only Theta significant
rsquare = stats00(1);

% Regression for both modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allripmodln(sigboth)', [ones(size(allmodln(sigboth)')) allmodln(sigboth)']);
xpts = min(allmodln):0.01:max(allmodln);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'b-','LineWidth',4);  % Theta vs Rip - Both SWR and Theta significant
rsquare = stats00(1);



% Do regression after shifting data to make intercept 0?
% ------------------------------------------------------
allripmodln_0 = allripmodln-mean(allripmodln);
allmodln_0 = allmodln-mean(allmodln);
[b0,bint0,r0,rint0,stats0] = regress(allripmodln_0',[ones(size(allmodln_0')) allmodln_0']);
bfit0 = b0(1)+b0(2)*xpts;
 
rshuffle = []; pshuffle = []; rsquare_shuffle = []; psquare_shuffle = []; b_shuffle = [];

% % Shuffling
% % ---------
for n=1:1000
    rorder = randperm(length(allmodln));
    randmodln = allmodln(rorder);
    % Get corrcoeff of shuffle
    [rsh,psh] = corrcoef(randmodln, allripmodln);
    r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
    % Get regression of shuffle after making intercept 0 / Or Not
    %shuffle_0 = shuffle - mean(shuffle);
    %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(allmodln_0')) allmodln_0']);
    [bsh,bintsh,rsh,rintsh,statssh] = regress(allripmodln', [ones(size(randmodln')) randmodln']);
    rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
    b_shuffle(n,:) = bsh;
end

% Significance from shuffle
prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
pshuf2 = length(find(r2(1,2)<r_shuffle))/n

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
























% -------------------------------% -------------------------------% -------------------------------
% -------------------------------% -------------------------------% -------------------------------



%Kappas vs. Ripple Modulation
%-------------------------------------
figure; hold on;
if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end

plot(allkappas, allripmodln, 'k.','MarkerSize',24); % Plot all
%plot(allkappas(sigtheta), allripmodln(sigtheta), 'k.','MarkerSize',24);
%plot(allkappas(sigrip), allripmodln(sigrip), 'k.','MarkerSize',24);
plot(allkappas(sigtheta), allripmodln(sigtheta), 'cs','MarkerSize',22,'LineWidth',2);
plot(allkappas(sigrip), allripmodln(sigrip), 'ro','MarkerSize',20,'LineWidth',2);
title(sprintf('%s: Theta vs %s Ripple Modulation', area, statename),'FontSize',tfont,'Fontweight','normal');
xlabel(['Kappa (Theta Conc Parm)'],'FontSize',xfont,'Fontweight','normal');
ylabel(sprintf('%s Ripple Modln',statename),'FontSize',yfont,'Fontweight','normal');
legend('All Units','Sig Theta Phlock',sprintf('Sig %s Ripple',statename));

[r1,p1] = corrcoef(allkappas,allripmodln) 
[r1rt,p1rt] = corrcoef(allkappas(sigboth),allripmodln(sigboth))
[r1t,p1t] = corrcoef(allkappas(sigtheta),allripmodln(sigtheta))
[r1r,p1r] = corrcoef(allkappas(sigrip),allripmodln(sigrip)) 

xaxis = 0:0.1:1.5;
plot(xaxis,zeros(size(xaxis)),'k--','LineWidth',2);

str = ''; rstr = ''; tstr = ''; rtstr = '';
if p1(1,2)<0.05, str = '*'; end
if p1(1,2)<0.01, str = '**'; end
if p1(1,2)<0.001, str = '***'; end
if p1r(1,2)<0.05, rstr = '*'; end
if p1r(1,2)<0.01, rstr = '**'; end
if p1r(1,2)<0.001, rstr = '***'; end
if p1t(1,2)<0.05, tstr = '*'; end
if p1t(1,2)<0.01, tstr = '**'; end
if p1t(1,2)<0.001, tstr = '***'; end
if p1rt(1,2)<0.05, rtstr = '*'; end
if p1rt(1,2)<0.01, rtstr = '**'; end
if p1rt(1,2)<0.001, rtstr = '***'; end


if strcmp(area,'CA1')
    %set(gca,'XLim',[0 2200]); set(gca,'YLim',[0 1600]);
    if ~strcmp(state,'sleep')
        text(1.1,900,sprintf('Signf in both: %d',length(sigboth)),'FontSize',30,'Fontweight','normal');
        text(1.1,800,sprintf('Signf in Run: %d',length(sigrip)),'FontSize',30,'Fontweight','normal');
        text(1.1,700,sprintf('Signf in Theta: %d',length(sigtheta)),'FontSize',30,'Fontweight','normal');
        text(1.1,600,sprintf('R = %0.2f%s',r1(1,2), str),'FontSize',30,'Fontweight','normal');
    end
else % PFC
    if ~strcmp(state,'sleep')
        set(gca,'XLim',[-0.02 1]);  
        %set(gca,'YLim',[-120 200]); % Modln peak
        %set(gca,'YLim',[-80 90]); % Modln
        text(0.05,40,sprintf('Total Valid Cells: %d', length(allsigshuf)),'FontSize',30,'Fontweight','normal');
        text(0.62,35,sprintf('Signf in both: %d',length(sigboth)),'FontSize',30,'Fontweight','normal');
        text(0.62,30,sprintf('Signf in Ripple: %d',length(sigrip)),'FontSize',30,'Fontweight','normal');
        text(0.62,25,sprintf('Signf in Theta: %d',length(sigtheta)),'FontSize',30,'Fontweight','normal');
        text(0.62,20,sprintf('R = %0.2f%s',r1(1,2), str),'FontSize',30,'Fontweight','normal');
        text(0.62,5,sprintf('Rr = %0.2f%s',r1r(1,2), rstr),'FontSize',30,'Fontweight','normal');
        text(0.62,10,sprintf('Rt = %0.2f%s',r1t(1,2), tstr),'FontSize',30,'Fontweight','normal');
        text(0.62,15,sprintf('Rrt = %0.2f%s',r1rt(1,2), rtstr),'FontSize',30,'Fontweight','normal');

    else
        set(gca,'XLim',[-0.02 1.1]);  
        %set(gca,'YLim',[-120 500]); % Modln peak
        %set(gca,'YLim',[-80 250]); % Modln
        text(0.2,40,sprintf('Total Valid Cells: %d', length(allsigshuf)),'FontSize',30,'Fontweight','normal');
        text(0.4,35,sprintf('Signf in both: %d',length(sigboth)),'FontSize',30,'Fontweight','normal');
        text(0.4,30,sprintf('Signf in Sleep SWR: %d',length(sigrip)),'FontSize',30,'Fontweight','normal');
        text(0.4,25,sprintf('Signf in Theta: %d',length(sigtheta)),'FontSize',30,'Fontweight','normal');
        text(0.4,20,sprintf('R = %0.2f%s',r1(1,2), str),'FontSize',30,'Fontweight','normal');
        text(0.4,5,sprintf('Rr = %0.2f%s',r1r(1,2), rstr),'FontSize',30,'Fontweight','normal');
        text(0.4,10,sprintf('Rt = %0.2f%s',r1t(1,2), tstr),'FontSize',30,'Fontweight','normal');
        text(0.4,15,sprintf('Rrt = %0.2f%s',r1rt(1,2), rtstr),'FontSize',30,'Fontweight','normal');

    end
end

figfile = [figdir,area,sprintf('_ThetaKappaVs%sRipModln',statename)];
if savefig1==1,   
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end



% Do a shuffle test - both a normal shuffle test, and a regression shuffle
% ------------------------------------------------------------------------
rshuf=[]; pshuf=[];
for i = 1:1000
    rorder = randperm(length(allkappas));
    randmodln = allkappas(rorder);
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
[b00,bint00,r00,rint00,stats00] = regress(allripmodln', [ones(size(allkappas')) allkappas']);
xpts = min(allkappas):0.01:max(allkappas);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
rsquare = stats00(1);


% Regression for only SWR modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allripmodln(sigrip)', [ones(size(allkappas(sigrip)')) allkappas(sigrip)']);
xpts = min(allkappas):0.01:max(allkappas);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'r-','LineWidth',4);  % Theta vs Rip - Only SWR significant
rsquare = stats00(1);


% Do regression after shifting data to make intercept 0
% ------------------------------------------------------
allripmodln_0 = allripmodln-mean(allripmodln);
allkappas_0 = allkappas-mean(allkappas);
[b0,bint0,r0,rint0,stats0] = regress(allripmodln_0',[ones(size(allkappas_0')) allkappas_0']);
bfit0 = b0(1)+b0(2)*xpts;
 
rshuffle = []; pshuffle = []; rsquare_shuffle = []; psquare_shuffle = []; b_shuffle = [];

% % Shuffling
% % ---------
for n=1:1000
    rorder = randperm(length(allkappas));
    randmodln = allkappas(rorder);
    % Get corrcoeff of shuffle
    [rsh,psh] = corrcoef(randmodln, allripmodln);
    r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
    % Get regression of shuffle after making intercept 0 / Or Not
    %shuffle_0 = shuffle - mean(shuffle);
    %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(allkappas_0')) allkappas_0']);
    [bsh,bintsh,rsh,rintsh,statssh] = regress(allripmodln', [ones(size(randmodln')) randmodln']);
    rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
    b_shuffle(n,:) = bsh;
end

% Significance from shuffle
prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
pshuf1 = length(find(r2(1,2)<r_shuffle))/n

% Get regression corresponding to 99 percentile
idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
plot(xpts,bfitsh,'k--','LineWidth',2);  % Theta vs Rip - 99% shuffle line

% % and 95 %tile
% idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,95));
% idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
% bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
% plot(xpts,bfitsh,'g--','LineWidth',2);  % Theta vs Rip - 95% shuffle line





figure; hold on;
plot(allcellfr,allripmodln,'bo');
plot(allcellfr(sigrip),allripmodln(sigrip),'rx');
title('Rip Modln vs Fir Rate')
[r,p] = corrcoef(allcellfr,allripmodln)
[rs,ps] = corrcoef(allcellfr(sigrip),allripmodln(sigrip))

% figure; hold on;
% plot(allcellfr,allrdmmodln,'bo');
% plot(allcellfr(sigrip),allrdmmodln(sigrip),'rx');
% title('Rdm Resp Modln vs Fir Rate')
% 
% 
% figure; hold on;
% plot(allcellfr,allmodln,'bo');
% plot(allcellfr(sigrip),allmodln(sigrip),'rx');
% title('Theta Modln vs Fir Rate')
% 
% 
% figure; hold on;
% plot(allcellfr,allkappas,'bo');
% plot(allcellfr(sigrip),allkappas(sigrip),'rx');
% title('Theta Kappas vs Fir Rate')






keyboard;




% -------------------------------
% -------------------------------
% -------------------------------
% -------------------------------
% -------------------------------


% -------------------------------
%Theta Kappas vs. Ripple Modulation - Radom Resp Data
%-------------------------------------
figure; hold on;
if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end

plot(allkappas, allrdmmodln, 'k.','MarkerSize',24); % Plot all if you want to

%plot(allkappas(sigtheta), allrdmmodln(sigtheta), 'k.','MarkerSize',24);
%plot(allkappas(sigrip), allrdmmodln(sigrip), 'k.','MarkerSize',24);
plot(allkappas(sigtheta), allrdmmodln(sigtheta), 'cs','MarkerSize',22,'LineWidth',2);
plot(allkappas(sigrip), allrdmmodln(sigrip), 'ro','MarkerSize',20,'LineWidth',2);
title(sprintf('%s: Theta vs %s Rdm Resp Modulation', area, statename),'FontSize',tfont,'Fontweight','normal');
xlabel(['Theta Modln'],'FontSize',xfont,'Fontweight','normal');
ylabel(sprintf('%s Ripple Modln',statename),'FontSize',yfont,'Fontweight','normal');
legend('All Units','Sig Theta Phlock',sprintf('Sig %s Ripple',statename));

[r2,p2] = corrcoef(allkappas,allrdmmodln)  
[r2rt,p2rt] = corrcoef(allkappas(allsig),allrdmmodln(allsig)) 
[r2t,p2t] = corrcoef(allkappas(sigtheta),allrdmmodln(sigtheta))
[r2r,p2r] = corrcoef(allkappas(sigrip),allrdmmodln(sigrip)) 

xaxis = 0.5:0.1:1;
plot(xaxis,zeros(size(xaxis)),'k--','LineWidth',2);

str = ''; rstr = ''; tstr = ''; rtstr = '';
if p2(1,2)<0.05, str = '*'; end
if p2(1,2)<0.01, str = '**'; end
if p2(1,2)<0.001, str = '***'; end
if p2r(1,2)<0.05, rstr = '*'; end
if p2r(1,2)<0.01, rstr = '**'; end
if p2r(1,2)<0.001, rstr = '***'; end
if p2t(1,2)<0.05, tstr = '*'; end
if p2t(1,2)<0.01, tstr = '**'; end
if p2t(1,2)<0.001, tstr = '***'; end
if p2rt(1,2)<0.05, rtstr = '*'; end
if p2rt(1,2)<0.01, rtstr = '**'; end
if p2rt(1,2)<0.001, rtstr = '***'; end


if strcmp(area,'CA1')
    if ~strcmp(state,'sleep')
        set(gca,'XLim',[0.55 1])
    
        text(0.56,700,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
    else
        set(gca,'XLim',[0.55 1])
        text(0.8,1450,sprintf('Signf in both: %d',length(sigboth)),'FontSize',30,'Fontweight','normal');
        text(0.8,1300,sprintf('Signf in Rip: %d',length(sigrip)),'FontSize',30,'Fontweight','normal');
        text(0.8,1200,sprintf('Signf in Theta: %d',length(sigtheta)),'FontSize',30,'Fontweight','normal');
        text(0.8,1100,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
    end
else
    if ~strcmp(state,'sleep'); 
        set(gca,'XLim',[0.5 1]); 
        %set(gca,'YLim',[-120 200]); % Modln peak
        %set(gca,'YLim',[-80 100]); % Modln
        text(0.82,20,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
        text(0.82,5,sprintf('Rr = %0.2f%s',r2r(1,2), rstr),'FontSize',30,'Fontweight','normal');
        text(0.72,5,sprintf('Rt = %0.2f%s',r2t(1,2), tstr),'FontSize',30,'Fontweight','normal');
        text(0.92,5,sprintf('Rrt = %0.2f%s',r2rt(1,2), rtstr),'FontSize',30,'Fontweight','normal');


    else
        set(gca,'XLim',[0.5 1.]);  
        %set(gca,'YLim',[-120 500]); % Modln peak
        % set(gca,'YLim',[-80 250]); % Modln
        text(0.75,20,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
        text(0.75,5,sprintf('Rr = %0.2f%s',r2r(1,2), rstr),'FontSize',30,'Fontweight','normal');
        text(0.65,5,sprintf('Rt = %0.2f%s',r2t(1,2), tstr),'FontSize',30,'Fontweight','normal');
        text(0.85,5,sprintf('Rrt = %0.2f%s',r2rt(1,2), rtstr),'FontSize',30,'Fontweight','normal');


    end
end

figfile = [figdir,area,sprintf('_ThetaModlnVs%sRdmModln',statename)]
if savefig1==1,
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end



 
% Do a shuffle test - both a normal shuffle test, and a regression shuffle
% ------------------------------------------------------------------------

rshuf=[]; pshuf=[];
for i = 1:1000
    rorder = randperm(length(allkappas));
    randmodln = allkappas(rorder);
    [rtmp, ptmp] = corrcoef(randmodln,allrdmmodln);
    rshuf(i) = rtmp(1,2); pshuf(i) = ptmp(1,2);
end

length(find(pshuf<0.05));
prctile95 = prctile(rshuf,95); prctile99 = prctile(rshuf,99);
if r2(1,2)>prctile95,
    Sig95 = 1,
else
    Sig95 = 0,
end
if r2(1,2)>prctile99,
    Sig99 = 1,
else
    Sig99 = 0,
end



% Regression
% -----------
[b00,bint00,r00,rint00,stats00] = regress(allrdmmodln', [ones(size(allkappas')) allkappas']);
xpts = min(allkappas):0.01:max(allkappas);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
rsquare = stats00(1);

% Regression for only SWR modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allrdmmodln(sigrip)', [ones(size(allkappas(sigrip)')) allkappas(sigrip)']);
xpts = min(allkappas):0.01:max(allkappas);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'r-','LineWidth',4);  % Theta vs Rip - Only SWR significant
rsquare = stats00(1);

% Regression for only theta modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allrdmmodln(sigtheta)', [ones(size(allkappas(sigtheta)')) allkappas(sigtheta)']);
xpts = min(allkappas):0.01:max(allkappas);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'c-','LineWidth',4);  % Theta vs Rip - Only Theta significant
rsquare = stats00(1);

% Regression for both modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allrdmmodln(sigboth)', [ones(size(allkappas(sigboth)')) allkappas(sigboth)']);
xpts = min(allkappas):0.01:max(allkappas);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'b-','LineWidth',4);  % Theta vs Rip - Both SWR and Theta significant
rsquare = stats00(1);



% Do regression after shifting data to make intercept 0?
% ------------------------------------------------------
allrdmmodln_0 = allrdmmodln-mean(allrdmmodln);
allkappas_0 = allkappas-mean(allkappas);
[b0,bint0,r0,rint0,stats0] = regress(allrdmmodln_0',[ones(size(allkappas_0')) allkappas_0']);
bfit0 = b0(1)+b0(2)*xpts;
 
rshuffle = []; pshuffle = []; rsquare_shuffle = []; psquare_shuffle = []; b_shuffle = [];

% % Shuffling
% % ---------
for n=1:1000
    rorder = randperm(length(allkappas));
    randmodln = allkappas(rorder);
    % Get corrcoeff of shuffle
    [rsh,psh] = corrcoef(randmodln, allrdmmodln);
    r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
    % Get regression of shuffle after making intercept 0 / Or Not
    %shuffle_0 = shuffle - mean(shuffle);
    %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(allkappas_0')) allkappas_0']);
    [bsh,bintsh,rsh,rintsh,statssh] = regress(allrdmmodln', [ones(size(randmodln')) randmodln']);
    rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
    b_shuffle(n,:) = bsh;
end

% Significance from shuffle
prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
pshuf2 = length(find(r2(1,2)<r_shuffle))/n

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








% -------------------------------
%Theta Modln vs. Ripple Modulation - Radom Resp Data
%-------------------------------------
figure; hold on;
if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end

plot(allmodln, allrdmmodln, 'k.','MarkerSize',24); % Plot all if you want to

%plot(allmodln(sigtheta), allrdmmodln(sigtheta), 'k.','MarkerSize',24);
%plot(allmodln(sigrip), allrdmmodln(sigrip), 'k.','MarkerSize',24);
plot(allmodln(sigtheta), allrdmmodln(sigtheta), 'cs','MarkerSize',22,'LineWidth',2);
plot(allmodln(sigrip), allrdmmodln(sigrip), 'ro','MarkerSize',20,'LineWidth',2);
title(sprintf('%s: Theta vs %s Rdm Resp Modulation', area, statename),'FontSize',tfont,'Fontweight','normal');
xlabel(['Theta Modln'],'FontSize',xfont,'Fontweight','normal');
ylabel(sprintf('%s Ripple Modln',statename),'FontSize',yfont,'Fontweight','normal');
legend('All Units','Sig Theta Phlock',sprintf('Sig %s Ripple',statename));

[r2,p2] = corrcoef(allmodln,allrdmmodln)  
[r2rt,p2rt] = corrcoef(allmodln(allsig),allrdmmodln(allsig)) 
[r2t,p2t] = corrcoef(allmodln(sigtheta),allrdmmodln(sigtheta))
[r2r,p2r] = corrcoef(allmodln(sigrip),allrdmmodln(sigrip)) 

xaxis = 0.5:0.1:1;
plot(xaxis,zeros(size(xaxis)),'k--','LineWidth',2);

str = ''; rstr = ''; tstr = ''; rtstr = '';
if p2(1,2)<0.05, str = '*'; end
if p2(1,2)<0.01, str = '**'; end
if p2(1,2)<0.001, str = '***'; end
if p2r(1,2)<0.05, rstr = '*'; end
if p2r(1,2)<0.01, rstr = '**'; end
if p2r(1,2)<0.001, rstr = '***'; end
if p2t(1,2)<0.05, tstr = '*'; end
if p2t(1,2)<0.01, tstr = '**'; end
if p2t(1,2)<0.001, tstr = '***'; end
if p2rt(1,2)<0.05, rtstr = '*'; end
if p2rt(1,2)<0.01, rtstr = '**'; end
if p2rt(1,2)<0.001, rtstr = '***'; end


if strcmp(area,'CA1')
    if ~strcmp(state,'sleep')
        set(gca,'XLim',[0.55 1])
    
        text(0.56,700,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
    else
        set(gca,'XLim',[0.55 1])
        text(0.8,1450,sprintf('Signf in both: %d',length(sigboth)),'FontSize',30,'Fontweight','normal');
        text(0.8,1300,sprintf('Signf in Rip: %d',length(sigrip)),'FontSize',30,'Fontweight','normal');
        text(0.8,1200,sprintf('Signf in Theta: %d',length(sigtheta)),'FontSize',30,'Fontweight','normal');
        text(0.8,1100,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
    end
else
    if ~strcmp(state,'sleep'); 
        set(gca,'XLim',[0.5 1]); 
        %set(gca,'YLim',[-120 200]); % Modln peak
        %set(gca,'YLim',[-80 100]); % Modln
        text(0.82,20,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
        text(0.82,5,sprintf('Rr = %0.2f%s',r2r(1,2), rstr),'FontSize',30,'Fontweight','normal');
        text(0.72,5,sprintf('Rt = %0.2f%s',r2t(1,2), tstr),'FontSize',30,'Fontweight','normal');
        text(0.92,5,sprintf('Rrt = %0.2f%s',r2rt(1,2), rtstr),'FontSize',30,'Fontweight','normal');


    else
        set(gca,'XLim',[0.5 1.]);  
        %set(gca,'YLim',[-120 500]); % Modln peak
        % set(gca,'YLim',[-80 250]); % Modln
        text(0.75,20,sprintf('R = %0.2f%s',r2(1,2), str),'FontSize',30,'Fontweight','normal');
        text(0.75,5,sprintf('Rr = %0.2f%s',r2r(1,2), rstr),'FontSize',30,'Fontweight','normal');
        text(0.65,5,sprintf('Rt = %0.2f%s',r2t(1,2), tstr),'FontSize',30,'Fontweight','normal');
        text(0.85,5,sprintf('Rrt = %0.2f%s',r2rt(1,2), rtstr),'FontSize',30,'Fontweight','normal');


    end
end

figfile = [figdir,area,sprintf('_ThetaModlnVs%sRdmModln',statename)]
if savefig1==1,
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end



 
% Do a shuffle test - both a normal shuffle test, and a regression shuffle
% ------------------------------------------------------------------------

rshuf=[]; pshuf=[];
for i = 1:1000
    rorder = randperm(length(allmodln));
    randmodln = allmodln(rorder);
    [rtmp, ptmp] = corrcoef(randmodln,allrdmmodln);
    rshuf(i) = rtmp(1,2); pshuf(i) = ptmp(1,2);
end

length(find(pshuf<0.05));
prctile95 = prctile(rshuf,95); prctile99 = prctile(rshuf,99);
if r2(1,2)>prctile95,
    Sig95 = 1,
else
    Sig95 = 0,
end
if r2(1,2)>prctile99,
    Sig99 = 1,
else
    Sig99 = 0,
end



% Regression
% -----------
[b00,bint00,r00,rint00,stats00] = regress(allrdmmodln', [ones(size(allmodln')) allmodln']);
xpts = min(allmodln):0.01:max(allmodln);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
rsquare = stats00(1);

% Regression for only SWR modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allrdmmodln(sigrip)', [ones(size(allmodln(sigrip)')) allmodln(sigrip)']);
xpts = min(allmodln):0.01:max(allmodln);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'r-','LineWidth',4);  % Theta vs Rip - Only SWR significant
rsquare = stats00(1);

% Regression for only theta modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allrdmmodln(sigtheta)', [ones(size(allmodln(sigtheta)')) allmodln(sigtheta)']);
xpts = min(allmodln):0.01:max(allmodln);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'c-','LineWidth',4);  % Theta vs Rip - Only Theta significant
rsquare = stats00(1);

% Regression for both modulated cells
% --------------------------------------------
[b00,bint00,r00,rint00,stats00] = regress(allrdmmodln(sigboth)', [ones(size(allmodln(sigboth)')) allmodln(sigboth)']);
xpts = min(allmodln):0.01:max(allmodln);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'b-','LineWidth',4);  % Theta vs Rip - Both SWR and Theta significant
rsquare = stats00(1);



% Do regression after shifting data to make intercept 0?
% ------------------------------------------------------
allrdmmodln_0 = allrdmmodln-mean(allrdmmodln);
allmodln_0 = allmodln-mean(allmodln);
[b0,bint0,r0,rint0,stats0] = regress(allrdmmodln_0',[ones(size(allmodln_0')) allmodln_0']);
bfit0 = b0(1)+b0(2)*xpts;
 
rshuffle = []; pshuffle = []; rsquare_shuffle = []; psquare_shuffle = []; b_shuffle = [];

% % Shuffling
% % ---------
for n=1:1000
    rorder = randperm(length(allmodln));
    randmodln = allmodln(rorder);
    % Get corrcoeff of shuffle
    [rsh,psh] = corrcoef(randmodln, allrdmmodln);
    r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
    % Get regression of shuffle after making intercept 0 / Or Not
    %shuffle_0 = shuffle - mean(shuffle);
    %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(allmodln_0')) allmodln_0']);
    [bsh,bintsh,rsh,rintsh,statssh] = regress(allrdmmodln', [ones(size(randmodln')) randmodln']);
    rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
    b_shuffle(n,:) = bsh;
end

% Significance from shuffle
prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
pshuf2 = length(find(r2(1,2)<r_shuffle))/n

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




% 


