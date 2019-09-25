
% Get output (gather files) from DFSsj_HPexpt_ThetaCovOnly, and and DFSsj_getripalignspiking
% and categorize Theta Cov by degrees of ripple modulation.Can also add theta modllation on top

clear;
savefig1=0;
savedir = '/data25/sjadhav/HPExpt/ProcessedData/';

thrsmodln = 20; % peak/mean change over bckground

% Get PFC ripmod and thetamod files
% ----------------------------------
area = 'PFC';
thetafile = [savedir 'HP_thetamod_',area,'_gather'];
state = ''; statename = 'Run'; 
%state = 'sleep'; statename = 'Sleep'; 
ripplefile = [savedir 'HP_ripplemod',state,'_',area,'_gather']; % ripple mod in awake or sleep
load(ripplefile, 'allripplemod','allripplemod_idx'); % load allripplemod and allripplemod_idx.
load(thetafile,'allthetamod','allthetamod_idx');

% Get Cov file with PFC ripmod cells - awake
% -------------------------------------------

% ------------------------
% 1) PFCripmod cells only
% ------------------------
gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFC',state,'ripmod_gather'];
load(gatherdatafile);

% Match Indexes
% --------------
cntidxs=0; cnt_mismatch=0; Qripmodln = [];

for i=1:length(Qdata)
    curridx = Qdata(i).index([1 2 4 5]); % anim-day-[skip ep]=tet-cell
    match = rowfind(curridx, allripplemod_idx);
    if match~=0,
        cntidxs = cntidxs+1;
        Qmod(cntidxs).idx = curridx;
        % Ripple
        Qripmodln(cntidxs) = abs(allripplemod(match).modln_peak);
        Qmod(cntidxs).sig_shuf = allripplemod(match).sig_shuf; % Use this to determine significance
        %Qmod(cntidxs).pshuf = allripplemod(match).pshuf;
        %Qmod(cntidxs).D = allripplemod(match).D; % Distance metric
        Qmod(cntidxs).ripmodln_peak = allripplemod(match).modln_peak; % % peak change above baseline
        Qmod(cntidxs).ripmodln = allripplemod(match).modln; % Mean change over baseline
        Qmod(cntidxs).ripmodln_shuf = allripplemod(match).modln_shuf; % %value of significance
        Qmod(cntidxs).sig_ttest = allripplemod(match).sig_ttest;
        
        % Theta CrossCov
        Qmod(cntidxs).Zsm = Qdata(i).Zsm;
        Qmod(cntidxs).Z = Qdata(i).Z;
        Qmod(cntidxs).peak = Qdata(i).peak;
        Qmod(cntidxs).trough = Qdata(i).trough;
        
        % Theta
        %         Qmod(cntidxs).sph = allthetamod(i).sph;
        %         Qmod(cntidxs).Nspk = allthetamod(i).Nspk;
        %         Qmod(cntidxs).kappa = allthetamod(i).kappa;
        %         Qmod(cntidxs).modln = allthetamod(i).modln;
        %         Qmod(cntidxs).meanphase = allthetamod(i).meanphase;
        %         Qmod(cntidxs).prayl = allthetamod(i).prayl;
    end
end


% Get mean std. cross cov for categories
% Threshold  from Cov gather files
%thrs = 3.66; % thrs = 2.58 (p=0.01) to thrs = 3.66 (p = 0.01/40)

%thrs=-Inf;

% a) > thrsmodln and <= thrsmodln
% --------------------------------

Qmean1 = []; Qmeansig1 = []; Qmeansig1only = []; Qmean2 = []; Qmeansig2 = []; Qmeansig2only = [];
cntsig = 0; cntsig1 = 0; cntsig2 = 0;
for i = 1:cntidxs
    
    currripmodln = Qripmodln(i);
    currZsm = Qmod(i).Zsm;
    currpeak = Qmod(i).peak; currtrough = Qmod(i).trough;
    npairs = size(currZsm,1);
    sigidx = find(currpeak >= thrs);
    if ~isempty(sigidx)
        cntsig = cntsig+1;
        nsigpairs = length(sigidx);
    end
    
    if currripmodln > thrsmodln
        Qmean1(i,:) = nansum(currZsm,1)./(sqrt(npairs));
        if ~isempty(sigidx)
            cntsig1 = cntsig1+1;
            Qmeansig1(i,:) = nansum(currZsm(sigidx,:),1)./(sqrt(nsigpairs));
            Qmeansig1only(cntsig1,:) = nansum(currZsm(sigidx,:),1)./(sqrt(nsigpairs)); % Discarding cells with no sign CCs for check. Not correct way.
            
        else
            Qmeansig1(i,:) = zeros(size(runcorrtime));
        end
    else
        Qmean2(i,:) = nansum(currZsm,1)./(sqrt(npairs));
        if ~isempty(sigidx)
            cntsig2 = cntsig2+1;
            Qmeansig2(i,:) = nansum(currZsm(sigidx,:),1)./(sqrt(nsigpairs));
            Qmeansig2only(cntsig2,:) = nansum(currZsm(sigidx,:),1)./(sqrt(nsigpairs)); % Discarding cells with no sign CCs for check. Not correct way.
            
        else
            Qmeansig2(i,:) = zeros(size(runcorrtime));
        end
    end % end currripmodln
    
end % end cntidxs


%b) Rip modln vs. Mean Standardized Cross Covariance - Peak
% -------------------------------------------------------

Qmeanp=[]; Qmeansigp=[]; Qmeansig1onlyp=[];
getidx = [];

for i = 1:cntidxs
    getidx(i,:) = Qmod(i).idx;
    currripmodln = Qripmodln(i);
    currZsm = Qmod(i).Zsm;
    currpeak = Qmod(i).peak; currtrough = Qmod(i).trough;
    npairs = size(currZsm,1);
    sigidx = find(currpeak >= thrs);
    %Qmeanp(i,:) = nansum(currZsm,1)./(sqrt(npairs));
    Qmeanp(i) = nansum(currpeak)./(sqrt(npairs));
    
    if ~isempty(sigidx)
        cntsig = cntsig + 1;
        nsigpairs = length(sigidx);
        Qmeansigp(i,:) = nansum(currpeak(sigidx))./(sqrt(nsigpairs));
        Qmeansigonlyp(cntsig1,:) = nansum(currpeak(sigidx))./(sqrt(nsigpairs)); % Discarding cells with no sign CCs for check. Not correct way.
    else
        Qmeansigp(i,:) = 0;
    end
end

cntidxs
cntsig1
cntsig2

x = find(Qmeanp > 4 & Qripmodln > 75 & Qripmodln < 100)
%x = find(Qmeanp > 4 && Qripmodln > 100)
getidx(x,:)



% ------------------------
% % 2) PFCall cells
% % ----------------------
gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFC',state,'all_gather'];

load(gatherdatafile);
% Match Indexes
% --------------
allcntidxs=0; cnt_mismatch=0; allQripmodln = [];

for i=1:length(Qdata)
    curridx = Qdata(i).index([1 2 4 5]); % anim-day-[skip ep]=tet-cell
    match = rowfind(curridx, allripplemod_idx);
    if match~=0,
        allcntidxs = allcntidxs+1;
        allQmod(allcntidxs).idx = curridx;
        % Ripple
        allQripmodln(allcntidxs) = abs(allripplemod(match).modln_peak);
        allQmod(allcntidxs).sig_shuf = allripplemod(match).sig_shuf; % Use this to determine significance
        %allQmod(allcntidxs).pshuf = allripplemod(match).pshuf;
        %allQmod(allcntidxs).D = allripplemod(match).D; % Distance metric
        allQmod(allcntidxs).ripmodln_peak = allripplemod(match).modln_peak; % % peak change above baseline
        allQmod(allcntidxs).ripmodln = allripplemod(match).modln; % Mean change over baseline
        allQmod(allcntidxs).ripmodln_shuf = allripplemod(match).modln_shuf; % %value of significance
        allQmod(allcntidxs).sig_ttest = allripplemod(match).sig_ttest;

        % Theta CrossCov
        allQmod(allcntidxs).Zsm = Qdata(i).Zsm;
        allQmod(allcntidxs).Z = Qdata(i).Z;
        allQmod(allcntidxs).peak = Qdata(i).peak;
        allQmod(allcntidxs).trough = Qdata(i).trough;
    end
end

% a) > thrsmodln and <= thrsmodln
% --------------------------------

%thrs=-Inf;

allQmean1 = []; allQmeansig1 = []; allQmeansig1only = []; allQmean2 = []; allQmeansig2 = []; allQmeansig2only = [];
allcntsig = 0; allcntsig1 = 0; allcntsig2 = 0;
for i = 1:allcntidxs
    
    currripmodln = allQripmodln(i);
    currZsm = allQmod(i).Zsm;
    currpeak = allQmod(i).peak; currtrough = allQmod(i).trough;
    npairs = size(currZsm,1);
    sigidx = find(currpeak >= thrs);
    if ~isempty(sigidx)
        allcntsig = allcntsig+1;
        nsigpairs = length(sigidx);
    end
    
    if currripmodln > thrsmodln
        allQmean1(i,:) = nansum(currZsm,1)./(sqrt(npairs));
        if ~isempty(sigidx)
            allcntsig1 = allcntsig1+1;
            allQmeansig1(i,:) = nansum(currZsm(sigidx,:),1)./(sqrt(nsigpairs));
            allQmeansig1only(allcntsig1,:) = nansum(currZsm(sigidx,:),1)./(sqrt(nsigpairs)); % Discarding cells with no sign CCs for check. Not correct way.
            
        else
            allQmeansig1(i,:) = zeros(size(runcorrtime));
        end
    else
        allQmean2(i,:) = nansum(currZsm,1)./(sqrt(npairs));
        if ~isempty(sigidx)
            allcntsig2 = allcntsig2+1;
            allQmeansig2(i,:) = nansum(currZsm(sigidx,:),1)./(sqrt(nsigpairs));
            allQmeansig2only(allcntsig2,:) = nansum(currZsm(sigidx,:),1)./(sqrt(nsigpairs)); % Discarding cells with no sign CCs for check. Not correct way.
            
        else
            allQmeansig2(i,:) = zeros(size(runcorrtime));
        end
    end % end currripmodln
    
end % end allcntidxs


%b) Rip modln vs. Mean Stdardized Cross Covariance - Peak
% -------------------------------------------------------

allQmeanp=[]; allQmeansigp=[]; allQmeansig1onlyp=[];

for i = 1:allcntidxs
    currripmodln = allQripmodln(i);
    currZsm = allQmod(i).Zsm;
    currpeak = allQmod(i).peak; currtrough = allQmod(i).trough;
    npairs = size(currZsm,1);
    sigidx = find(currpeak >= thrs);
    %allQmeanp(i,:) = nansum(currZsm,1)./(sqrt(npairs));
    allQmeanp(i) = nansum(currpeak)./(sqrt(npairs));
    
    if ~isempty(sigidx)
        allcntsig = allcntsig + 1;
        nsigpairs = length(sigidx);
        allQmeansigp(i,:) = nansum(currpeak(sigidx))./(sqrt(nsigpairs));
        allQmeansigonlyp(allcntsig1,:) = nansum(currpeak(sigidx))./(sqrt(nsigpairs)); % Discarding cells with no sign CCs for check. Not correct way.
    else
        allQmeansigp(i,:) = 0;
    end
end

allcntidxs
allcntsig1
allcntsig2



% Plot
% ----
figdir = '/data25/sjadhav/HPExpt/Figures/03Nov/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',40);
tfont = 40;
xfont = 40;
yfont = 40;

% ---------------------------------
% A) Plot Qmean sig
% -------------------------------

% For PFC ripmodcells
% -------------------------------
Qmeansig1pop = nanmean(Qmeansig1,1);
Qmeansig2pop = nanmean(Qmeansig2,1);
figure; hold on;
plot(runcorrtime, Qmeansig1pop,[ 'r-'],'LineWidth',3);
plot(runcorrtime, Qmeansig2pop,[ 'b-'],'LineWidth',3);
low = min([Qmeansig1pop, Qmeansig2pop]); high = max([Qmeansig1pop, Qmeansig2pop]);
ylabel(sprintf('Mean Std. CrossCov'),'FontSize',xfont);
set(gca,'XLim',[-0.4 0.4]);
xlabel('Time (sec)','FontSize',xfont);
title([sprintf('PFC RimodCells %s. Mean Q. Thrs = %d', statename, thrsmodln)],'FontSize',xfont);
legend('High RipModln','LowRipModln');

low = -0.1; high = 0.75;
line([0 0], [low high],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([0.2 0.2], [low high],'Color',[0.5 0.5 0.5],'LineWidth',1);
line([-0.2 -0.2], [low high],'Color',[0.5 0.5 0.5],'LineWidth',1);
figfile = [figdir,'PFC_RipmodCells_',statename,'_CrossCov_RipMod_HighLow']
if savefig1==1,
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end

% For PFC all cells
% -------------------------------
allQmeansig1pop = nanmean(allQmeansig1,1);
allQmeansig2pop = nanmean(allQmeansig2,1);
figure; hold on;
plot(runcorrtime, allQmeansig1pop,[ 'r-'],'LineWidth',3);
plot(runcorrtime, allQmeansig2pop,[ 'b-'],'LineWidth',3);
low = min([allQmeansig1pop, allQmeansig2pop]); high = max([Qmeansig1pop, Qmeansig2pop]);
ylabel(sprintf('Mean Std. CrossCov'),'FontSize',xfont);
set(gca,'XLim',[-0.4 0.4]);
xlabel('Time (sec)','FontSize',xfont);
title([sprintf('All PFC cells %s. Thrs = %d', statename, thrsmodln)],'FontSize',xfont);
legend('High RipModln','LowRipModln');

low = -0.1; high = 0.75;
line([0 0], [low high],'Color',[0.5 0.5 0.5],'LineWidth',2);
line([0.2 0.2], [low high],'Color',[0.5 0.5 0.5],'LineWidth',1);
line([-0.2 -0.2], [low high],'Color',[0.5 0.5 0.5],'LineWidth',1);
figfile = [figdir,'PFC_AllCells_',statename,'_CrossCov_Allcells_HighLow']
if savefig1==1,
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end



% ---------------------------------
% B) Plot Rip Modln vs Qmeanp
% ---------------------------------

% For PFC ripmodcells
% -------------------------------
figure; hold on;
plot(Qmeanp, Qripmodln,'ro','MarkerSize',20,'LineWidth',2);
legend('PFCRipmodcells');
ylabel(['SWR Modln'],'FontSize',yfont,'Fontweight','normal');
xlabel(['Peak in Std. Cross Cov'],'FontSize',xfont,'Fontweight','normal');
title([sprintf('%s SWR modln vs CC peak', statename)],'FontSize',tfont);


% Corrln
str = ''; 
[r,p] = corrcoef(Qmeanp, Qripmodln)
if p(1,2)<=0.05, str = ''; end
if p(1,2)<0.01, str = '**'; end
if p(1,2)<0.001, str = '***'; end
text(5,100,sprintf('R = %0.2f%s',r(1,2), str),'FontSize',30,'Fontweight','normal');

% Stats for Ripmod cells
% ----------------------
% Do a shuffle test - both a normal shuffle test, and a regression shuffle
rshuf=[]; pshuf=[];
for i = 1:1000
    rorder = randperm(length(Qmeanp));
    randmodln = Qmeanp(rorder);
    [rtmp, ptmp] = corrcoef(randmodln,Qripmodln);
    rshuf(i) = rtmp(1,2); pshuf(i) = ptmp(1,2);
end
length(find(pshuf<0.05));
prctile95 = prctile(rshuf,95); prctile99 = prctile(rshuf,99);
if r(1,2)>prctile95,
    Sig95 = 1,
else
    Sig95 = 0,
end
if r(1,2)>prctile99,
    Sig99 = 1,
else
    Sig99 = 0,
end
% Regression
% -----------
[b00,bint00,r00,rint00,stats00] = regress(Qripmodln', [ones(size(Qmeanp)) Qmeanp]);
xpts = min(Qmeanp):0.01:max(Qmeanp);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'r-','LineWidth',4);  % Theta vs Rip
rsquare = stats00(1);
% % Regression Shuffling
% % --------------------
for n=1:1000
    rorder = randperm(length(Qmeanp));
    randmodln = Qmeanp(rorder);
    randmodln = randmodln';
    % Get corrcoeff of shuffle
    [rsh,psh] = corrcoef(randmodln, Qripmodln);
    r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
    % Get regression of shuffle after making intercept 0 / Or Not
    %shuffle_0 = shuffle - mean(shuffle);
    %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(Qmeanp_0')) Qmeanp_0']);
    [bsh,bintsh,rsh,rintsh,statssh] = regress(Qripmodln', [ones(size(randmodln')) randmodln']);
    rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
    b_shuffle(n,:) = bsh;
end
% Significance from shuffle
prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
pshuf2 = length(find(r(1,2)<r_shuffle))/n
% Get regression corresponding to 99 percentile
idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
plot(xpts,bfitsh,'r--','LineWidth',2);  % Theta vs Rip - 99% shuffle line

figfile = [figdir,'PFC_RipmodCells_',statename,'_CrossCov_Vs_RipMod']
if savefig1==1,
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end



% For PFC all cells
% -------------------------------
figure; hold on;
plot(allQmeanp, allQripmodln,'k.','MarkerSize',24);
legend('PFCallcells');
ylabel(['SWR Modln'],'FontSize',yfont,'Fontweight','normal');
xlabel(['Peak in Std. Cross Cov'],'FontSize',xfont,'Fontweight','normal');
title([sprintf('%s SWR modln vs CC peak', statename)],'FontSize',tfont);


allstr = ''; 
[allr,allp] = corrcoef(allQmeanp, allQripmodln)
if allp(1,2)<=0.05, str = ''; end
if allp(1,2)<0.01, str = '**'; end
if allp(1,2)<0.001, str = '***'; end
text(5,80,sprintf('allR = %0.2f%s',allr(1,2), str),'FontSize',30,'Fontweight','normal');

% Do a shuffle test - both a normal shuffle test, and a regression shuffle
rshuf=[]; pshuf=[];
for i = 1:1000
    rorder = randperm(length(allQmeanp));
    randmodln = allQmeanp(rorder);
    [rtmp, ptmp] = corrcoef(randmodln,allQripmodln);
    rshuf(i) = rtmp(1,2); pshuf(i) = ptmp(1,2);
end
length(find(pshuf<0.05));
prctile95 = prctile(rshuf,95); prctile99 = prctile(rshuf,99);
if allr(1,2)>prctile95,
    Sig95 = 1,
else
    Sig95 = 0,
end
if allr(1,2)>prctile99,
    Sig99 = 1,
else
    Sig99 = 0,
end
% Regression
% -----------
[b00,bint00,r00,rint00,stats00] = regress(allQripmodln', [ones(size(allQmeanp)) allQmeanp]);
xpts = min(allQmeanp):0.01:max(allQmeanp);
bfit00 = b00(1)+b00(2)*xpts;
plot(xpts,bfit00,'k-','LineWidth',4);  % Theta vs Rip
rsquare = stats00(1);
% % Regression Shuffling
% % --------------------
for n=1:1000
    rorder = randperm(length(allQmeanp));
    randmodln = allQmeanp(rorder);
    randmodln = randmodln';
    % Get corrcoeff of shuffle
    [rsh,psh] = corrcoef(randmodln, allQripmodln);
    r_shuffle(n) = rsh(1,2); p_shuffle(n) = psh(1,2);
    % Get regression of shuffle after making intercept 0 / Or Not
    %shuffle_0 = shuffle - mean(shuffle);
    %[bsh,bintsh,rsh,rintsh,statssh] = regress(shuffle_0', [ones(size(allQmeanp_0')) allQmeanp_0']);
    [bsh,bintsh,rsh,rintsh,statssh] = regress(allQripmodln', [ones(size(randmodln')) randmodln']);
    rsquare_shuffle(n) = statssh(1); preg_shuffle(n) = statssh(3);
    b_shuffle(n,:) = bsh;
end
% Significance from shuffle
prctile(rsquare_shuffle,99); prctile(r_shuffle,99); %figure; hist(r_shuffle,50); hist(rsquare_shuffle,50);
pshuf2 = length(find(allr(1,2)<r_shuffle))/n
% Get regression corresponding to 99 percentile
idxs=find(rsquare_shuffle>=prctile(rsquare_shuffle,99));
idx=idxs(find(rsquare_shuffle(idxs)==min(rsquare_shuffle(idxs))));
bfitsh = b_shuffle(idx,1)+b_shuffle(idx,2)*xpts;
plot(xpts,bfitsh,'k--','LineWidth',2);  % Theta vs Rip - 99% shuffle line

figfile = [figdir,'PFC_AllCells_',statename,'_CrossCov_Vs_RipMod'];
if savefig1==1,
    print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
end


% -------------------------  Filter Format Done -------------------------



% ----------------------------------
% Whether to gather data or to load previously gathered data
% --------------------------------------------------------------------
gatherdata = 0; savegatherdata = 0;

% CA1all files
% ----
%gatherdatafile = [savedir 'HP_thetacov_CA1allPFCripmod_gather']; area = 'CA1allPFC'; kind = 'allvsripmod'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1allPFCripunmod_gather']; area = 'CA1allPFC'; kind = 'allvsripunmod'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1allPFCthetamod_gather']; area = 'CA1allPFC'; kind = 'allvsthetamod'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1allPFCthetaunmod_gather']; area = 'CA1allPFC'; kind = 'allvsthetaunmod'; state = '';

% Compare CA1 (theta modulated only) and PFC (theta modulated vs unmodulated)
% -------------------------------------------------------------------------
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCthetamod_gather']; area = 'CA1thetamodPFCthetamod'; kind = 'thetamod'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCthetaunmod_gather']; area = 'CA1thetamodPFCthetaunmod'; kind = 'thetaunmod'; state = '';
% %Both unmodulated
%gatherdatafile = [savedir 'HP_thetacov_CA1thetaunmodPFCthetaunmod_gather']; area = 'CA1thetaunmodPFCthetaunmod'; kind = 'unmod'; state = '';


% % Compare CA1 (theta modulated only) and PFC ripandtheta mod vs unmod
% ------------------------------------------------------------------------
%gatherdatafile = [savedir 'HP_thetacov_CA1thetaPFCboth_gather']; area = 'CA1thetaPFCboth'; kind = 'both'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1thetaPFCunboth_gather']; area = 'CA1thetaPFCboth'; kind = 'both'; state = '';

% % CA1 (theta modulated only) and PFC ripmod vs ripunmod. Also Compare to Sleep Ripmod computed in other Script
% ------------------------------------------------------------------------
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCripmod_gather']; area = 'CA1thetamodPFCripmod'; kind = 'thetamodripmod'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCripunmod_gather']; area = 'CA1thetamodPFCripunmod'; kind = 'thetamodripunmod'; state = '';
% ------
% Sleep RipMod
% ------
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCsleepripmod_gather']; area = 'CA1thetamodPFCsleepripmod'; kind = 'thetamodsleepripmod'; state = '';
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCsleepripunmod_gather']; area = 'CA1thetamodPFCsleepripunmod'; kind = 'thetamodsleepripunmod'; state = '';


% % Theta and Ripple Mod mix
% --------------------------
% PFC only theta mod, but not ripple mod
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFConlythetamod_gather']; area = 'CA1thetamodPFConlythetamod'; kind = 'thetamod'; state = '';
% PFC only ripple mod, but not theta mod
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFConlyripmod_gather']; area = 'CA1thetamodPFConlyripmod'; kind = 'thetamod'; state = '';

% % CA1 ThetaMod vs PFC all - for soarting by ripmodln later
% ----------------------------------------------------------
%gatherdatafile = [savedir 'HP_thetacov_CA1thetamodPFCall_gather']; area = 'CA1thetamodPFCall'; kind = 'thetamod'; state = '';












% ------------------------------------------------------------------
% COMBINING PLOTS ACROSS FILES
% ------------------------------------------------------------------

