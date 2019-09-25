%function showCorrel
%function showCorrel(rerun,spsel)
function showCorrel(spsel)
%  
% Compute and analyze cross correlations between pairs of cells.
% Shantanu: Sen Chengs from Lorens directory

%load ripple/SCI_data
%SCIAUX= AUX;
%clear sel

addpath /home/chengs/theta/common

rerun= 1;
plotall= 0;
cmpCCG= 0;
showSeqComp= 1;
showMax= 0;
showExcessCorr= 0;
showExcessCorrBar= 0;
showThetaPower= 0;
showAsymIndex= 0;
showArea= 0;
showMean= 0;
showStd= 0;
showPeakLoc= 0;
showBigMax= 0;
%plottype= 'corr-r';
plottype= 'cdf';
%plottype= 'bar';

pairsel= 'overlap';  % 'none', 'overlap', 'all'
spsel= 'placefield_nonripple';   % spike selection {'all', 'traj', 'sat', 'traj_sat', 'run'}
%spsel= 'placefield'; 
%spsel= 'all'; 
opt.scaling= 'coef';  % 'coef', 'none'

opt.dt= 0.002;      % [sec] timestep, code only works with 2ms
opt.T= 2;         % [sec] +/- width of CCG window
opt.cutfreq= 25;    % cutoff frequency [Hz] 
opt.minspikes= 0;  % minimum number of spikes to include
opt.ref= 1;
%opt.timeid= 'acc';

id= 'pf9';
opt.id= id;

plots= {};
%plots{1}.range= [1:8]; plots{1}.outname= [id '_' spsel ];
%plots{1}.plotcol=[.7*[1 1 1]; zeros(3,3); [1 0 0]; [0 1 0]; [0 0 1]; [0 0 0]];
%plots{1}.plot= plottype;

%plots{1}.range= [5]; 
%plots{1}.plotcol=[[1 0 0]];
%plots{1}.plot= plottype;

%plots{1}.range= [2,5]; plots{1}.outname= [id '_' pairsel '_' spsel '_day1'];
%plots{1}.plotcol=[[0 0 0]; [1 0 0]];
%plots{1}.plot= plottype;
%plots{2}.range= [3,6]; plots{2}.outname= [id '_' pairsel '_' spsel '_day2'];
%plots{2}.plotcol=[[0 0 0]; [0 1 0]];
%plots{2}.plot= plottype;
%plots{3}.range= [4,7]; plots{3}.outname= [id '_' pairsel '_' spsel '_day3'];
%plots{3}.plotcol=[[0 0 0]; [0 0 1]];
%plots{3}.plot= plottype;
%plots{4}.range= [1,2:4]; plots{4}.outname= [id '_' pairsel '_' spsel '_famArm'];
%plots{4}.plotcol=[0.7*[1 1 1]; [1 0 0]; [0 1 0]; [0 0 1]];
%plots{4}.plot= plottype;

%plots{5}.range= [5:7]; plots{5}.outname= [id '_' pairsel '_' spsel '_novelArm'];
%plots{5}.plotcol=[[1 0 0]; [0 1 0]; [0 0 1]];
%plots{5}.plot= plottype;

plots{6}.range= [8,5:7]; plots{6}.outname= [id '_' pairsel '_' spsel '_famArm123'];
plots{6}.plotcol=[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]];
plots{6}.plot= plottype;
%plots{6}.range= [1,8,5:7]; plots{6}.outname= [id '_' pairsel '_' spsel '_famArm123'];
%plots{6}.plotcol=[.7*[1 1 1]; [0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]];
%plots{6}.plot= plottype;
%plots{7}.range= [8,2:4]; plots{7}.outname= [id '_' pairsel '_' spsel '_famArmInt'];
%plots{7}.plotcol=[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]];
%plots{7}.plot= plottype;
%plots{8}.range= [1,5:7]; plots{8}.outname= [id '_' pairsel '_' spsel '_fam123'];
%plots{8}.plotcol=[.7*[1 1 1]; [1 0 0]; [0 1 0]; [0 0 1]];
%plots{8}.plot= plottype;

%plots{8}.range= [1:4]; plots{8}.outname= [id '_' pairsel '_' spsel '_famArm'];
%plots{8}.plotcol=[.4*[1 1 1]; [1 0 0]; [0 1 0]; [0 0 1]; [1 1 1]];
%plots{8}.plot= plottype;

%plots{9}.range= [1,8]; plots{9}.outname= [id '_' pairsel '_' spsel '_famfamArm'];
%plots{9}.plotcol=[.7*[1 1 1]; [0 0 0]; [0 1 0]; [0 0 1]];
%plots{9}.plot= plottype;

%plots{8}.range= [8]; plots{8}.outname= [id '_' pairsel '_' spsel '_famArmAll'];
%plots{8}.plotcol=[[0 0 0]];
%plots{8}.plot= plottype;

%sel{1}.selectid= [id '-famArm'];
%sel{1}.pairsel= pairsel;
%sel{1}.spsel= spsel;
%sel{1}.t= [1:150];
%sel{1}.label= 'fam-150';

sel{1}.selectid= [id '-fam'];
sel{1}.pairsel= pairsel;
sel{1}.spsel= spsel;
sel{1}.label= 'fam';

sel{2}.selectid= [id '-famArm'];
sel{2}.pairsel= pairsel;
sel{2}.spsel= spsel;
sel{2}.novelDay= 1;
sel{2}.label= 'fam-1';

sel{3}.selectid= [id '-famArm'];
sel{3}.pairsel= pairsel;
sel{3}.spsel= spsel;
sel{3}.novelDay= 2;
sel{3}.label= 'fam-2';

sel{4}.selectid= [id '-famArm'];
sel{4}.pairsel= pairsel;
sel{4}.spsel= spsel;
sel{4}.novelDay= 3;
sel{4}.label= 'fam-3';

sel{5}.selectid= [id '-novelArm'];
sel{5}.pairsel= pairsel;
sel{5}.spsel= spsel;
sel{5}.novelDay= 1;
sel{5}.label= 'novel-1';

sel{6}.selectid= [id '-novelArm'];
sel{6}.pairsel= pairsel;
sel{6}.spsel= spsel;
sel{6}.novelDay= 2;
sel{6}.label= 'novel-2';

sel{7}.selectid= [id '-novelArm'];
sel{7}.pairsel= pairsel;
sel{7}.spsel= spsel;
sel{7}.novelDay= 3;
sel{7}.label= 'novel-3';

sel{8}.selectid= [id '-famArm'];
sel{8}.pairsel= pairsel;
sel{8}.spsel= spsel;
sel{8}.label= 'famArm';

%sel{1}.selectid= [id '-novelArm'];
%sel{1}.pairsel= pairsel;
%sel{1}.spsel= 'all';
%sel{1}.novelDay= 1;
%sel{1}.label= 'split';

%sel{2}.selectid= [id '-novelArm'];
%sel{2}.pairsel= pairsel;
%sel{2}.spsel= 'nonripple';
%sel{2}.novelDay= 1;
%sel{2}.label= 'nonripple';

newsel= {}; nS= 0;
for is=1:length(sel)
    if ~isempty(sel{is}); nS= nS+1; newsel{nS}= sel{is}; end
end
sel= newsel;

setRoot;
olddir= pwd;
%cd([root '/work']);
opt.root= root;
nsel=length(sel);


for is=1:nsel
    [XC{is}, AUX{is}]= auxGetCCG(sel{is}, opt);
end


%fprintf(1, 'Analyzing pairs=''%s'', spikes=''%s'', scaling=''%s'', dt= %dms, T=
%ds\n', pairsel, spsel, opt. scaling, round(dt*1000), round(T*1000));
for is=1:nsel
    fprintf(1, 'number of pairs %s: %d used, %d rejected (%d total)\n', ...
            sel{is}.label, AUX{is}.npairs-AUX{is}.nrej, AUX{is}.nrej, AUX{is}.npairs);
end

if showBigMax | showPeakLoc | showSeqComp | showAsymIndex | showArea | plotall | cmpCCG
    load([root '/data/ccg_filter'])
    B= ccg_filter.tf.num;
    A= ccg_filter.tf.den;
%    B= fir1(4, opt.cutfreq*2*opt.dt);
%    A= 1;
    for is=1:nsel
        FXC{is}= filtfilt(B, A, XC{is});
%        FXC{is}= FXC{is}-ones(length(AUX{is}.lags),1)*mean(XC{is});
    %    DT= AUX{is}.lags*opt.dt;
    %    ind= find(DT>=-.070 & DT<=.070);
    %    FXC{is}= FXC{is}(ind,:);
    end
end

if showSeqComp
    % calculated distance between placefields
    global dPF lp
    if rerun
        dPF= {}; lp= {};
        for is=1:nsel
            tmp= auxGetDistancePF(AUX{is}.pairlist, opt);
            dPF{is}= tmp(AUX{is}.ix);
        end
    end
end

if(showSeqComp | showPeakLoc)
    global lp
    if rerun
        % location of center peak 
        [lp]= auxCalcPeakLoc(FXC, opt, AUX, XC);

        for is=1:nsel
            ind= find(isfinite(lp{is}));
            lp{is}= lp{is}(ind);
            if(showSeqComp) dPF{is}= dPF{is}(ind); end
            AUX{is}.ipf= AUX{is}.ix(ind);
            AUX{is}.ip= AUX{is}.pairlist.ip(AUX{is}.ipf);
        end
    end
    save('SCI_data', 'sel', 'AUX', 'lp', 'dPF');
end


if showExcessCorr | showExcessCorrBar
    global ExcessCorr
    if rerun
        ExcessCorr= auxExcessCorr(XC, opt, AUX);
    end
end

if showMax
    global mXC
    if rerun
        maxy= 0; mXC= {};
        for is=1:nsel
            DT= AUX{is}.lags*opt.dt;
            ind= find(DT>=-.070 & DT<=.070);
        %    ind=[1:length(AUX{is}.lags)];
            mtmax= max(XC{is}(ind,:))';
            maxy= max([maxy; mtmax]);
%            mtmin= min(XC{is}(ind,:))';
            mXC{is}= mtmax;
        %    mXC{is}= (mtmax-mtmin);
%             ind= mtmax>1e-6;
%            mXC{is}= (mtmax(ind)-mtmin(ind))./mtmax(ind);
        %    mXC{is}= (mtmax(ind)-mtmin(ind))./(mtmax(ind)+mtmin(ind));
        %    mXC{is}= (mtmax-mtmin)./mtmax;
        %    mXC{is}= mtmin(ind);
        end
    end
end

if showThetaPower
    global thetaPower
    if rerun; thetaPower= auxCalcPS(XC, opt); end
end

%  asymmetry in CCG
if showAsymIndex
    global asym
    if rerun; 
        asym= auxCalcAsymIndex(FXC, opt, AUX); 
        for is=1:nsel
            ind= find(isfinite(asym{is}));
            asym{is}= asym{is}(ind);
        end
    end
end

%  area
if showArea
    global area
    if rerun; area= auxCalcCCGArea(FXC, opt, AUX); end
end

if showBigMax 
    global bigmax ixBM
    if rerun; 
        bigmax= auxGetBigMax(XC, opt, AUX); 
        ixBM= {};
        for is=1:nsel; ixBM{is}= AUX{is}.ix; end
    end
end

%  mean
if showMean
    global m bigmax ixBM
    if rerun; 
        m= auxCalcMean(XC, opt, AUX); 
        for is=1:nsel
            ip= ixBM{is}(find(bigmax{is}<0));
            [c, ia, ix]= intersect(ip, AUX{is}.ix);
            m{is}(ix)= -m{is}(ix);

%            x= dPF{is}; X= [x, ones(length(x),1)];
%            y= lp{is}'; b= regress(y, X); sresid= y- X*b;
%            sind= find(abs(sresid)<40);
%            [c, ix]= intersect(AUX{is}.ix, SCIAUX{is}.ip(sind));
%            m{is}= m{is}(ix);

            ind= find(isfinite(m{is}));
            m{is}= m{is}(ind);
        end
    end
end

%  std
if showStd
    global stdev
    if rerun; 
        stdev= auxCalcStd(XC, opt, AUX); 
        for is=1:nsel
            ind= find(isfinite(stdev{is}));
            stdev{is}= stdev{is}(ind);
        end
    end
end


for iplot=1:length(plots)
    if isempty(plots{iplot}); continue; end
    range= plots{iplot}.range;
    plots{iplot}.ref= opt.ref; plots{iplot}.id= opt.id;

    % regression slope 
    %rs= auxCalcRegSlope(FXC, opt, AUX);
    %for is=1:nsel
    %    tmp= auxGetDistancePF(AUX{is}.pairlist, opt);
    %    tmp= auxGetOverlapPF(AUX{is}.pairlist, opt);
    %    tmp= tmp(AUX{is}.ix);
    %    ind= find(tmp<15);
    %    rs2{is}= rs{is}(ind);
    %end
    %auxCmpDist(rs, sel, 'slope', opt);
    %auxPlotCDF(rs, sel, 'slope', '')

    if(showSeqComp)
        plots{iplot}.xvar= 'distance'; plots{iplot}.xstr= 'distance (cm)';
        plots{iplot}.yvar= 'PeakSep'; plots{iplot}.ystr= 'CCG peak (ms)';
        auxPlotTwoVar({dPF{range}}, {lp{range}}, plots{iplot}, {sel{range}})
    end

    if showExcessCorr
        plots{iplot}.norm= 1; plots{iplot}.size= [2 1.8];
%        plots{iplot}.legend= 1; 
%        plots{iplot}.xrange= [-0.01 0.03];
        plots{iplot}.nbins= 15; 
        plots{iplot}.title= 'excessCorr'; plots{iplot}.label= 'excess corr';
%        plots{iplot}.xrange= [-0.007 0.06];
%        plots{iplot}.yrange=[0 0.7];
        auxCmpDist2({ExcessCorr{range}}, {sel{range}}, plots{iplot});
%        auxTestDist({ExcessCorr{range}}, {sel{range}}, plots{iplot});
    end

    % max of power spectrum
    %[m,f]= auxCalcPSMax(XC, opt);
    %auxCmpDist(f, sel, 'MaxFreq', opt);
    %auxPlotCDF(f, sel, 'MaxFreq', opt.scaling)
    %auxCmpDist(m, sel, 'MaxPower', opt);
    %auxPlotCDF(m, sel, 'MaxPower', opt.scaling)

    if(showMax)
        plots{iplot}.legend= 0; plots{iplot}.size= [1.5 1];
        plots{iplot}.title= 'peakHeight';
        plots{iplot}.label= 'peak height';
        auxCmpDist2({mXC{range}}, {sel{range}}, plots{iplot});
    end

    if showThetaPower
        plots{iplot}.legend= 0; plots{iplot}.size= [1.5 1];
        plots{iplot}.title= 'thetaPower';
        plots{iplot}.label= 'theta power';
        auxCmpDist2({thetaPower{range}}, {sel{range}}, plots{iplot});
    end

    %  asymmetry in CCG
    if showAsymIndex
        plots{iplot}.legend= 0; plots{iplot}.size= [1.5 1];
        plots{iplot}.title= 'asym';
        plots{iplot}.label= 'CCG peak (ms)';
        auxCmpDist2({asym{range}}, {sel{range}}, plots{iplot});
    end

    %  area
    if showArea
        plots{iplot}.legend= 0; plots{iplot}.size= [1.5 1];
        plots{iplot}.title= 'area'; plots{iplot}.label= 'CCG area';
%        plots{iplot}.xrange= [0 400];
        auxCmpDist2({area{range}}, {sel{range}}, plots{iplot});
    end

    %  mean
    if showMean
        plots{iplot}.legend= 1; plots{iplot}.size= [1.5 1];
        plots{iplot}.nbins= 30;
        plots{iplot}.title= 'mean';
        plots{iplot}.label= 'mean CCG (ms)';
%        keyboard
        auxTestDist({m{range}}, {sel{range}}, plots{iplot});
        auxCmpDist2({m{range}}, {sel{range}}, plots{iplot});
    end

    %  std
    if showStd
        plots{iplot}.legend= 0; plots{iplot}.size= [1.5 1];
        plots{iplot}.title= '2nd';
        plots{iplot}.label= 'RMS time lag (ms)';
        auxCmpDist2({stdev{range}}, {sel{range}}, plots{iplot});
    end

    %  peak location
    if showPeakLoc
        plots{iplot}.legend= 0; plots{iplot}.size= [1.5 1];
        plots{iplot}.nbins= 30;
        plots{iplot}.title= 'peakLoc';
        plots{iplot}.label= 'peak CCG (ms)';
        auxTestDist({lp{range}}, {sel{range}}, plots{iplot});
        auxCmpDist2({lp{range}}, {sel{range}}, plots{iplot});
    end

    if showBigMax
        plots{iplot}.legend= 0; plots{iplot}.size= [1.5 1];
        plots{iplot}.nbins= 30;
        plots{iplot}.title= 'bigmax';
        plots{iplot}.label= 'big peak CCG (ms)';
        auxCmpDist2({bigmax{range}}, {sel{range}}, plots{iplot});
    end

    % plot all CCG
    if(cmpCCG)
        opt.cmp= 1;
    end
    if(plotall | cmpCCG)
        opt.plotcol= plots{iplot}.plotcol;
        auxPlotCCG({XC{range}}, {FXC{range}}, {AUX{range}}, opt, {sel{range}});
    end

end

if showExcessCorrBar
    ind= [[2:4]' [5:7]'];
    for d=1:3
        for i=1:2
            j= ind(d,i);
            Zm(d,i)= mean(ExcessCorr{j});
            Zdev(d,i)= std(ExcessCorr{j})/ sqrt(length(ExcessCorr{j}));
        end
        [h,p(d,1)]= kstest2(ExcessCorr{ind(d,1)}, ExcessCorr{ind(d,2)});
%        [h,p(d,1)]= ttest2(ExcessCorr{ind(d,1)}, ExcessCorr{ind(d,2)});
    end
    p

    opt.ylabel= 'excess corr'; opt.title= 'excessCorr'; 
    opt.plotsize=[1.5 2/3]; 
%    opt.ylim= [0 0.7];
    opt.xticklabels= {'day1'  'day2'  'day3'}; 
    opt.plotcol={zeros(3,3), [[1 0 0]; [0 1 0]; [0 0 1]]};
    auxBar(Zm, Zdev, opt, p);
end

function auxCmpDist(x, sel, var, opt)
fprintf(1, '\ncomparing %s\n', var);
nsel= length(x);
for is=1:nsel
    if is==opt.ref | isempty(x{is}); continue; end
    p= ranksum(x{opt.ref}, x{is});
    fprintf(1, 'Ranksum test %s vs. %s p=%.3g.\n', sel{opt.ref}.label, sel{is}.label, p);
    [h,p]= kstest2(x{opt.ref}, x{is});
    fprintf(1, 'KS test %s vs. %s p=%.3g.\n', sel{opt.ref}.label, sel{is}.label, p);
end

function auxHist(x, sel, var, outname)
nsel= length(x);
minx=1e10; maxx= -1e10;
for is=1:nsel; maxx= max([x{is}, maxx]); minx= min([x{is}, minx]); end
b= linspace(minx, maxx, 11);
lstr= {};
for is=1:nsel
    h(:,is)= histc(x{is}', b);
    h(:,is)= h(:,is)/sum(h(:,is));
    lstr{is}= sel{is}.label;
end
figure
%plot(b(1:end-1)+(b(2)-b(1))/2, h(1:end-1,:));
bar(b+(b(2)-b(1))/2,h);
legend(lstr, 0);
xlabel(var);
myprint('large', ['hist_' var '_' outname]);

function auxPlotCDF(x, sel, var, outname)
figure
plotcol= {'k', 'r', 'b', 'g', 'm', 'y', 'c', 'r', 'k'};
lstr= {};
nsel= length(x);
il= 0;
for is=1:nsel
    if isempty(x{is}); continue; end
    h= cdfplot(x{is});
    set(h, 'Color', plotcol{is});
    il=il+1;
    lstr{il}= sel{is}.label;
    hold on
end
legend(lstr, 0);
xlabel(var);
ylabel('cum. fraction');
title('');
figname= [var '_' outname];
set(gcf, 'Name', figname);
myprint(2.5*[1 2/3], figname);

function [lp]= auxCalcPeakLoc(x, opt, aux, rawx)
nsel= length(x);
for is=1:nsel
    M= (length(aux{is}.lags)-1)/2;
    np= size(x{is}, 2);
    for ip=1:np
        lp{is}(ip)= nan; 
%        lp{is}(ip)= 0; 

        y= x{is}(:,ip);
        t= aux{is}.lags*opt.dt;

        ind= find(t>=-0.10 & t<=0.10)';
%        ind= find(t>-0.05 & t<=0.05)';

        if sum(rawx{is}(ind, ip)) < 30 % for seq compresion calc
%        if sum(rawx{is}(ind, ip)) < .2
            continue
        end
        % find peaks
        ipeaks= ind(findPeaks(y(ind)'));
        if isempty(ipeaks); 
            continue
        else
            [m,i]= max(y(ipeaks));
%            fprintf(1,'max= %.2f, mean= %.2f, std= %.2f\n', m, ...
%                mean(y(ind)), std(y(ind)));
%            if m < 5;  continue; end
            if m < mean(y(ind))+1.5*std(y(ind)); continue; end
            lp{is}(ip)= (ipeaks(i)-M-1)*opt.dt*1000;
        end

if(0)
            clf
            hold on
            plot(t, y);
            plot(lp{is}(ip)/1000, m, 'bo');
            set(gca,'xlim', [-.13 .13]);
        end

    end
end
%keyboard


function [pi]= auxExcessCorr(xc, opt, aux);

gstd= 0.25;
nstd=round(gstd/opt.dt);
g= gaussian(nstd, 5*nstd+1);

nstd=round(0.005/opt.dt);
gf= gaussian(nstd, 6*nstd+1);
%b= fir1(200,.001);
%figure(1)
%freqz(b,1);

%    load([root '/data/ccg_theta_low_filter'])
%    load([root '/data/ccg_thetafilter'])
%B= ccg_filter.tf.num;
%A= ccg_filter.tf.den;

for is=1:length(xc)
    dt= aux{is}.lags*opt.dt;
    center= find(aux{is}.lags==0);
    nx= size(xc{is},2);
    for ix=1:nx
        x= xc{is}(:,ix);

%        fx= fxc{is}(:,ix);
%        fx= filtfilt(B,A, x);
        fx= filtfilt(gf, 1, x);

        sx= filtfilt(g,1, x);
%        sx= filtfilt(thetafilter.tf.num, thetafilter.tf.den, x);
%        sx= sx+filtfilt(g,1, x);

%        if aux{is}.ix(ix)== 111;
%        if aux{is}.ix(ix)== 83;
if 0

        figure(1); clf
        maxt= .350;
        hold on
        plot(1000*dt,x, 'Color', .7*[1 1 1])
        plot(1000*dt,fx, 'k')
        plot(1000*dt,sx,'k')
%        maxy= max(x(dt>=-maxt & dt<=maxt));
        maxy= 0.06;
        set(gca, 'xlim', 1000*[-maxt maxt]);
        set(gca, 'ylim', [0 maxy*1.05]);
        plot(1000*dt(center)*[1 1],[0 maxy*1.05], 'k-');
        hold off
        xlabel('time lag (ms)');
        ylabel('cross correl.');
        sprintf(num2str(fx(center)-sx(center)))
        keyboard
        myprint(1.5*[1 2/3], 'bla')
%        pause
end

        pi{is}(ix)= fx(center)-sx(center);
    end
end


function auxCut(x, sel, cut, var)
nsel= length(x);
n= zeros(nsel, 2);
frac= zeros(nsel, 2);
for is=1:nsel
    n(is,1)= sum(x{is}<cut);
    n(is,2)= sum(x{is}>cut);
    frac(is,:)= n(is,:)/sum(n(is,:));
end
frac

% bootstrap analysis
reps= 10000;
ref= 1;
n1= length(x{ref});
f1= zeros(reps, 2);
for r=1:reps
    ind= floor(n1*rand(n1,1))+1;
    xb= x{ref}(ind);
    tmp(1)= sum(xb<cut);
    tmp(2)= sum(xb>=cut);
    f1(r,:)= tmp/n1;
end

for is=1:nsel
    p= sum(f1(:,2)>= frac(is,2))/reps;
    fprintf(1, 'bootstrap test %s vs. %s p=%.3f.\n', sel{1}.label, sel{is}.label, p);
end

function auxPlotCCG(xc, fxc, aux, opt, sel)
if ~opt.cmp
    nsel= length(xc);
else 
    nsel= 1;
end
olddir= pwd;

showMS= 1;
showFXC= 1;

for is=1:nsel
%    workdir= sprintf('%s/work/CCG/%s/%s/%s', ...
%        opt.root, sel{is}.pairsel, sel{is}.spsel, sel{is}.selectid);
%    cd(workdir)
    np= size(xc{is},2);
    pl= aux{is}.pairlist;
    T= aux{is}.lags(end)*opt.dt;
    if showMS
        t= 1000*opt.dt*aux{is}.lags';
        ind= find(t>=-350& t<=350);
    else
        t= opt.dt*aux{is}.lags';
        ind= find(t>=-1& t<=1);
    end
%        ind= 1:length(t);
    nt= length(t);
%    if is~=4; continue; end %@@
    for ix=1:np
        ip= aux{is}.ix(ix);

%        keyboard
        if ~any(ip== [44]); continue; end %@@
%        clf; hold on
        figure; hold on
        if showFXC
%            plot(t(ind), x, 'Color', .7*[1 1 1]);
%            plot(t(ind), fxc{is}(ind,ix), 'Color', opt.plotcol(is,:));
            if ~opt.cmp
                x= fxc{is}(ind,ix);
            else
                x= [fxc{2}(ind,ix) fxc{1}(ind,ix)];
                hline= plot(t(ind), x);
                set(hline(2), 'color', opt.plotcol(is,:));
                set(hline(1), 'color', [.7 .7 .7]);
            end
        else
            x= xc{is}(ind,ix);
        end

        %@@
%        fx= t(ind); fi= find(fx>=-100 & fx <=100); fx= fx(fi);
%        fy= x(fi);
%        h= fill([flipud(fx);fx], [.1*ones(length(fi),1);fy], .7*[1 1 1]);
%        set(h,'EdgeColor', .7*[1 1 1]);

        if ~opt.cmp
            plot(t(ind), x, 'Color', opt.plotcol(is,:));
        end
        set(gca, 'XLim', [t(ind(1)), t(ind(end))]);
        maxx= max(max(x));
        if maxx==1; maxy= 0.2; else maxy= (maxx+1e-10)*1.05; end
        set(gca, 'YLim', [0 maxy]);
        plot([0 0], [0 max(1, maxy)], 'k-');
%        plot([0 0], [0 25], 'k-'); %@@

        num= pl.cellnum(ip,:);
        tstr= sprintf('%s: %d %d %d %d', pl.rat{ip}, num(1:4));
        fnstr= sprintf('%s-%d-%s-%d-%d-%d-%d', sel{is}.label, ip, pl.rat{ip}, num(1:4));
        if(length(num)>4)
            tstr= sprintf('%s, %d %d', tstr, num(7:8));
            fnstr= sprintf('%s_%d-%d', fnstr, num(7:8));
            ylabel('cross correl.');
        else
            ylabel('auto correl.');
        end
        fnstr= sprintf('%s_%s', fnstr, sel{is}.spsel);
        if isfield(pl, 'traj')
            tstr= [tstr '; traj '];
            trajs=  pl.traj{ip};
            for itraj=1:length(trajs)
                tstr= sprintf('%s %d', tstr, trajs(itraj));
            end
        end
%        title(tstr);
        if showMS
            xlabel('time lag (ms)');
        else
            xlabel('time lag (s)');
        end
        set(gcf, 'Name', fnstr);
%        myprint('small', fnstr, [], 1, []);
        myprint('small', fnstr);
    end
end
cd(olddir)

function [m,f]= auxCalcPSMax(x, opt)
nsel= length(x);
for is=1:nsel
    np= size(x{is}, 2);
    for ip=1:np
        y= x{is}(:,ip);
%        y= y-mean(y);
        % normalize
%        if mean(y.^2) > 0
%            y= y/sqrt(mean(y.^2));
%        end
%        if mean(y) > 0
%            y= y/mean(y);
%        end
%        [Px,freq] = pwelch(y,ip),[],[],[],1/opt.dt);
        [Px, freq] = pmtm(y,[],[],1/opt.dt,'eigen');
        % find peaks
        i= findPeaks(Px');
        [tmp, imin]= min(abs(freq(i)-8));
        m{is}(ip)= Px(i(imin));
        f{is}(ip)= freq(i(imin));
    end
end

function m= auxCalcPS(x, opt)
nsel= length(x);
for is=1:nsel
    disp(is);
    np= size(x{is}, 2);
%    for ip=[27,35]  % 27, 35
    for ip=1:np
        y= x{is}(:,ip);
%        if sum(y) > 1e-2
%            y= y/mean(y);
%        end
%        y= y-mean(y);
        % normalize
%        if mean(y.^2) > 0
%            y= y/sqrt(mean(y.^2));
%        end
%        [Px,freq] = pwelch(x{is}(:,ip),[],[],[],1/opt.dt);
        [Px, freq] = pmtm(y,[],[],1/opt.dt,'eigen');

%        hold off
%        semilogy(freq, Px);
%        set(gca, 'XLim', [0,100]);
%        i= findPeaks(Px');
%        hold on
%        semilogy(freq(i), Px(i), '.');
%        pause

        % normalize
%        ind= find(1<freq & freq<100);
%        if sum(Px(ind))>0;
%            Px= Px/sum(Px(ind));
%        end

        % get power in theta band
%        ind= find(freq<6);
        ind= find(6<freq & freq<10);
%        ind= find(11<freq & freq<100);
        m{is}(ip)= sum(Px(ind));
    end
end


function d= auxGetDistancePF(pl,opt)
np= length(pl.rat);
d= zeros(np,1);
traj_fac=-[1 -1 1 -1];

if(0)
    for ip=1:np
        d(ip)= mean(pl.pf{ip}(3:4))- mean(pl.pf{ip}(1:2));
        d(ip)= traj_fac(pl.traj{ip}+1)*d(ip);
    end
else 
    olddir= pwd;
    ind=[1:4; 5:8];
    indPF=[1:2; 3:4];
    for ip=1:np
        if ip==1 | ~strcmp(pl.rat{ip}, pl.rat{ip-1})
            cd(sprintf('%s/%s/data2/', opt.root, pl.rat{ip}));
            load rate
            setLocalOptions
        end
        for j=1:2
            num= pl.cellnum(ip,ind(j,:));
            rdata= rate{num(1)}{num(2)}{num(3)}{num(4)}.data;
            traj= pl.traj{ip};
            x= pl.pf{ip}(indPF(j,1)):1:pl.pf{ip}(indPF(j,2));
            r= interp1(rdata(:,1), rdata(:,traj+2), x, 'linear');
            cent(j)= sum(r.*x)/sum(r);
        end
        d(ip)= traj_fac(pl.traj{ip}+1)*diff(cent);
    end
    cd(olddir);
end


function auxPlotTwoVar(x,y,opt,sel)
nx= length(x);
if nx~= length(y); 
    error('  The two variables have different number of sets.');
end
plotcol= opt.plotcol;
for i=1:nx
    if ~strcmp(opt.plot, 'none'); figure; hold on; end
    fprintf(1, '%s, \t', sel{i}.label);
    opt.plotcol= plotcol(i,:);
    opt.xrange= [-40 40]; opt.yrange= [-105 105];
    [p(i), r(i), n(i)]= linrel(x{i}, y{i}, opt, 0.05);
    if ~strcmp(opt.plot, 'none'); 
        figname= ['seq_comp_' opt.id '_' sel{i}.label];
        set(gcf, 'Name', figname);
        myprint(1.5*[1 2/3], figname);
    end
%    myprint(2*[1 2/3], figname);
end
%myprint('small', [opt.yvar '-vs-' opt.xvar]);

z= 0.5*log((1+r)./(1-r));

if n(opt.ref)<4; error('need more than 3 data points for testing corr'); end
for i=1:nx
    if i== opt.ref; continue; end
    if n(i)<4; warning('skipping test (less than 4 datapoints)'); end
    Z= (z(i)-z(opt.ref))/sqrt(1/(n(i)-3)+1/(n(opt.ref)-3));
    if Z>0; Z=-Z; end
    pval= normcdf(Z)*2;
    fprintf(1, '  "%s" vs. "%s": p=%.3f\n',  sel{opt.ref}.label, sel{i}.label, pval);
end


function d= auxGetOverlapPF(pl,opt)
np= length(pl.rat);
d= zeros(np,1);
for ip=1:np
    x= reshape(pl.pf{ip},2,2);
    ol= [max(x(1,:)) min(x(2,:))];
%    d(ip)= 2*diff(ol)/(sum(diff(x)));
    d(ip)= diff(ol);

%    d(ip)= abs(diff(mean(x))); % distance between centers
%    d(ip)= diff(x(1,:)); % distance between centers

%    d(ip)= abs(diff(mean(x))) / mean(diff(x)); % distance between centers, norm
%    by width
end

function [rs, Rsq]= auxCalcRegSlope(x, opt, aux)
nsel= length(x);
for is=1:nsel
    t= aux{is}.lags'*opt.dt;
    nt=length(t);
    np= size(x{is}, 2);
    for ip=1:np
        [b,bint,r,rint,stats]= regress(x{is}(:,ip), [t ones(nt,1)]);
%        fprintf(1,'%.4g\n', stats(3));
%        rs{is}(ip,:)= abs(b(1));
        rs{is}(ip,:)= b(1);
        Rsq{is}(ip,:)= stats(1);

%        plot(t,x{is}(:,ip))
%        hold on
%        plot(t,b(1)*t+b(2), 'k');
%        hold off
%        pause
    end
end

%function [rs, Rsq]= auxCalcRegSlope(x, opt, aux)
%figure
%nsel= length(x);
%for is=1:nsel
%    t= aux{is}.lags'*opt.dt;
%    nt=length(t);
%    np= size(x{is}, 2);
%    for ip=1:np
%        y= x{is}(:,ip);
%        i=findPeaks(-y');
%        n=length(i);
%        [b,bint,r,rint,stats]= regress(y(i), [t(i) ones(n,1)]);
%        rs{is}(ip,:)= abs(b(1));
%        rs{is}(ip,:)= b(1);
%        Rsq{is}(ip,:)= stats(1);

%        plot(t,y);
%        hold on
%        plot(t,b(1)*t+b(2), 'k');
%        plot(t(i), y(i), 'ro');
%        hold off
%        pause
%    end
%end

function area= auxCalcCCGArea(x, opt, aux)
nsel= length(x);
for is=1:nsel
    np= size(x{is}, 2);
    dt= aux{is}.lags*opt.dt;
    ind= find(dt>-0.01 & dt<0.01);
    for ip=1:np
        area{is}(ip,:)= sum(x{is}(ind,ip));
    end
end

function asym= auxCalcAsymIndex(x, opt, aux)
nsel= length(x);
for is=1:nsel
    np= size(x{is}, 2);
    for ip=1:np
        y= x{is}(:,ip);
        t= aux{is}.lags*opt.dt;
        ind= find(t>-0.1 & t<=0);
        r1= sum(y(ind));

        ind= find(t>0 & t<=0.1);
        r2= sum(y(ind));

        if r1+r2 > 0; 
            asym{is}(ip,:)= (r2-r1)/(r1+r2);
        else
            asym{is}(ip,:)= nan;
        end
    end
end
%keyboard

function m= auxCalcMean(x, opt, aux)
nsel= length(x);
for is=1:nsel
    np= size(x{is}, 2);
    for ip=1:np
        y= x{is}(:,ip);
        t= aux{is}.lags'*opt.dt;
        ind= find(t>=-0.1 & t<=0.1);
        den= sum(y(ind));
        if den >= 10
            yn= y(ind)/ den;
            m{is}(ip,:)= 1000*sum(yn.*t(ind));
        else 
            m{is}(ip,:)= nan;
        end
    end
end

function stdev= auxCalcStd(x, opt, aux)
nsel= length(x);
for is=1:nsel
    nx= size(x{is}, 2);
    for ix=1:nx
        y= x{is}(:,ix);
        t= 1000*aux{is}.lags'*opt.dt;
        ind= find(t>=-100 & t<= 100);
        den= sum(y(ind));
        if den > 10
            yn= y(ind)/ den;
    % analyze power spectrum : peak height, width, peak location
%    [ix, lx, hy, wx]= anaPeaks(y(ind),t(ind), .50, 'both');
%            [mh, i]= max(hy);
%            stdev{is}(ix,:)= lx(i);
%            m= sum(yn.*t(ind));
%            stdev{is}(ix,:)= sqrt(sum((t(ind)-m).^2.*yn));
            stdev{is}(ix,:)= sqrt(sum((t(ind)).^2.*yn));
        else 
            stdev{is}(ix,:)= nan;
        end

        %@@
%        ip= aux{is}.ix(ix);
%        if ip== 28 | ip== 50; fprintf(1, 'is= %d, ip= %d, rms= %f\n', is, ip, stdev{is}(ix,:)); end
    end
end

function m= auxGetBigMax(x, opt, aux)
nsel= length(x);
for is=1:nsel
    np= size(x{is}, 2);
    for ip=1:np
        [mtmp, i]= max(x{is}(:,ip));
        t= aux{is}.lags'*opt.dt;
        m{is}(ip,:)= 1000*sum(t(i));
    end
end
%keyboard

