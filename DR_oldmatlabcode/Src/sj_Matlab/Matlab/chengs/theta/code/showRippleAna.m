function showRippleAna
%function showRippleAna

plotall= 1;
pairsel= 'overlap';  % 'none', 'overlap', 'all'
spsel= 'ripple';   % spike selection {'all', 'traj', 'sat', 'traj_sat', 'run'}
%                 'traj_sat2', 'placefield', 'ripple', 'non_run'
opt.scaling= 'coef';
opt.dt= 0.002;      % [sec] timestep
opt.T= 0.200;         % [sec] +/- width of CCG window
opt.cutfreq= 40;    % cutoff frequency [Hz] 
opt.minspikes= 40;  % minimum number of spikes to include
opt.ref= 1;


sel{1}.selectid= 'pf9-famArm';
sel{1}.pairsel= pairsel;
sel{1}.spsel= spsel;
sel{1}.label= 'famArm';

sel{2}.selectid= 'pf9-novelArm';
sel{2}.pairsel= pairsel;
sel{2}.spsel= spsel;
sel{2}.novelDay= 1;
sel{2}.label= 'day-1';

sel{3}.selectid= 'pf9-novelArm';
sel{3}.pairsel= pairsel;
sel{3}.spsel= spsel;
sel{3}.novelDay= 2;
sel{3}.label= 'day-2';

sel{4}.selectid= 'pf9-novelArm';
sel{4}.pairsel= pairsel;
sel{4}.spsel= spsel;
sel{4}.novelDay= 3;
sel{4}.label= 'day-3';

sel{5}.selectid= 'pf9-novelArm';
sel{5}.pairsel= 'non_overlap';
sel{5}.spsel= spsel;
sel{5}.label= 'novel non-overlap';

setRoot;
opt.root= root;
nsel=length(sel);

acc.length= [];
acc.frac= [];
for is=1:nsel
    [VAR{is}, AUX{is}]= auxAnalyzeSpikes(sel{is}, opt);
%    RATE{is}= VAR{is}.Rate;
    SPIKES_BURST{is}= VAR{is}.Nspikes./VAR{is}.Nbursts;
    SPIKES_EVENT{is}= VAR{is}.Nspikes./VAR{is}.Nevents(AUX{is}.cl2ep);
    BURSTS_EVENT{is}= VAR{is}.Nbursts./VAR{is}.Nevents(AUX{is}.cl2ep);
    NSPIKES{is}= VAR{is}.Nspikes;
    NEVENTS{is}= VAR{is}.Nevents;
    T{is}= VAR{is}.T;
    LENGTH{is}= VAR{is}.T./NEVENTS{is};
    FRAC{is}= VAR{is}.fracT;

    acc.length= [acc.length LENGTH{is}];
    acc.frac= [acc.frac FRAC{is}];
end

% analyze behavior
%auxCmpDist(T, sel, [sel{1}.spsel 'Time'],opt);
%auxPlotCDF(T, sel, [sel{1}.spsel 'Time'], '')
%auxCmpDist(NEVENTS, sel, [sel{1}.spsel 'Events'],opt);
%auxPlotCDF(NEVENTS, sel, [sel{1}.spsel 'Events'], '')
%auxCmpDist(SPIKES_BURST, sel, 'SpikesPerBurst', opt)
%auxPlotCDF(SPIKES_BURST, sel, 'SpikesPerBurst', '')
%auxCmpDist(SPIKES_EVENT, sel, 'SpikesPerEvent', opt)
%auxPlotCDF(SPIKES_EVENT, sel, 'SpikesPerEvent', '')
%auxCmpDist(BURSTS_EVENT, sel, 'BurstsPerEvent', opt)
%auxPlotCDF(BURSTS_EVENT, sel, 'BurstsPerEvent', '')
%auxCmpDist(LENGTH, sel, [sel{1}.spsel 'Length'], opt)
%keyboard
opt.title= 'ripple length';
opt.outname= [sel{1}.spsel 'Length'];
opt.label= 'Length'; opt.plot= 'pdf';
auxCmpDist2(LENGTH, sel, opt)
%auxPlotCDF(LENGTH, sel, [sel{1}.spsel 'Length'], '')
fprintf(1, 'ripple length= %.2f +/- %.2f ms\n', 1000*mean(acc.length), 1000*std(acc.length));
auxCmpDist(FRAC, sel, [sel{1}.spsel 'Frac'], opt)
auxPlotCDF(FRAC, sel, [sel{1}.spsel 'Frac'], '')
fprintf(1, 'ripple frac= %.2f +/- %.2f%%\n', 100*mean(acc.frac), 100*std(acc.frac));


function auxCmpDist(x, sel, var, opt)
fprintf(1, '\ncomparing %s\n', var);
nsel= length(x);
for is=1:nsel
    if is==opt.ref; continue; end
    p= ranksum(x{opt.ref}, x{is});
    fprintf(1, 'Ranksum test %s vs. %s p=%.3f.\n', sel{opt.ref}.label, sel{is}.label, p);
    [h,p]= kstest2(x{opt.ref}, x{is});
    fprintf(1, 'KS test %s vs. %s p=%.3f.\n', sel{opt.ref}.label, sel{is}.label, p);
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
plotcol= {'k', 'r', 'b', 'g', 'm', 'y', 'c'};
lstr= {};
nsel= length(x);
for is=1:nsel
    h= cdfplot(x{is});
    set(h, 'Color', plotcol{is});
    lstr{is}= sel{is}.label;
    hold on
end
%legend(lstr, 0);
xlabel(var);
ylabel('cdf');
title('');
myprint('small', ['cdf_' var '_' outname]);

function lp= auxCalcPeakLoc(x, dt, lags)
M= (length(lags)-1)/2;
nsel= length(x);
for is=1:nsel
    np= size(x{is}, 2);
    for ip=1:np
        % find peaks
        ipeaks= findpeaks(x{is}(:,ip)');
        i= min(abs(ipeaks-M-1));
        lp{is}(ip)= i*dt;
    end
end
%keyboard

function pr= auxCalcPeakRatios(x, dt, lags)
M= (length(lags)-1)/2;
% interval of preceding, central, and following theta cycle
tb=[-.1875, -.0625, .0625, .1875];
ti= floor(tb/dt)+M+1;
for j=1:3
    ind{j}= [ti(j):ti(j+1)];
end
nsel= length(x);
for is=1:nsel
    np= size(x{is}, 2);
    for ip=1:np
        % find peaks
        for j=1:3
            [m(j), l(j)]= max(x{is}(ind{j},ip));
        end
        if all(m<10e-4)
            pr{is}(ip)= nan;
        elseif m(1)+m(3)<10e-10
%                pr{is}(ip)= nan;
                pr{is}(ip)= 10;
        else
            pr{is}(ip)= 2*m(2)/(m(1)+m(3));
        end
    end
%    pr{is}(isnan(pr{is}))= max(pr{is});
%    pr{is}(isnan(pr{is}))= 10;
    
    pr{is}= pr{is}(isfinite(pr{is}));
end


function pr= auxCalcPeakIndex(x, dt, lags)
M= (length(lags)-1)/2;
% interval of preceding, central, and following theta cycle
tb=[-.1875, -.0625, .0625, .1875];
ti= floor(tb/dt)+M+1;
for j=1:3
    ind{j}= [ti(j):ti(j+1)];
end
nsel= length(x);
for is=1:nsel
    np= size(x{is}, 2);
    for ip=1:np
        % find peaks
        for j=1:3
            [m(j), l(j)]= max(x{is}(ind{j},ip));
        end
        av= (m(1)+m(3))/2;
        pr{is}(ip)= (m(2)-av)/(m(2)+av);
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

function auxPlotCCG(xc, aux, opt, sel)
figure
nsel= length(xc);
olddir= pwd;

for is=1:nsel
    workdir= sprintf('%s/work/CCG/%s/%s/%s', ...
        opt.root, sel{is}.pairsel, sel{is}.spsel, sel{is}.selectid);
    cd(workdir)
    np= size(xc{is},2);
    pl= aux{is}.pairlist;
    T= aux{is}.lags(end)*opt.dt;
    for ix=1:np
        ip= aux{is}.ix(ix);

%        if ~((is==1 & ip==12) | (is==2 & ip==5)); continue; end
        plot(opt.dt*aux{is}.lags, xc{is}(:,ix));
        set(gca, 'XLim', [-T, T]);
        num= pl.cellnum(ip,:);
        tstr= sprintf('%s: %d %d %d %d', pl.rat{ip}, num(1:4));
        if(length(num)>4)
            tstr= sprintf('%s, %d %d', tstr, num(7:8));
        end
        if isfield(pl, 'traj')
            tstr= [tstr '; traj '];
            trajs=  pl.traj{ip};
            for itraj=1:length(trajs)
                tstr= sprintf('%s %d', tstr, trajs(itraj));
            end
        end
        title(tstr);
        xlabel('time lag (sec)');
        ylabel('cross correl');
        myprint('large', sprintf('CCG-%s-%d-%s-%d-%d-%d-%d__%d-%d', ...
            sel{is}.label, ip, pl.rat{ip}, num([1:4,7,8])));
    end
end
cd(olddir)

function [var, aux]=auxAnalyzeSpikes(sel, opt)

global behavdata
workdir= sprintf('%s/work/CCG/misc', opt.root);
fname= sprintf('%s/beh_spikes_%s_%s.mat', workdir, sel.selectid, sel.spsel);

%if ~exist(fname, 'file')
if 1
%    test whether directory exists
    if ~exist(workdir); error(['directory ' workdir ' does not exist']); end
    nsel= length(sel);
    fprintf(1, 'Calculating cell select ''%s'', spike select= ''%s''...\n', ...
        sel.selectid, sel.spsel);
    [cl, nc]= collectCellList(sel.selectid);
    aux.celllist= cl;
    aux.ncells= nc;
    ie=0;
    for ic=1:nc
        fprintf(1, '  %3d/%3d\n', ic, nc);
        rat= cl.rat{ic};
        num= cl.cellnum(ic,:); d=num(1); e=num(2); tet=num(3); c=num(4);
        setRoot
        data2dir= fullfile(root,rat,'data2');
        loadVar(data2dir, 'behavdata', d);
        if ic==1 | ~strcmp(rat, cl.rat{ic-1}) ...
            | d~= cl.cellnum(ic-1,1)| e~= cl.cellnum(ic-1,2)
            bd= behavdata{d}{e};
            ie= ie+1;
            switch sel.spsel
            case 'all'
                ind= [1:length(bd.time)];
            case 'ripple'
                ind= find(bd.ripple);
            case 'non_ripple'
                ind= find(~bd.ripple);
            end
            var.fracT(ie)= length(ind)/length(bd.time);
            var.T(ie)= length(ind)*0.002;
            var.Nevents(ie)= sum(diff(ind)>1)+1;
            [lo,hi]= findcontiguous(ind);
            var.rip_length{ie}= (hi-lo+1)*0.002;
        end
    end
    save(fname, 'var', 'aux');
else
    load(fname);
end
if isfield(sel, 'novelDay') & sel.novelDay
     ind= find(aux.celllist.day==sel.novelDay);
    aux.ncells= length(ind);
else
    ind=[1:aux.ncells];
end

% clean-up
aux.ix= ind;
var.Nspikes= var.Nspikes(ind);
var.Nbursts= var.Nbursts(ind);
var.Rate= var.Rate(ind);
[ie,i,j]= unique(aux.cl2ep(ind));
var.fracT= var.fracT(ie);
var.T= var.T(ie);
var.Nevents= var.Nevents(ie);
aux.cl2ep= j;
%keyboard
