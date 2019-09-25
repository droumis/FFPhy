function showSync
%function showSync
%  Summarize  cross correlation result showing change in synchronous firing.

rerun= 0;

id= 'pf9';
pairsel= 'overlap';
spsel= {'placefield', 'placefield_nonripple', 'nonripple', 'all', 'nonripple_nonrun', 'nonrun', 'ripple'};
spsellabels= {'PF', 'all-rip', 'all', 'rest-rip', 'rest', 'ripples'};
selectid= {[id '-famArm'], [id '-novelArm'], [id '-fam']};
selectlabel= {'fam arm', 'novel arm', 'fam conf'};

setRoot;
fname= [root '/work/CCG/' pairsel '/synch-' id '.mat'];

if rerun | ~exist(fname, 'file')
    opt.scaling= 'coef';
    opt.dt= 0.002;      % [sec] timestep
    opt.T= 2;         % [sec] +/- width of CCG window
    opt.cutfreq= 40;    % cutoff frequency [Hz] 
    opt.minspikes= 0;
    opt.ref= 1;


    nsp= length(spsel);
    nsel= length(selectid);

    for p=1:nsp
        for s=1:nsel
            for d=1:3
                sel{p}{s}{d}.selectid= selectid{s};
                sel{p}{s}{d}.spsel= spsel{p};
                sel{p}{s}{d}.pairsel= pairsel;
                if strfind(selectid{s}, 'Arm');
                    sel{p}{s}{d}.novelDay= d;
                end
            end
        end
    end

    opt.root= root;

    load([root '/work/ccg_filter']);
    B= ccg_filter.tf.num;
    A= ccg_filter.tf.den;
    %B= fir1(4, opt.cutfreq*2*opt.dt);
    %A= 1;
    for p=1:nsp
        for s=1:nsel
            for d=1:3
                [XC{p}{s}{d}, AUX{p}{s}{d}]= auxGetCCG(sel{p}{s}{d}, opt);
                FXC{p}{s}{d}= filtfilt(B, A, XC{p}{s}{d});
            end
        end
    end
    save(fname, 'opt', 'sel', 'XC', 'FXC', 'AUX');
else % rerun
    load(fname);
    nsp= length(sel);
    nsel= length(sel{1});
end % rerun


if(0)
    % find peak ratios
pr= {};
for p=1:nsp
    for s=1:nsel
        for d=1:3
            pr{p}{s}{d}= auxCalcPeakRatios(FXC{p}{s}{d}, opt.dt, ...
                    AUX{p}{s}{d}.lags);
        end
    end
end

auxPlot(pr, 'peak ratio', 'peakRatio', spsellabels, selectlabel, pairsel, id);
end

if(0)
    % find maximum correlation
    peak= {};
    for p=1:nsp
        for s=1:nsel
            for d=1:3
                DT= AUX{p}{s}{d}.lags*opt.dt;
                ind= find(DT>=-.050 & DT<=.50);
%                ind= find(DT>=.050 & DT<=.180);
                peak{p}{s}{d}= max(FXC{p}{s}{d}(ind,:))';
            end
        end
    end
    auxPlot(peak, 'peak', 'peak', spsellabels, selectlabel, pairsel, id);
end

if(1)
%    controlVar= {{'spikes/'; 'ripple'}, {'burst/';'ripple'},...
%        {'spikes/';'burst'}, 'number', {'length'; '(ms)'}, {'frac';'time'}};
    controlVar= {{'spikes/'; 'ripple'}, {'burst/';'ripple'}, {'spikes/';'burst'}};
    for s=1:nsel
        for d=1:3
            [tmp, tmpaux]= auxAnalyzeSpikes(sel{7}{s}{d}, opt);

            control{1}{s}{d}= tmp.Nspikes./ tmp.Nevents(tmpaux.cl2ep);
            control{2}{s}{d}= tmp.Nbursts./ tmp.Nevents(tmpaux.cl2ep);
            control{3}{s}{d}= control{1}{s}{d}./control{2}{s}{d};
            control{4}{s}{d}= tmp.Nevents;
            control{5}{s}{d}= tmp.T*1000./tmp.Nevents;
            control{6}{s}{d}= tmp.fracT;
        end
    end
    auxPlot(control, [], 'rippleControl', controlVar, selectlabel, pairsel, id);

    fprintf(1, '\n(ranksum test, KS test, t-test)\n');
    for v=1:length(controlVar)
        if iscell(controlVar{v})
            fprintf(1, '%s\t', strcat(controlVar{v}{:}));
        else
            fprintf(1, '%s\t', strcat(controlVar{v}));
        end
        for d=1:3
            prank= ranksum(control{v}{2}{d},control{v}{3}{d});
            [h,pKS]= kstest2(control{v}{2}{d},control{v}{3}{d});
            [h,pT]= ttest2(control{v}{2}{d},control{v}{3}{d});
            fprintf(1, 'day%d: (%.2g, %.2g, %.2g)\t', ...  
                d, prank, pKS, pT);
        end
        fprintf(1, '\n');
    end

end


function auxPlot(pr, varname, figlabel, labels1, labels2, pairsel, id)

%plotcol= {'k', 'r', 'c'};
plotcol= {'k', 'r', 0.8*[1 1 1];...
    'k', 'g', 0.8*[1 1 1];...
    'k', 'b', 0.8*[1 1 1]};
%plotcol= {'k', 'r', 0.8*[1 1 1]};
if isempty(varname); 
    ncol= length(labels1);
%    plotcol= {'r', 'g', 'b'};
else
    ncol= 1;
%    plotcol= {'k', 'r', 0.8*[1 1 1]};
end
nsp= length(pr);
nsel= length(pr{1});
m= zeros(nsp,nsel,3);
n= zeros(nsp,nsel,3);
sd= zeros(nsp,nsel,3);
dev= zeros(nsp,nsel,3);
for p=1:nsp
    for s=1:nsel
        for d=1:3
            sdata= sort(pr{p}{s}{d});
%            m(p,s,d)= median(sdata);
            m(p,s,d)= mean(sdata);
            N= length(sdata);
            n(p,s,d)= length(sdata);
            sd(p,s,d)= std(sdata);
%            se(p,s,d)= sd(p,s,d)/sqrt(n(p,s,d));
%            dev(p,s,d)= (sdata(ceil(N*0.90))-sdata(ceil(N*0.1)))/2;
            dev(p,s,d)= sd(p,s,d)/sqrt(n(p,s,d));
        end
    end
end


maxy= nan;
for d=1:3
%    maxy= max([maxy, max(max(m(1:5,:,d)))]);
    maxy= max([maxy, max(max(m(:,:,d)+dev(:,:,d)))]);
end

figure
figname= [figlabel '-' id '-' pairsel];
set(gcf, 'Name', figname);

if ncol==1
    for d=1:3
        subplot(3,1,d);
        hold on
        h= bar(m(:,:,d));
        set(gca, 'YLim', [0 1.05*maxy]);
        for s=1:nsel
            set(h(s), 'FaceColor', plotcol{d,s});
%            xd= get(get(h(s), 'Children'), 'XData');
            xd= get(h(s), 'XData');
            errorbar(mean(xd(2:3,:)), m(:,s,d), dev(:,s,d), '.', 'Color', 0.4*[1 1 1]);
    %        errorbar(mean(xd(2:3,:)), m(:,s,d), se(:,s,d));
        end
        ylabel(varname);
        if d==3; 
            set(gca, 'XTick', [1:nsp])
            set(gca, 'XTickLabel', labels1)
        else
            set(gca, 'XTickLabel', {})
        end
    %    ht= text(.5, 0.9*maxy, ['day ' num2str(d)]);
    %    set(ht, 'VerticalAlignment', 'top');
        hold off
    end
    myprint('small', figname);
else
%    for d=1:3
%        for v=1:ncol
%            maxy= max(max(squeeze(m(v,:,:))));
%            subplot(3,ncol,(d-1)*ncol+v);
%            hold on
%            h= bar(m(v,:,d));
%            set(gca, 'YLim', [0 1.05*maxy]);
%            ylabel(labels1{v})
%            set(gca, 'XTickLabel', {})
%        end
%    end
    for v=1:ncol
        if isempty(varname); 
            maxy= max(max(squeeze(m(v,:,:)+dev(v,:,:))));
        end
        subplot(ncol,1,v);
        hold on
        h= bar(squeeze(m(v,:,:))');
        for s=1:nsel
%            tmp= get(get(h(s), 'Children'), 'XData');
            tmp= get(h(s), 'XData');
            xd{s}= mean(tmp(2:3,:));
        end
        bw= diff(tmp(2:3,1));
        cla
        for d=1:3
            for s=1:nsel
                hs= bar(xd{s}(d)+[0 1],[m(v,s,d) 0], bw);
                set(hs, 'FaceColor', plotcol{d,s});
                errorbar(xd{s}(d), m(v,s,d), dev(v,s,d), 'Color', 0.4*[1 1 1]);
            end
        end
        set(gca, 'YLim', [0 1.05*maxy]);
        set(gca, 'XLim', [0.5 3.5]);
        ylabel(labels1{v})
        if v==ncol
            set(gca, 'XTick', [1:nsel])
            set(gca, 'XTickLabel', {'day 1', 'day 2', 'day 3'})
%            set(gca, 'XTickLabel', labels2)
        else
            set(gca, 'XTickLabel', {})
        end
    end
    myprint('mini', figname);
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
np= size(x, 2);
pr= [];
for ip=1:np
    % find peaks
    for j=1:3
        [m(j), l(j)]= max(x(ind{j},ip));
    end
    if all(m<10e-4)
        pr(ip)= nan;
    elseif m(1)+m(3)<10e-10
%                pr(ip)= nan;
            pr(ip)= 10;
    else
        pr(ip)= 2*m(2)/(m(1)+m(3));
    end
end
%    pr(isnan(pr))= max(pr);
%    pr(isnan(pr))= 10;

pr= pr(isfinite(pr));


function [var, aux]=auxAnalyzeSpikes(sel, opt)

workdir= sprintf('%s/work/CCG/misc', opt.root);
fname= sprintf('%s/beh_spikes_%s_%s.mat', workdir, sel.selectid, sel.spsel);

if ~exist(fname, 'file')
%if 1
%    test whether directory exists
    if ~exist(workdir); error(['directory ' workdir ' does not exist']); end
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
        if ic==1 | ~strcmp(rat, cl.rat{ic-1}) | d~= cl.cellnum(ic-1,1)
            load(sprintf('%s/%s/data2/spikedata%.2d.mat', opt.root, rat, d));
            load(sprintf('%s/%s/data2/behavdata%.2d.mat', opt.root, rat, d));
        end
        if ic==1 | ~strcmp(rat, cl.rat{ic-1}) ...
            | d~= cl.cellnum(ic-1,1)| e~= cl.cellnum(ic-1,2)
            bd= behavdata{d}{e};
            ie= ie+1;
            switch sel.spsel
            case 'all'
                ind= [1:length(bd.time)];
            case 'ripple'
                ind= find(bd.ripple);
            case 'nonripple'
                ind= find(~bd.ripple);
            end
            var.fracT(ie)= length(ind)/length(bd.time);
            var.T(ie)= length(ind)*0.002;
            var.Nevents(ie)= sum(diff(ind)>1)+1;
        end

        sd= spikedata{d}{e}{tet}{c};

        opt.minspikes=0;
        [t, rej]= auxSelectSpikes({sd}, bd, sel.spsel, [], [], opt);
        var.Nspikes(ic)= length(t{1});
        var.Nbursts(ic)= sum(diff(t{1})>0.010)+1;
        var.Rate(ic)= var.Nspikes(ic)/var.T(ie);
        aux.cl2ep(ic)= ie;
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
