function  cmpStats 
%function  cmpStats 
%
% Multi purpose function to visualize statistics calculated on adaptive
% filtering results.

global RATSEL
RATSEL= [];  % '', 'kyl', 'ter', 'sta', 'fel'

rerun= 0;
plotsingle= 0;
plotVsSets= 0;
saveData= 0;
plotMultiMean= 0;
plotMultiHist= 0;
plotMultiDist= 0;
plotDist= 0;
plotmean= 0;
plotSessionMean= 1;
plotPDF= 0;
derivative= 0;
plotcorr= 0;
plotpair= 0;
plotscatter= 0;

%timeind= [60 110 300];
timeind= [1 180];
global plotrange plotcol fillcol
plotcol= hsv(20);
plotcol(1:6,:)=[[0 0 0]; [1 0 0]; [0 0 1]; [0 1 0]; [.7 0 1]; [0 1 1]];

id= 'pf9';

%plots{1}.range=[2:4, 5:7]; plots{1}.outname= [id];
%plots{1}.plotcol={zeros(3,3), [[1 0 0]; [0 1 0]; [0 0 1]]};

%plots{2}.range= [2,5]; plots{2}.outname= [id '_day1'];
%plots{2}.plotcol=[[0 0 0]; [1 0 0]];
%plots{3}.range= [3,6]; plots{3}.outname= [id '_day2'];
%plots{3}.plotcol=[[0 0 0]; [0 1 0]];
%plots{4}.range= [4,7]; plots{4}.outname= [id '_day3'];
%plots{4}.plotcol=[[0 0 0]; [0 0 1]];


%plots{5}.range= [1,2,5]; plots{5}.outname= [id '_fam1'];
%plots{5}.plotcol=[.7*[1 1 1]; [0 0 0]; [1 0 0]];
plots{6}.range= [1,2:4]; plots{6}.outname= [id '_famArm'];
plots{6}.plotcol={[.7*[1 1 1]; [1 0 0]; [0 1 0]; [0 0 1]]};
%plots{6}.plotcol={[.7*[1 1 1]; zeros(3,3)]};
%plots{7}.range= [8,5:7]; plots{7}.outname= [id '_famArm123'];
%plots{7}.plotcol=[[0 0 0]; [1 0 0]; [0 1 0]; [0 0 1]];
%plots{7}.range= [5:7]; plots{7}.outname= [id '_novelArm'];
%plots{7}.plotcol=[[1 0 0]; [0 1 0]; [0 0 1]];

%plots{7}.range= [1:3]; plots{7}.outname= [id '_novelArm'];
%plots{7}.plotcol=[[1 0 0]; [0 1 0]; [0 0 1]];

%plots{8}.range= [1:4]; plots{8}.outname= [id '_novelArm'];
%plots{8}.plotcol=[[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]];
%plots{8}.plotcol=[[1 0 0];];

subtract= 'none'; % subtract from each cell across time: 'none', 'first' ,'mean'
maxTimeDiff= -1;  % error if time difference larger than this, "-1" means no err

adaptid= 'nonripple_xp6_t005';

%timeid= 'pass';
%timelabel= 'passes';
%maxt= 21;
timeid= 'occ';
timelabel= 'occupancy (s)';
maxt= 180;

%% setting up the datasets

%sets{1}.selectid= [id '-famArm'];
%sets{1}.adaptid= adaptid;
%sets{1}.title= 'famArm';
%sets{1}.label= 'fam arm';

%sets{2}.selectid= [id '-novelArm'];
%sets{2}.adaptid= adaptid;
%sets{2}.title= 'novel';
%sets{2}.label= 'novel';

%sets{3}.selectid= [id '-fam'];
%sets{3}.adaptid= adaptid;
%sets{3}.title= 'fam';
%sets{3}.label= 'fam conf';

sets{1}.selectid= [id '-fam'];
sets{1}.adaptid= adaptid;
sets{1}.title= 'fam';
sets{1}.label= 'fam conf';

sets{2}.selectid= [id '-famArm'];
sets{2}.adaptid= adaptid;
sets{2}.novelDays= 1;
sets{2}.title= 'fam1';
sets{2}.label= 'fam 1';

sets{3}.selectid= [id '-famArm'];
sets{3}.adaptid= adaptid;
sets{3}.novelDays= 2;
sets{3}.title= 'fam2';
sets{3}.label= 'fam 2';

sets{4}.selectid= [id '-famArm'];
sets{4}.adaptid= adaptid;
sets{4}.novelDays= 3;
sets{4}.title= 'fam3';
sets{4}.label= 'fam 3';

sets{5}.selectid= [id '-novelArm'];
sets{5}.adaptid= adaptid;
sets{5}.novelDays= 1;
sets{5}.title= 'novel1';
sets{5}.label= 'novel 1';

sets{6}.selectid= [id '-novelArm'];
sets{6}.adaptid= adaptid;
sets{6}.novelDays= 2;
sets{6}.title= 'novel2';
sets{6}.label= 'novel 2';

sets{7}.selectid= [id '-novelArm'];
sets{7}.adaptid= adaptid;
sets{7}.novelDays= 3;
sets{7}.title= 'novel3';
sets{7}.label= 'novel 3';

sets{8}.selectid= [id '-famArm'];
sets{8}.adaptid= adaptid;
sets{8}.title= 'famArm';
sets{8}.label= 'fam arm';


newsets= {}; nS= 0;
for is=1:length(sets)
    if ~isempty(sets{is}); nS= nS+1; newsets{nS}= sets{is}; end
end
sets= newsets;

for iS=1:length(sets); NSEL{iS}={}; end

%for iS=1:length(sets); NSEL{iS}=nan; end
%NSEL{5}= [100]; % good examples of single cells, pf9
%plots{2}.range= [5]; plots{2}.outname= [id '_novel1_100'];
%plots{2}.plotcol=[[1 0 0]];

%NSEL{2}= [101]; NSEL{5}= [40]; % good examples of single cells, pf9
%plots{2}.range= [2,5]; plots{2}.outname= [id '_day1_single_101_40'];

%% setting up the statistics to be extracted

opt{1}.name= 'LinCorrOffset';
opt{1}.title= 'coefAbs';
opt{1}.label= 'correl coef';
opt{1}.n= 1;
opt{1}.ndeg= 3;
opt{1}.timeid= timeid;
opt{1}.timelabel= timelabel;
opt{1}.maxt= maxt;

%opt{2}.name= 'LinCorrShift';
%opt{2}.title= 'coefOpt';
%opt{2}.label= 'correl coef';
%opt{2}.n= 1;
%opt{2}.ndeg= 3;
%opt{2}.timeid= timeid;
%opt{2}.timelabel= timelabel;
%opt{2}.maxt= maxt;

%opt{2}.name= 'PeakRate';
%opt{2}.title= 'peakrate';
%opt{2}.label= 'peak rate (Hz)';
%opt{2}.timeid= timeid;
%opt{2}.timelabel= timelabel;
%opt{2}.ndeg= 3;

%opt{2}.name= 'Integral2d';
%opt{2}.title= 'Volume';
%opt{2}.label= 'integ. volume';
%opt{2}.ndeg= 3;
%opt{2}.timeid= timeid;
%opt{2}.maxt= maxt;
%opt{2}.mean= 1;
%opt{2}.norm= 1;

%opt{1}.name= 'MeanTheta';
%opt{1}.title= 'MeanThetaSp';
%opt{1}.label= 'mean theta (rad)';
%opt{1}.circ= 1;
%opt{1}.spikes= 1;
%opt{1}.ndeg= 3;
%opt{1}.timeid= timeid;
%opt{1}.timelabel= timelabel;
%opt{1}.maxt= maxt;

%opt{2}.name= 'ThetaDispersion';
%opt{2}.title= 'CircVar';
%opt{2}.label= 'circ var';
%opt{2}.spikes= 1;
%opt{2}.ndeg= 3;
%opt{2}.timeid= timeid;
%opt{2}.timelabel= timelabel;
%opt{2}.maxt= maxt;


CircNames= {'MeanChange', 'circvar', 'MeanTheta', 'change weighted by min', 'change weighted  by mean'};
CircLabels= {'MeanChange', 'circular var', 'MeanTheta', 'change weighted by min', 'change weighted  by mean'};
%for j=1:length(CircNames)
%for j=2;
for j=[];
    ind=2;
    opt{j}.name= 'CircMeasures';
    opt{j}.n= ind; % [mean change, mean dispersion, mean angle, change weighted by min, change w.  by mean]
    opt{j}.ndeg= 3;
    opt{j}.title= CircNames{ind};
    opt{j}.label= CircLabels{ind};
    if(ind==3) opt{j}.circ= 1; end
    opt{j}.timeid= timeid;
    opt{j}.timelabel= timelabel;
    opt{j}.maxt= maxt;
end


MomentNames= {'meanx', 'meany', 'varx', 'vary', 'medianx', 'mediany', 'mutual'};
MomentLabels= {'meanx', 'meany', 'varx', 'vary', 'medianx', 'mediany', 'mutual info (bits)'};
mapind(4)=1;
mapind(5)=3;
mapind(6)=7;
%for j=1:length(MomentNames)
%for j=4:6
%for j=3;
for j=[];
%    ind=mapind(j);
ind=7;
    opt{j}.name= 'Moments2d';
    opt{j}.n= ind; % [mean change, mean dispersion, mean angle, change weighted by min, change w.  by mean]
    opt{j}.ndeg= 3;
    opt{j}.title= MomentNames{ind};
    opt{j}.label= MomentLabels{ind};
    if(ind==2) opt{j}.circ= 1; end
    opt{j}.timeid= timeid;
    opt{j}.timelabel= timelabel;
    opt{j}.maxt= maxt;
end


global T V

% change titles
for iO=1:length(opt)
    switch subtract
    case 'first'
        opt{iO}.title= [opt{iO}.title '-firstsub'];
    case 'mean' % subtract cell's mean across time
        opt{iO}.title= [opt{iO}.title '-meansub'];
    end
end

if(rerun)
    setRoot;
%    olddir= pwd;
%    cd([root '/work/'])


    %%%%%%%%%% collect stats %%%%%%%%%% 
    V={}; T={};
    for iS=1:length(sets)
        load([root '/data/analist-' sets{iS}.selectid])
        nC= length(analist.rat);
        if isempty(NSEL) | isempty(NSEL{iS}); 
            nsel=[1:nC];  
        elseif ~isfinite(NSEL{iS})
            nsel= [];
        else
            nsel= NSEL{iS};
        end
        if ~isfield(sets{iS}, 'novelDays') | isempty(sets{iS}.novelDays)
            selectNovelDays= 0;
        else 
            selectNovelDays= 1;
        end
        if iS==1 | ~strcmp(sets{iS}.selectid, sets{iS-1}.selectid)
            tm= allTimes(timeid, sets{iS}.selectid, maxt);
        end
        nPF= 0;
        for iO=1:length(opt)
            V{iS}{iO}= nan*ones(nC, maxt);
            T{iS}{iO}= nan*ones(nC, maxt);
        end
        oldrat= ''; oldd= 0;

        pf= allPlaceFields(sets{iS}.selectid);

        for iC= 1:nC
            isel= find(iC== nsel);
            if isempty(isel); continue; end

            rat= analist.rat{iC};
            if ~isempty(RATSEL) & ~strcmp(rat, RATSEL); continue; end
            
            num= analist.cellnum(iC,:); d=num(1); e=num(2); tet=num(3); c=num(4);
            if ~strcmp(oldrat, rat) 
                oldrat= rat; 
%                cd([root '/' rat '/results/' sets{iS}.adaptid]);
%                setLocalOptions;
%                fmaux.selectid= sets{iS}.selectid;
%                fmaux.select=[ fmaux.data2dir '/select-' fmaux.selectid];
%                load(['stats-' sets{iS}.selectid]);
                load([root '/' rat '/results/' sets{iS}.adaptid '/stats-' sets{iS}.selectid]);
            end
            traj= analist.traj(iC);


            if selectNovelDays & ~ismember(analist.day(iC), sets{iS}.novelDays)
%            if selectNovelDays & ~ismember(day(iC), sets{iS}.novelDays)
                continue;
            end

            nPF= nPF+1;
    %        fprintf(1,'rat= %s, n= %d\n', rat, iC);

            for iO=1:length(opt)
                if isempty(opt{iO}); continue; end

                if ~isfield(opt{iO}, 'spikes') | ~opt{iO}.spikes
            %            statstime= stats{d}{e}{tet}{c}.x{analist.iana(iC)}.time;
                    statstime= stats{d}{e}{tet}{c}.([opt{iO}.name '_time']){analist.iana(iC)};
                    switch opt{iO}.name
                    case {'LinCorr', 'LinCorrShift'}
                        fac= -[1 -1 1 -1];
                        svar= fac(traj+1)*stats{d}{e}{tet}{c}.(opt{iO}.name){analist.iana(iC)}(:,1);
%                        svar=
%                        abs(stats{d}{e}{tet}{c}.(opt{iO}.name){analist.iana(iC)}(:,1));
                    case {'LinCorrOffset'}
%                        fac= -[1 -1 1 -1];
%                        svar=
%                        fac(traj+1)*stats{d}{e}{tet}{c}.(opt{iO}.name){analist.iana(iC)(:,1)};
                        svar=
                        abs(stats{d}{e}{tet}{c}.(opt{iO}.name){analist.iana(iC)}(:,1));
                    case 'Moments2d'
                        svar= stats{d}{e}{tet}{c}.Moments2d{analist.iana(iC)}(:,opt{iO}.n);
                        if(opt{iO}.n== 1) 
                            fac= [1 -1 1 -1];
                            svar= fac(traj+1)*svar; 
                            if(fac(traj+1)<0) svar= svar+240; end
                        end
                    case 'CircMeasures'
                        svar=
                        stats{d}{e}{tet}{c}.CircMeasures{analist.iana(iC)}(:,opt{iO}.n);
                    otherwise
                        svar= stats{d}{e}{tet}{c}.(opt{iO}.name){analist.iana(iC)};
                    end
                else % spike analyses
                    if ~strcmp(rat, oldrat) | d~=oldd
                        oldrat= rat; oldd= d;
                        load(sprintf('%s/%s/data2/spikedata%.2d.mat', root, rat, d));
                    end
                    statstime= tm{iC};
                    svar= auxSpikeAna(spikedata{d}{e}{tet}{c}, pf(iC,:), opt{iO},statstime);
                end

                t= tm{iC};
                nt= min(length(t), maxt);
                T{iS}{iO}(iC,1:nt)= t(1:nt);
%                    if any(isnan(t)); keyboard; end
                for it=1:nt
                    istats= min(find(t(it) <= statstime+1e-5));
                    if isempty(istats); continue; end
                    if(maxTimeDiff > 0 & abs(statstime(istats)-t(it)) > maxTimeDiff) 
                        error('could not find good time'); 
                    end
                    V{iS}{iO}(iC, it)= svar(istats);
                end
                switch subtract
                case 'first'
                    stmp= V{iS}{iO}(iC, :);
                    V{iS}{iO}(iC, :)= stmp-stmp(min(find(isfinite(svar))));
%                    V{iS}{iO}(iC, :)= abs(stmp-stmp(min(find(isfinite(svar)))));
                case 'mean' % subtract cell's mean across time
                    stmp= V{iS}{iO}(iC, :);
                    if isfield(opt{iO},'circ') & opt{iO}.circ
                        % circular data
                        V{iS}{iO}(iC, :)= stmp-circstat(stmp);
                    else
                        V{iS}{iO}(iC, :)= stmp-nanmean(stmp);
                    end

                end

                if(iS==2 & iC==110); keyboard; end
            end % for iC

        end % for iO, statistics 
        fprintf(1, '%3d placefields selected in %s\n', nPF, sets{iS}.title);
    end % for iS, sets


    for iO=1:length(opt)
        if isfield(opt{iO}, 'norm') & opt{iO}.norm
            opt{iO}.title= [opt{iO}.title '-norm'];
        end
        if isfield(opt{iO}, 'mean') & opt{iO}.mean
            V{iS}{iO}= nanmean(V{iS}{iO}')';
        end
    end

%    cd(olddir);
end

%%%%%%%%%% show single var stats %%%%%%%%%% 

for iS= 1:length(sets)
    for iO=1:length(opt)
        vall= V{iS}{iO};
        tall= T{iS}{iO};

        opt{iO}.figname= [opt{iO}.name '-' opt{iO}.timeid '-' sets{iS}.adaptid '-' sets{iS}.selectid];

        switch opt{iO}.timeid
        case 'passes'
            opt{iO}.xstr= 'pass';
        case 'occ'
            opt{iO}.xstr= 'occ. (sec)';
        case 'acc'
            opt{iO}.xstr= 'acc. (sec)';
        case 'minocc'
            opt{iO}.xstr= 'occ. (min)';
        case 'minacc'
            opt{iO}.xstr= 'acc. occ. (min)';
        case 'first'
            opt{iO}.xstr= 'first';
        end

    %    if(opt{iO}.plotsingle) 
        if(plotsingle)  auxPlotSingle(vall, opt{iO}, sets{iS}); end
        if(plotmean)    auxPlotMean(vall, opt{iO}); end
        if(plotPDF)     auxPlotPDF(vall, opt{iO}); end

        %%%%%%%%%% calc and show derivatives %%%%%%%%%% 

    %    if isfield(opt{iO}, 'derivative') & opt{iO}.derivative
        if derivative
            vall= auxCalcDerivative(vall, opt{iO});
            opt{iO}.figname= ['delta-' opt{iO}.figname];
            opt{iO}.name= (['\Delta ' opt{iO}.name]);
        %    opt{iO}.maxt=5;
            auxPlotCDF(vall,opt{iO},12);
            auxPlotIND(vall,opt{iO},12);
        end
    end % for iO, statistics

    %%%%%%%%%% show two var stats %%%%%%%%%% 

    if(plotcorr | plotpair)  
        [corr,p]= auxPlotCorr(V{iS}, opt, sets{iS}.title); 
    end

    if(plotpair)  auxPlotPair(V{iS}, opt, corr, p); end

    if(plotscatter)  auxPlotScatter(V{iS}, opt, iS); end

end % for iS, sets

for iplot=1:length(plots)
    if isempty(plots{iplot}); continue; end
    range= plots{iplot}.range;
    if isempty(range); range= [1:length(sets)]; end

    for iO=1:length(opt)
        Otmp= opt{iO};
        Otmp.outname= plots{iplot}.outname;
%        Otmp.outname= [Otmp.title '_' plots{iplot}.outname];
%        Otmp.title= [Otmp.title '_' plots{iplot}.outname];
        Otmp.plotcol= plots{iplot}.plotcol;
        Vtmp= {}; Stmp= {}; itmp= 0; 
        for iS=range
            itmp=itmp+1;
            Vtmp{itmp}= V{iS}{iO};
            Stmp{itmp}= sets{iS};
        end
        if(plotVsSets)  auxPlotVsSets(Vtmp, Otmp, Stmp, timeind); end
        if(plotMultiMean)  auxPlotMean(Vtmp, Otmp, Stmp); end
        if(plotSessionMean)  auxPlotSessionMean(Vtmp, Otmp, Stmp); end
        if(saveData)  auxSaveData(Vtmp, Otmp, Stmp); end
        if(plotMultiHist)  auxPlotHist(Vtmp, Otmp, Stmp, timeind); end
        if(plotMultiDist)  auxPlotMultiDist(Vtmp, Otmp, Stmp); end
        if(plotDist)  
            for it= timeind
                iv= 0; sets_t= {};
                for iv=1:length(Vtmp)
                    Vt{it}{iv}= Vtmp{iv}(:,it);
                    Vt{it}{iv}= Vt{it}{iv}(isfinite(Vt{it}{iv}));
                end
                Otmp.outname= ['t' num2str(it)];
                Otmp.plot= 'hist'; Otmp.norm= 1; Otmp.size= [1 2/3];
                Otmp.legend= 0; Otmp.ref= 1;
                auxCmpDist2(Vt{it}, Stmp, Otmp); 
            end
%                Bt{1}=Vt{1};
%                Bt{1}=Vt{1}{1};
%                Bt{2}=Vt{180}{1};
%                Stmp{2}= Stmp{1}; 
%                auxCmpDist2(Bt, Stmp, Otmp); 
        end
    end
end


function auxPlotSingle(vall, opt, set)
fh=figure;
ind= find(sum(isfinite(vall),2))';
n= length(ind);
nx= ceil(sqrt(n)); ny= ceil(n/nx);
for i=1:n
    iC=ind(i);
    subplot(nx, ny, i);
    ph= plot(vall(iC,:));
    title(num2str(iC))
%    axis([0 180 0 1])
%    xlabel(opt.timelabel);
%    ylabel(opt.label);
%    set(gca, 'XLim', [1 180])
    axis tight
%    axis off
end
%orient landscape
%print(gcf, '-depsc2', [opt.figname '-single']);
%set(fh, 'PaperPosition', [0 0 4 2]);
figname= [opt.title '-' set.title '-' opt.timeid];
set(fh, 'Name', figname);
%myprint('mini', figname);

function auxPlotMultiDist(V, opt, sets)

opt.xstr= opt.timelabel;
opt.ystr= opt.label;
for is=1:length(sets)
    figure
    [i,x]= find(isfinite(V{is}));
    y= V{is}(find(isfinite(V{is})));
    [pval, r, n]= linrel(x,y,opt,0.05);
end
figname= ['multiDist-' opt.title '-' opt.timeid];
set(gcf, 'Name', figname);

function auxPlotMean(V, opt, sets)

nmin= 1;

multi= 1;
if ~iscell(V)
    V= {V};
    multi= 0;
end
nV= length(V);
maxt= opt.maxt;

circ= 0;
if(isfield(opt, 'circ') & opt.circ) circ= 1; end

fillcol= rgb2hsv(opt.plotcol);
fillcol(:,2)= .4;
fillcol= hsv2rgb(fillcol);
ind= find(all(fillcol==0,2));
fillcol(ind,:)= .7*ones(length(ind),3);
clear ind

figh= figure;
hold on

if ~circ
    for iV=1:nV
        nval= sum(isfinite(V{iV}));

        nval(find(nval==0))= nan;

        mvar{iV}= nanmean(V{iV});
        stdvar= nanstd(V{iV});
        sevar= stdvar./sqrt(nval);

        ind{iV}= find(isfinite(mvar{iV}) & isfinite(stdvar) & nval >= nmin);
        %fh= fill([1:maxt maxt:-1:1], [mvar{iV}(1:maxt)-stdvar(1:maxt),
        %fliplr(mvar{iV}(1:maxt)+stdvar(1:maxt))], fillcol);
        %fh= fill([ind{iV} fliplr(ind{iV})], [mvar{iV}(ind{iV})-stdvar(ind{iV}),
        %fliplr(mvar{iV}(ind{iV})+stdvar(ind{iV}))], fillcol);
%        if(nV<=2 | iV<=1)
        if(nV<=2)
            fh= fill([ind{iV} fliplr(ind{iV})], [mvar{iV}(ind{iV})-sevar(ind{iV}), fliplr(mvar{iV}(ind{iV})+sevar(ind{iV}))], fillcol(iV,:));
            set(fh, 'EdgeColor', fillcol(iV,:));
        end
    %    plot(ind{iV}, mvar{iV}(ind{iV})-sevar(ind{iV}), 'Color', opt.plotcol(iV,:), '-']);
    %    plot(ind{iV}, mvar{iV}(ind{iV})+sevar(ind{iV}), 'Color', opt.plotcol(iV,:), '-']);
        %plot(1:maxt, mvar{iV}(1:maxt), 'k', 'LineWidth', 3);


        %% regression
%        fi=find(isfinite(V{iV}));
%        Xtmp= ones(size(V{iV},1),1)*[1:maxt]; X= Xtmp(fi); Y= V{iV}(fi);
%        [b,bint,r,rint,stats] = regress(Y,[X ones(length(fi),1)]);
%        fprintf(1,'%s:\tm= %.2g, b= %.2g, r^2= %.2f, p= %.2g\n', sets{iV}.title, b(1), b(2), stats(1), stats(3));

    end

    for iV=1:nV
%        if iV<3; continue; end
        hp(iV)= plot(ind{iV}, mvar{iV}(ind{iV}));
        set(hp(iV), 'Color', opt.plotcol(iV,:));
    end

    set(gca, 'Xlim', [0 maxt]);
%    set(gca, 'ylim', [0.25 0.47]);
%axis tight
%    line([90, 90], get(gca, 'YLim'), 'Color', 'k', 'LineStyle', ':');
 

else % circ
    dp= {};
    for iV=1:nV
        nval= sum(isfinite(V{iV}));

        [mvar{iV}, dp{iV}]= circstat(V{iV});
        ind{iV}= find(isfinite(mvar{iV}) & nval >= nmin);

        indtmp= ind{iV};
        mtmp= mvar{iV}(ind{iV}); 
        % subtract population's mean acros time
        moff= circstat(mvar{iV}(ind{iV}));
        mtmp= mod(mtmp-moff+pi, 2*pi)-pi;
        mtmp= mod(mtmp+pi, 2*pi)-pi;
%        indtmp= [indtmp; indtmp];
%        mtmp= mod(mtmp -moff, 2*pi)-2*pi;
%        mtmp= [mtmp; 2*pi+mtmp];
        hptmp= plot(indtmp', mtmp', 'o', ...
            'MarkerSize', 4, 'MarkerFaceColor', opt.plotcol(iV,:),...
            'Color', opt.plotcol(iV,:));
        hp(iV)= hptmp(1);
%        set(gca, 'ylim', 2*pi*[-1,1]);
    end
end

xlabel(opt.timelabel);
ylabel(opt.label);
if multi; 
%    if nV<3
%        lstr= {};
%        for(iV=1:nV) lstr{iV}= sets{iV}.title; end
%        legend(hp, lstr, 0)
%    end

    figname= [opt.title '_' opt.outname '-' opt.timeid]; %%@@
%    figname= [opt.title];

    global RATSEL
    if ~isempty(RATSEL); figname= [figname '_' RATSEL]; end

    set(figh, 'Name', figname);
    if(nV<=3 | iV<=1)
        myprint( 2*[1 2/3], figname);
    else
        myprint('small', figname);
    end
%    print(gcf, '-depsc2', figname);
else
    set(figh, 'Name', opt.figname);
    myprint('landscape', [opt.figname '-mean']);
%    print(gcf, '-depsc2', figname);
end

if circ
    figh=figure;
    hold on
    for iV=1:nV
        plot(ind{iV}, 1-dp{iV}(ind{iV}), 'Color', opt.plotcol(iV,:), 'LineWidth', 4);
    end
    legend(lstr, 0)
    figname= ['multiDisp-' opt.title '-' opt.timeid];
    set(figh, 'Name', figname);
    print(gcf, '-depsc2', figname);
end

%for iV=1:nV
%    figure
%    [i,j]= find(isfinite(V{iV}));
%    i= find(isfinite(V{iV}));
%    x= V{iV}(i);
%    linrel(j, x);
%end

function auxPlotPDF(vall, opt)
% 
[h, xbin]= hist(vall);
for it=1:opt.maxt
    h(:,it)= h(:,it)/sum(isfinite(vall(:,it)));
end
figure
if(opt.maxt >=3) 
    ind= [1 ceil(opt.maxt/2) opt.maxt];
elseif(opt.maxt ==2) 
    ind= [1 2];
else
    ind= 1;
end
plot(xbin, h(:, ind), 'LineWidth', 2);
legend(num2str(ind'));
xh= xlabel(opt.name);
ylabel('fraction');
print(gcf, '-depsc2', [opt.figname '-pdf']);

if(opt.maxt < 20)
    figure
    plot(xbin, h(:,1:opt.maxt), 'LineWidth', 2);
    legend(num2str([1:opt.maxt]'))
    xlabel(opt.name);
    ylabel('fraction');
end
    
fprintf(1, 'Testing distributions at beginning and end:\n');
%cmpSamples(vall(:,1), vall(:,opt.maxt));
cmpSamples(vall(:,1), vall(:,150));

function auxPlotCDF(vall, opt, lump)
% 
if nargin<3; lump=1; end
figure
nC= size(vall,1);
%mt= min(7,opt.maxt);
mt= floor(opt.maxt/lump);
if(mt*lump~= opt.maxt) error('mt*lump~= opt.maxt'); end
h= sort(reshape(vall, nC*lump, mt));
n= nan*zeros(nC*lump, mt);
for im=1:mt
    nind= length(find(isfinite(h(:,im))));
    n(1:nind, im)= [1:nind]'/nind;
end
plot(h,n, 'LineWidth', 2);
%legend(num2str([1:mt]'))
xlabel(opt.name);
ylabel('cum. fraction');
print(gcf, '-depsc2', [opt.figname '-cdf']);

%cmpSamples(vall(:,1), vall(:,mt));

function h= auxCalcDerivative(vall, opt)

nC= size(vall,1);
pad= floor(opt.ndeg/2)+1;
h= nan*ones(size(vall));

for iC=1:nC
    % polynomial fit to f(t)
    valid= find(isfinite(vall(iC,:)));
    if length(valid)<= opt.ndeg
        vall(iC,valid)= nan;
        continue;
    end
    for j=pad+1:length(valid)-pad
        ind= valid(j-pad:j+pad);
        v= vall(iC,ind);
        p= polyfit(ind, v, opt.ndeg);

        % get derivative of polynomial
        dp= [opt.ndeg:-1:1] .* p(1:end-1);

        % calc derivative of polynomial
        h(iC,valid(j))= polyval(dp,valid(j));
    end

    %figure; hold on
    %tmp=ind(1):.1:ind(end);
    %plot(tmp, polyval(p,tmp), 'k');
    %plot(tmp, polyval(dp,tmp),'r');

end

function auxPlotIND(vall,opt,lump);
figure
if nargin<3; lump=1; end
mt= floor(opt.maxt/lump);
nC= size(vall,1);
h= sort(reshape(vall, nC*lump, mt));
x= nanmean(abs(h));
%x= x/ max(x);
bar(x)

xlabel('"t"');
ylabel(['< ' opt.name ' >']);
print(gcf, '-depsc2', [opt.figname '-IND']);

function vall= auxNorm(vall);
maxt= size(vall,2);
n= nanmean(vall')'*ones(1,maxt);
vall= vall./n;

function [corr, p]= auxPlotCorr(V, opt, selid)
if length(V)~= 2; error('This function can only handle two variables'); end
nC= size(V{1},1);
maxt= size(V{1},2);
corr= nan*ones(nC,1);
p= nan*ones(nC,1);
for iC=1:nC
    ind= find(isfinite(V{1}(iC,:)) & isfinite(V{2}(iC,:)));
    if length(ind) < 3;
        corr(iC)= nan;
        p(iC)= nan;
    else
        [tmp tmpp]= corrcoef(V{1}(iC,ind), V{2}(iC,ind));
        corr(iC)= tmp(1,2);
        p(iC)= tmpp(1,2);
    end
end
%figure
%subplot(3,1,1)
%hist(corr,20)
%xlabel('corr');
%subplot(3,1,2)
%hist(p,20)
%xlabel('p-value');

nvalid= sum(isfinite(corr));
%sigind= find(p < 0.05);
sigind= find(p < 0.01);
nsig= length(sigind);

%subplot(3,1,3)
%hist(corr(sigind),20)
%xlabel('sig corr');

figure
x= linspace(-1,1+1e-10,21);
h= histc(corr(sigind),x);
h(:,2)= histc(corr,x);
h(:,2)= diff(h,1,2);
hb= bar(x+(x(2)-x(1))/2,h, 'stacked');
set(hb(2),'FaceColor', 'w')
set(gca, 'XLim', [-1 1])
xlabel('correlation');
ylabel('count');

nind= sigind(find(corr(sigind)<0));
pind= sigind(find(corr(sigind)>0));
nneg=  sum(corr(sigind)<0);
percneg=  nneg/nsig*100;


%mean(corr(sigind))
%mean(corr(nind))
%mean(corr(pind))

fprintf(1, '%.1f %% of corrcoef sig\n', nsig/nvalid*100);
fprintf(1, '   of those %.1f %% pos, %.1f %% neg\n', 100-percneg, percneg);


%set(gcf, 'PaperPosition', [0 0 6 4.5]);
print(gcf, '-depsc2', ['corr-' opt{1}.title '-' opt{2}.title '-' selid '-' opt{1}.timeid]);

for iV=1:2
%for iV=1:-1
    figure; 
    subplot(2,2,1)
    plot(nanmean(V{iV}(sigind,:)))
    title('sig')
    subplot(2,2,2)
    plot(nanmean(V{iV}(nind,:)))
    title('sig and neg')
    subplot(2,2,3)
    plot(nanmean(V{iV}(pind,:)))
    title('sig and pos')
    subplot(2,2,4)
    plot(nanmean(V{iV}))
    title('all')
end

function auxPlotPair(V, opt, corr, p)
figure
nC= size(V{1},1);
nx= ceil(sqrt(nC)); ny= ceil(nC/nx);
for iC=1:nC
    if(p(iC) < 0.01)
        plotsym= '-';
    else
        plotsym= ':';
%        continue
    end
    ind= find(isfinite(V{1}(iC,:)) & isfinite(V{2}(iC,:)));
    if isempty(ind); continue; end;

    subplot(nx, ny, iC);
    plot([V{1}(iC,ind)' V{2}(iC,ind)'], plotsym);
    title([num2str(iC) ', ' num2str(corr(iC))])
%    axis tight
    set(gca, 'ylim', [0 2]);
    axis off
end

%orient landscape
print(gcf, '-depsc2', ['pair-' opt{1}.figname '-' opt{2}.name]);

function auxPlotPairOnly(V, opt,nsel)
figure
nS= length(nsel)
nx= ceil(sqrt(nS)); ny= ceil(nS/nx);
plotsym= '-';
for iS=1:nS
    iC= nsel(iS);
    subplot(nx, ny, iS);
    ind= find(isfinite(V{1}(iC,:)) & isfinite(V{2}(iC,:)));
    plot([V{1}(iC,ind)' V{2}(iC,ind)'], plotsym);
    axis tight
    axis off
end

%orient landscape
print(gcf, '-depsc2', ['pair-' opt{1}.figname '-' opt{2}.name]);

function auxPlotScatter(V, opt, is)
%if(is==1) 
    figure; 
%end
global plotcol
%v1= nanmean(V{1}')'; v2= nanmean(V{2}')';
v1= V{1}(:,1); v2= V{2}(:,1);
popt.xstr= opt{1}.name;
popt.ystr= opt{2}.name;
popt.plotcol= plotcol(is,:);
linrel(v1,v2,popt);

%nC= size(V{1},1);
%maxt= size(V{1},2);
%ind= find(isfinite(V{1}) & isfinite(V{2}));
%nind= length(ind);
%v1= V{1}(ind); v2= V{2}(ind);

%plot(v1, v2, '.');
%xlabel(opt{1}.name);
%ylabel(opt{2}.name);
%[tmp, tmpp]= corrcoef(v1,v2);
%C= tmp(1,2);
%p= tmpp(1,2);
%title(sprintf('corr= %.2f, p= %.4g', C, p));
%[b,bint,r,rint,stats] = regress(v2,[v1 ones(nind,1)]);
%title(sprintf('r^2= %.2f, p= %.4g', stats(1), stats(3)));

myprint('small', ['scatter-' opt{1}.figname '-' opt{2}.name]);

function svar= auxSpikeAna(sp, pf, opt, t);
nt= length(t);
svar= nan*zeros(1,nt);

indall= find(pf(1)<sp.linpos & sp.linpos<pf(2));
if(0)
    ind= {};
    t= [0, t];
    for it=1:nt
        ind{it}= find(t(it)<sp.time(indall) & sp.time(indall)<t(it+1));
    end

    switch opt.name
    case {'MeanTheta'}
        for it=1:nt;
            th= sp.phase(ind{it});
            svar(it)= mod(circstat(th),2*pi);
        end
    case {'ThetaDispersion'}
        for it=1:nt;
            th= sp.phase(ind{it});
            [m,svar(it)]= circstat(th);
        end
    otherwise
        error(['unknown variable ''' opt.name '''']);
    end
else

    switch opt.name
    case {'MeanTheta'}
        svar(1)= mod(circstat(sp.phase(indall)),2*pi);
    case {'ThetaDispersion'}
        [m,svar(1)]= circstat(sp.phase(indall));
    otherwise
        error(['unknown variable ''' opt.name '''']);
    end
end

function auxPlotHist(v, opt, sets, t)

alpha=0.05;
nt= length(t);
nC= size(v{1},1);
maxt= size(v{1},2);
if sum(t>maxt) | sum(t<1); error('requested invalid time'); end
nx= ceil(sqrt(nt)); ny= ceil(nt/nx);
nv= length(v);
testval= 0;
if(strcmp(opt.title, 'MeanDispersion')) testval= 1; end

if(0)
for iv=1:nv
    fh= figure;
    for it=1:nt
        subplot(nx,ny,it);
        [n,x]= hist(v{iv}(:,t(it)),20);
        bar(x,n);
    %    plot(x,n, 'k-', 'LineWidth', 2)
        xlabel(opt.label);
        ylabel('count');
        axis tight

%        [h, p]= ttest(v{iv}(:,t(it)), testval, alpha, 0);
%        title(sprintf('ti= %d, p= %.4g', t(it), p));
    end

    figname= ['hist-' opt.title '-' sets{iv}.title '-' opt.timeid];
    set(fh, 'Name', figname);
%    set(gcf, 'PaperPosition', [0 0 4 3]);
    print(gcf, '-depsc2', figname);
end
else

lstr= {}; lstrcat= sets{1}.title;
for(iv=1:nv) 
    lstr{iv}= sets{iv}.title; 
    if(iv>1) lstrcat= [lstrcat '-' lstr{iv}]; end
end

for it=1:nt
    figure
    vtmp= zeros(nC,nv);
    for iv=1:nv
        vtmp(:,iv)= v{iv}(:,t(it));
    end
    [n,x]= hist(vtmp, 20);

    subplot(nx,ny,it);
    bar(x,n);
    xlabel(opt.label);
    ylabel('count');
    legend(lstr);

    figname= sprintf('hist-%s-%s-%s', opt.title, lstrcat, opt.timeid);
    set(gcf, 'Name', figname);
%    set(gcf, 'PaperPosition', [0 0 4 3]);
%    if nv==2
%        m1= v{1}(:,t(it)); m1=m1(isfinite(m1));
%        m2= v{2}(:,t(it)); m2=m2(isfinite(m2));
%        p= ranksum(m1, m2, alpha);
%        title(sprintf('ranksum test, p=%.2g\n', p)); 
%    end
end
myprint('large', figname);
end

fprintf(1, 'variable  %s\n', opt.title);

ref= 1;
for it=1:nt
    fprintf(1, 't= %d\n', t(it));
    m1= v{ref}(:,t(it)); m1=m1(isfinite(m1));
    for iv=1:nv
        if(iv==ref); continue; end
        m2= v{iv}(:,t(it)); m2=m2(isfinite(m2));
%        [h, p]= ttest2(m1,m2, alpha, 0);
%        fprintf(1, 'two sample t-test F(%s)=F(%s), p=%.4g\n', sets{ref}.title, sets{iv}.title, p); 
            p= ranksum(m1, m2, alpha);
            fprintf(1, 'ranksum test F(%s)=F(%s), p=%.4g\n', sets{ref}.title, sets{iv}.title, p); 
    end

end


function auxPlotVsSets(v, opt, sets, t)

global plotcol

nt= length(t);
nv= length(v);

m= zeros(nv,nt);
sd= zeros(nv,nt);
se= zeros(nv,nt);
n= zeros(nv,nt);

x= zeros(nv,1);
for iv=1:nv
    for it=1:nt
        vars= v{iv}(:,t(it));
        vars= vars(isfinite(vars));
        n(iv,it)= length(vars);

        m(iv,it)= mean(vars);
        sd(iv,it)= std(vars);
        se(iv,it)= sd(iv,it)/ sqrt(n(iv,it));

    end
    if isfield(sets{iv}, 'num');
        x(iv)= sets{iv}.num;
    end
end

figh= figure;
if isfield(sets{iv}, 'num');
    errorbar(x*ones(1,nt), m, se);
    pad= min(diff(x));
    set(gca, 'XLim', [x(1)-pad x(end)+pad]);
    xlabel([sets{1}.numlabel ' (ms)'])
    legend(num2str(t'));
else
    h= bar(m');
    hold on
    lstr= {};
    for iv=1:nv; 
        set(h(iv), 'FaceColor', plotcol(iv,:)); 
        lstr{iv}= sets{iv}.title;
        xd= get(get(h(iv), 'Children'), 'XData');
        for it=1:size(xd,2)
%            plot(mean(xd(2:3,it)), .3, 'o');
            errorbar(mean(xd(2:3,it)), m(iv,it), se(iv,it), 'Color', 0.4*[1 1 1]);
        end
    end
    hold off
%    legend(lstr, 'Location', 'BestOutside', 'Orientation', 'horizontal');
    set(gca, 'XTickLabel', t);
    xlabel(opt.timelabel);
end
ylabel(opt.label);
figname= ['sum-' opt.title ];
set(figh, 'Name', figname);
myprint('small', figname);


function auxSaveData(V, opt, sets)

for i=1:length(V)
    coefAbs.ind= find(isfinite(V{i}(:,1)));
    coefAbs.val= nanmean(V{i}(coefAbs.ind,1:60)')';
    coefAbs.sets= sets{i};
    save(['coefAbs_' sets{i}.title], 'coefAbs');
end

function auxPlotSessionMean(V, opt, sets)
global RATSEL
if ~isempty(RATSEL); opt.outname= [opt.outname '_' RATSEL]; end

    %%@@
for i=1:length(V)
    tmpV= nanmean(V{i}');
    Zd{i}= tmpV(isfinite(tmpV));
end
opt.xticklabels= {'train.' 'day1'  'day2'  'day3'}; 
opt.plotcol= opt.plotcol{1};
opt.plot= 'bar';
auxCmpDist2(Zd, sets, opt);


if length(V)==6
    for d=1:3
        for ia=1:2
            j= 3*(ia-1)+d;
            tmpV= nanmean(V{j}');
            Zd{d}{ia}= tmpV(isfinite(tmpV));
        end
    end

    for d=1:3
        for ia=1:2
            Zm(d,ia)= nanmean(Zd{d}{ia});
            Zdev(d,ia)= nanstd(Zd{d}{ia})/sqrt(sum(isfinite(Zd{d}{ia})));
            if ia>1
    %             [h,p(d,ia-1)]= ttest2(Zd{d}{ia-1}, Zd{d}{ia});
                 [p(d,ia-1), h]= ranksum(Zd{d}{ia-1}, Zd{d}{ia});
    %             [h,p(d,ia-1)]= kstest2(Zd{d}{ia-1}, Zd{d}{ia});
             end
        end
    %    p(d,:)= p(d,:)* length(p(d,:));
    %    opt.outname
    end
    opt.xticklabels= {'day1'  'day2'  'day3'}; 
    auxBar(Zm, Zdev, opt, p);
    p
else
    for i=1:length(V)
        tmpV= nanmean(V{i}');
        Zd{i}= tmpV(isfinite(tmpV));
    end

    for i=1:length(V)
        Zm(i,1)= nanmean(Zd{i});
        Zdev(i,1)= nanstd(Zd{i})/sqrt(sum(isfinite(Zd{i})));
%        if ia>1
    %             [h,p(d,ia-1)]= ttest2(Zd{d}{ia-1}, Zd{d}{ia});
%                 [p(d,ia-1), h]= ranksum(Zd{d}{ia-1}, Zd{d}{ia});
    %             [h,p(d,ia-1)]= kstest2(Zd{d}{ia-1}, Zd{d}{ia});
%             end
    %    p(d,:)= p(d,:)* length(p(d,:));
    %    opt.outname
    end
    opt.xticklabels= {'train.' 'day1'  'day2'  'day3'}; 
    figure
    auxBar(Zm, Zdev, opt);
%    figure
%    opt.plotcol= opt.plotcol{1};
%    opt.ref= 3;
%    auxCmpDist2(Zd, sets, opt);
end
