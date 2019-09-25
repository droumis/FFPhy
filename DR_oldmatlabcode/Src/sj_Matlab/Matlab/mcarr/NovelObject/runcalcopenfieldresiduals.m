%% CALCULATE FIRING RATE ACROSS ALL TIME

% Animal Selection
animals = {'CorrianderNO','Cyclops','Dunphy','Fafnir','Godzilla','Grendel','Cml','Nico'};

% Epoch selection
epochfilter{1} = 'isequal($type,''run'') & isequal($session,''familiar'')';
epochfilter{2} = 'isequal($type,''run'') & isequal($session,''novel'')';
epochfilter{3} = 'isequal($type,''run'') & isequal($session,''supernovel'')';

%Define time filter
timefilter = {{'get2dstate', '($velocity > 1)'},{'getriptimes', '($nripples == 0)', [], 'cellfilter', 'isequal($area, ''CA1'')'}};
iterator = 'multicellanal';

%RUN FOR CA1
% Cell Filter
cellfilter = 'isequal($area, ''CA1'') & $meanrate < 7 & $numspikes > 50';

%RUN FOR CA3
% Cell Filter
cellfilter = 'isequal($area, ''CA3'') & $meanrate < 7 & $numspikes > 50';

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);
f = setfilterfunction(f, 'getopenfieldrates', {'spikes','pos','task'});
f = runfilter(f);

g = calcopenfieldresiduals(f);

%% Look at firing rate changes by session using residuals
all = []; fam = []; nov = []; sup = []; famobj = []; novobj = [];
for an = 1:length(g)
    for d = 1:length(g(an).resid)
        tmp = g(an).resid{d};
        all = [all; nanmean(tmp)];
        tmp = g(an).resid{d}(g(an).session{d}==1);
        fam = [fam; nanmean(tmp(tmp~=0))];
        tmp = g(an).resid{d}(g(an).session{d}==2);
        nov = [nov; nanmean(tmp(tmp~=0))];
        tmp = g(an).resid{d}(g(an).session{d}==3);
        sup = [sup; nanmean(tmp(tmp~=0))];
        
        tmp = g(an).resid{d}(g(an).session{d}==2 & g(an).quad{d}==2);
        tmp2 = g(an).resid{d}(g(an).session{d}==3 & g(an).quad{d}==2);
        famobj = [famobj; nanmean(tmp(tmp~=0)) nanmean(tmp2(tmp2~=0))];

        tmp = g(an).resid{d}(g(an).session{d}==2 & g(an).quad{d}==1);
        tmp2 = g(an).resid{d}(g(an).session{d}==3 & g(an).quad{d}==1);
        novobj = [novobj; nanmean(tmp(tmp~=0)) nanmean(tmp2(tmp2~=0))];
    end
end
all(isnan(all)) = []; fam(isnan(fam)) = []; nov(isnan(nov)) = []; sup(isnan(sup))=[];
famobj(isnan(famobj)) = []; novobj(isnan(novobj)) = [];

figure
hold on
bar(1:4,[mean(all) mean(fam) mean(nov) mean(sup)],'b')
errorbar2(1:4,[mean(all) mean(fam) mean(nov) mean(sup)],[stderr(all) stderr(fam) stderr(nov) stderr(sup)],'k')
set(gca,'xtick',1:4,'xticklabel',[{'All'},{'Fam'},{'Nov'},{'Sup'}])

%Run statistics
group = [ones(size(fam)); 2*ones(size(nov)); 3*ones(size(sup))];
[a b stats]=anovan([fam; nov; sup],group,'display','off');
multcompare(stats);

figure
hold on
bar([1:4 6:7],[mean(all) mean(fam) mean(nov) mean(sup) mean(famobj) mean(novobj)],'b')
errorbar2([1:4 6:7],[mean(all) mean(fam) mean(nov) mean(sup) mean(famobj) mean(novobj)],[stderr(all) stderr(fam) stderr(nov) stderr(sup) stderr(famobj) stderr(novobj)],'k')
set(gca,'xtick',[1:4 6:7],'xticklabel',[{'All'},{'Fam'},{'Nov'},{'Sup'},{'FamObj'},{'NovObj'}])


p = ranksum(famobj,novobj);
[h p] = ttest2(famobj,novobj);

%% Look at firing rate (residuals) over time

fam = []; ftime = []; nov = []; ntime = []; sup = []; stime = [];
for an = 1:length(g)
    for d = 1:length(g(an).resid)
        fam = [fam; g(an).resid{d}(g(an).session{d}==1)];
        tmp = g(an).time{d}(g(an).session{d}==1);
        if ~isempty(tmp)
            ftime = [ftime tmp-tmp(1)];
        end
        nov = [nov; g(an).resid{d}(g(an).session{d}==2)];
        tmp = g(an).time{d}(g(an).session{d}==2);
        if ~isempty(tmp)
            ntime = [ntime tmp-tmp(1)];
        end    
        sup = [sup; g(an).resid{d}(g(an).session{d}==3)];
        tmp = g(an).time{d}(g(an).session{d}==3);
        if ~isempty(tmp)
            stime = [stime tmp-tmp(1)];
        end
    end
end

bin = 60:60:600;
af = accumarray(lookup(ftime,bin,1),fam,[length(bin) 1],@(x) nanmean(x),NaN);
binn = 60:60:300;
an = accumarray(lookup(ntime,binn,1),nov,[length(binn) 1],@(x) nanmean(x),NaN);
as = accumarray(lookup(stime,binn,1),sup,[length(binn) 1],@(x) nanmean(x),NaN);

af_err = accumarray(lookup(ftime,bin,1),fam,[length(bin) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1),NaN);
an_err = accumarray(lookup(ntime,binn,1),nov,[length(binn) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1),NaN);
as_err = accumarray(lookup(stime,binn,1),sup,[length(binn) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1),NaN);

af = accumarray(lookup(ftime,bin,-1),fam,[length(bin) 1],@(x) nanmean(x(x~=0)),NaN);
an = accumarray(lookup(ntime,binn,-1),nov,[length(binn) 1],@(x) nanmean(x(x~=0)),NaN);
as = accumarray(lookup(stime,binn,-1),sup,[length(binn) 1],@(x) nanmean(x(x~=0)),NaN);

%% Plot example place field from CA1 and CA3

%CA1
figure
imagesc(g(2).binx{1},g(2).biny{1},g(2).rate{1}')
axis xy
colorbar('ylim',[0 12],'ytick',12)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_CA1globalrate_example.png', m, d, y);
print('-dpng', savestring)


load '/data6/monster/Cyc/Cycpos01.mat'
load '/data6/monster/Cyc/Cycspikes01.mat'

figure
subplot(1,3,1)
plot(pos{1}{2}.data(:,2),pos{1}{2}.data(:,3),'k',spikes{1}{2}{1}{1}.data(:,2),spikes{1}{2}{1}{1}.data(:,3),'r.')
set(gca,'xlim',[min(g(2).binx{1}) max(g(2).binx{1})],'Ylim',[min(g(2).biny{1}) max(g(2).biny{1})])
subplot(1,3,2)
plot(pos{1}{4}.data(:,2),pos{1}{4}.data(:,3),'k',spikes{1}{4}{1}{1}.data(:,2),spikes{1}{4}{1}{1}.data(:,3),'r.')
set(gca,'xlim',[min(g(2).binx{1}) max(g(2).binx{1})],'Ylim',[min(g(2).biny{1}) max(g(2).biny{1})])
subplot(1,3,3)
plot(pos{1}{6}.data(:,2),pos{1}{6}.data(:,3),'k',spikes{1}{6}{1}{1}.data(:,2),spikes{1}{6}{1}{1}.data(:,3),'r.')
set(gca,'xlim',[min(g(2).binx{1}) max(g(2).binx{1})],'Ylim',[min(g(2).biny{1}) max(g(2).biny{1})])

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_CA1placefield_example.pdf', m, d, y);
print('-dpdf', savestring)

%CA3
figure
imagesc(g(2).binx{10},g(2).biny{10},g(2).rate{10}')
axis xy
colorbar('ylim',[0 19],'ytick',19)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_CA3globalrate_example.png', m, d, y);
print('-dpng', savestring)

load '/data6/monster/Cyc/Cycpos06.mat'
load '/data6/monster/Cyc/Cycspikes06.mat'

figure
subplot(1,3,1)
plot(pos{6}{2}.data(:,2),pos{6}{2}.data(:,3),'k',spikes{6}{2}{5}{3}.data(:,2),spikes{6}{2}{5}{3}.data(:,3),'r.')
set(gca,'xlim',[min(g(2).binx{10}) max(g(2).binx{10})],'Ylim',[min(g(2).biny{10}) max(g(2).biny{10})])
subplot(1,3,2)
plot(pos{6}{4}.data(:,2),pos{6}{4}.data(:,3),'k',spikes{6}{4}{5}{3}.data(:,2),spikes{6}{4}{5}{3}.data(:,3),'r.')
set(gca,'xlim',[min(g(2).binx{10}) max(g(2).binx{10})],'Ylim',[min(g(2).biny{10}) max(g(2).biny{10})])
subplot(1,3,3)
plot(pos{6}{6}.data(:,2),pos{6}{6}.data(:,3),'k',spikes{6}{6}{5}{3}.data(:,2),spikes{6}{6}{5}{3}.data(:,3),'r.')
set(gca,'xlim',[min(g(2).binx{10}) max(g(2).binx{10})],'Ylim',[min(g(2).biny{10}) max(g(2).biny{10})])

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/NovelObjectPaper/%d_%d_%d_CA3placefield_example.pdf', m, d, y);
print('-dpdf', savestring)

