%% CALCULATE FIRING RATES

% Animal Selection
animals = {'CorrianderNO','Cyclops','Dunphy','Fafnir','Godzilla','Grendel','Cml','Nico'};

% Epoch selection
epochfilter{1} = 'isequal($session,''familiar'')';
epochfilter{2} = 'isequal($session,''novel'')';
epochfilter{3} = 'isequal($session,''supernovel'')';

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

ca1 = calcopenfieldratestats(f);
ca3 = calcopenfieldratestats(f);

%% Save filters
%save('/data6/monster/analysis/ca1_globalopenfieldrates.mat','f')
%save('/data6/monster/analysis/ca1_residuals.mat','g')

%save('/data6/monster/analysis/ca3_globalopenfieldrates.mat','f')
%save('/data6/monster/analysis/ca3_residuals.mat','g')

%% Initialize variables
m = []; c = []; lf = []; ln = []; ls = []; peak = [];

for an = 1:length(ca1)
    m = [m; ca1(an).mean];
    c = [c; ca1(an).corr];
    peak = [peak; ca1(an).peak];
    for d = 1:length(ca1(an).location)
        if ~isempty(ca1(an).location{d})
            lf = [lf; ca1(an).location{d}(1,:)];
            ln = [ln; ca1(an).location{d}(2,:)];
            ls = [ls; ca1(an).location{d}(3,:)];
        else
            lf = [lf; [NaN NaN]];
            ln = [ln; [NaN NaN]];
            ls = [ls; [NaN NaN]];
        end
    end
end

%% Look at mean rate, novel and supernovel
nov = m(:,1); sup = m(:,2); nloc = ln; sloc = ls;
ln(isnan(nov),:) = []; nov(isnan(nov)) = []; ls(isnan(sup),:) = []; sup(isnan(sup)) = [];

figure
hold on
bar(1:2,[mean(nov) mean(sup)],'b')
errorbar(1:2,[mean(nov) mean(sup)],[stderr(nov) stderr(sup)],'k')
set(gca,'xtick',1:2,'xticklabel',[{'Nov'},{'Sup'}])

%Run statistics
signrank(nov,1); %p<0.05
signrank(sup,1); %p<0.01

%% Look at mean rate by distance to objects

%% Look at corrleation

%% Look at correlation by distance to objects


