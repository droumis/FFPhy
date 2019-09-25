%% Calculate place fields for slow speeds
% Animal Selection
animals = {'Conley','Corriander','Dudley','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilter = [];
for i = 1:10
    epochfilter{i} = ['($exposure == ',num2str(i),') & isequal($description,''TrackA'')'];
end

% Time Filter
% timefilter = {{'getriptimes', '($nripples== 0)', [], 'cellfilter', 'isequal($area, ''CA1'')'},...
%     {'get2dstate', 'abs($velocity) < 5'}};
timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', 'isequal($area, ''CA1'')'},...
    {'getgammatimes', '($ngamma >0)', [], 'low','tetfilter', 'isequal($area,''CA1'')','minthresh',3,'inclusive',1,'min_separation',0.5}};

% Cell Filter
cellfilter = 'isequal($area, ''CA1'') & $numspikes > 200 & $peakrate>=5 & $meanrate < 7';
iterator = 'singlecellanal';

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);

f = setfilterfunction(f, 'filtercalclinfields', {'spikes','linpos'},'binsize',2.5);

f = runfilter(f);

slow = f;

%% Calculate place fields for fast times

% Time Filter
%timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', 'isequal($area, ''CA1'')'},...
%    {'get2dstate', 'abs($velocity) > 5'}};
timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', 'isequal($area, ''CA1'')'},...
    {'getgammatimes', '($ngamma >0)', [], 'high','tetfilter', 'isequal($area,''CA1'')','minthresh',3,'inclusive',1,'min_separation',0.5}};

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);

f = setfilterfunction(f, 'filtercalclinfields', {'spikes','linpos'},'binsize',2.5);

f = runfilter(f);

fast = f;

%% Calculate place fields for all times

% Time Filter
timefilter = {{'getriptimes', '($nripples == 0)', [], 'cellfilter', 'isequal($area, ''CA1'')'}};

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter, 'iterator', iterator);

f = setfilterfunction(f, 'filtercalclinfields', {'spikes','linpos'},'binsize',2.5);

f = runfilter(f);

all = f;

%% Plot examples
for an = 1:length(all)
    for d = 1:length(all(an).output)
        for e = 1:length(all(an).output{d})
            figure(1)
            max_rate = 0;
            for traj = 1:length(all(an).output{d}(e).trajdata)
                max_rate = max(max_rate,max(all(an).output{d}(e).trajdata{traj}(:,5))); 
            end
            for traj = 1:length(all(an).output{d}(e).trajdata)
                subplot(length(all(an).output{d}(e).trajdata),1,traj)
                hold on
                plot(all(an).output{d}(e).trajdata{traj}(:,1),all(an).output{d}(e).trajdata{traj}(:,5),'k')
                plot(fast(an).output{d}(e).trajdata{traj}(:,1),fast(an).output{d}(e).trajdata{traj}(:,5),'g')
                plot(slow(an).output{d}(e).trajdata{traj}(:,1),slow(an).output{d}(e).trajdata{traj}(:,5),'r')
                set(gca,'ylim',[0 max_rate*2]);
            end
            pause(5)
            clf
        end
    end
end

%% Compute deviation from over all place field
deviation = cell(length(all),1);
for an = 1:length(all)
    deviation{an} = cell(length(all(an).output),1);
    for d = 1:length(all(an).output)
        deviation{an}{d} = nan(length(all(an).output{d}),2);
        for e = 1:length(all(an).output{d})
            tmp = nan(4,3);
            for traj = 1:length(all(an).output{d}(e).trajdata)
                invalid = isnan(all(an).output{d}(e).trajdata{traj}(:,5)+...
                        slow(an).output{d}(e).trajdata{traj}(:,5) +...
                        fast(an).output{d}(e).trajdata{traj}(:,5) );
                    
                tmp(traj,1) = trapz(all(an).output{d}(e).trajdata{traj}(~invalid,1),...
                    all(an).output{d}(e).trajdata{traj}(~invalid,5));
                tmp(traj,2) = trapz(slow(an).output{d}(e).trajdata{traj}(~invalid,1),...
                    slow(an).output{d}(e).trajdata{traj}(~invalid,5));
                tmp(traj,3) = trapz(fast(an).output{d}(e).trajdata{traj}(~invalid,1),...
                    fast(an).output{d}(e).trajdata{traj}(~invalid,5)); 
            end
            deviation{an}{d}(e,1) = nansum(tmp(:,2))/nansum(tmp(:,1));
            deviation{an}{d}(e,2) = nansum(tmp(:,3))/nansum(tmp(:,1));
        end
    end
end

%Now combine across animals
day_deviation = cell(10,1);
for an = 1:length(all)
    for d = 1:length(all(an).output)
        day_deviation{d} = stack(day_deviation{d},deviation{an}{d});
    end
end

bin = 0:0.05:2;
slow_deviation = zeros(length(day_deviation),length(bin));
fast_deviation = zeros(length(day_deviation),length(bin));
for d = 1:length(day_deviation)
    slow_deviation(d,:) = hist(day_deviation{d}(:,1),bin)./sum(~isnan(day_deviation{d}(:,1)));
    fast_deviation(d,:) = hist(day_deviation{d}(:,2),bin)./sum(~isnan(day_deviation{d}(:,2)));
end

%% Compute overlap of place fields
overlap = cell(10,1);
count = cell(10,1);
for an = 1:length(all)
    for d = 1:length(all(an).output)
        if isempty(count{d})
            count{d} = 1;
            overlap{d} = nan(100,1);
        end
        for e = 1:length(all(an).output{d})
            overlap{d}(count{d}) = calcoverlap(fast(an).output{d}(e).trajdata,...
                slow(an).output{d}(e).trajdata);
            count{d} = count{d}+1;
        end
    end
end
for d = 1:length(overlap)
    overlap{d}(isnan(overlap{d})) = [];
end