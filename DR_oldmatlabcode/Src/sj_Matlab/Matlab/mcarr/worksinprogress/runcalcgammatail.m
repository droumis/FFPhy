%% FOR CA1

% Animal selection
animals = {'Bond','Frank','Ten'};

% Epoch selection
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';
%epochfilter{2} = 'isequal($type,''sleep'')';

% Tetrode selection
tetfilter = 'isequal($area, ''CA1'')& $maxcell==1';

%Iterator Selection
iterator = 'epocheegnonreferenceanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter, 'eegtetrodes',tetfilter,'iterator', iterator,'excludetime',timefilter);
f = setfilterfunction(f, 'calcgammatail', {'lowgamma','ripple','ripples','cellinfo'});
f = runfilter(f);

%% PLOT DISTRIBUTION OF GAMMA RETURNING TO BASELINE

g = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        g = [g; f(an).output{d}(:)];
    end
end
g = g(~isnan(g));

figure
bar(-0.1:0.005:0.1,hist(g,-0.1:0.005:0.1)./length(g))
set(gca,'xlim',[-0.1 0.1],'ylim',[0 0.1])
ylabel('Proportion')
xlabel('Time that gamma power returns to baseline relative to end of SWR (s)')

%On average, gamma returns to baseline 1ms after ripple power returns to
%baseline. This effect was not significantly different for all animals
%by signrank (p<0.001, p<0.05, p>0.1) or ttest.

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_distribution_gamma_offset.pdf', m, d, y);
print('-dpdf', savestring)

g = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        g = [g; ttest(f(an).output{d}(:))];
    end
end

%% FOR CA3
% Animal selection
animals = {'Bond','Frank','Ten'};

% Epoch selection
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';
%epochfilter{2} = 'isequal($type,''sleep'')';

% Tetrode selection
tetfilter = {'isequal($area, ''CA3'') & $maxcell==1','isequal($area, ''CA1'')&$maxcell==1'};

%Iterator Selection
iterator = 'epocheegnonreferenceanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter, 'eegtetrodepairs',tetfilter,'iterator', iterator,'excludetime',timefilter);
f = setfilterfunction(f, 'calcgammatail', {'lowgamma','ripple','ripples','cellinfo'});
f = runfilter(f);

%% PLOT DISTRIBUTION OF GAMMA RETURNING TO BASELINE

g = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        g = [g; f(an).output{d}(:)];
    end
end
g = g(~isnan(g));

figure
bar(-0.1:0.005:0.1,hist(g,-0.1:0.005:0.1)./length(g))
set(gca,'xlim',[-0.1 0.1],'ylim',[0 0.1])
ylabel('Proportion')
xlabel('Time that gamma power returns to baseline relative to end of SWR (s)')

%On average, gamma returns to baseline 0ms after ripple power returns to
%baseline. This effect was not significantly different by signrank 
%(p>0.2) or ttest (p>0.3).

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_CA3_distribution_gamma_offset.pdf', m, d, y);
print('-dpdf', savestring)

g = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        g = [g; ttest(f(an).output{d}(:))];
    end
end
