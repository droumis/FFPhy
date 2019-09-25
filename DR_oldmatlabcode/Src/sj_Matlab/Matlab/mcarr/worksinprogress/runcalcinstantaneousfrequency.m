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
f = setfilterfunction(f, 'calcinstantaneousfrequency', {'eeg','ripples','cellinfo'});
f = runfilter(f);

%% PLOT distribution of frequencies observed for SWRs

g = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).rip)
                a = accumarray(f(an).output{d}(e).rip, f(an).output{d}(e).freq,[],@(x) mean(x),NaN);
                g = [g; a]; 
            end
        end
    end
end
g = g(~isnan(g));

figure
hold on
subs = 10:4:50;
val = lookup(g,subs,1);
boxplot(g,'colors','k','orientation','horizontal','notch','on','whisker',0,'symbol','k')
plot(subs,hist(val,1:length(subs))./length(val))
set(gca,'xlim',[10 50],'ylim',[0 1.1],'ytick',0:0.1:0.5,'ytickLabel',0:0.1:0.5)
xlabel('Frequencies observed during SWRs in CA1')
ylabel('Proportion')
box off


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_distribution_frequencies_during_SWRs.pdf', m, d, y);
print('-dpdf', savestring)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_CA3_distribution_frequencies_during_SWRs.pdf', m, d, y);
print('-dpdf', savestring)

% CA1: 2% of SWRs have frequencies that are best described as <20Hz, the mean
% SWR frequency is 29.6 +- 0.17 and there is no evidence of a bimodal
% distribution. It seems that during SWRs there is a single frequency
% closer to 30Hz rather than two sepearable 20Hz and 50Hz oscillations
% median frequency = 29.3Hz, 25th prctile = 25.6Hz, 75th prctile = 33Hz

% CA3: 1.5% of SWRs have frequencies that are best described as <20Hz, the mean
% SWR frequency is 29.3 +- 0.16 and there is no evidence of a bimodal
% distribution. It seems that during SWRs there is a single frequency
% closer to 30Hz rather than two sepearable 20Hz and 50Hz oscillations
% median frequency = 29.1Hz, 25th prctile = 26Hz, 75th prctile = 33Hz