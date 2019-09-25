% Animal selection
animals = {'Bond','Ten','Frank'};

% Epoch selection
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';
epochfilter{2} = 'isequal($type,''sleep'')';

% Tetrode pair selection
tetfilter = {'isequal($area, ''CA1'') & isequal($hemisphere,''left'')','isequal($area, ''CA1'') & isequal($hemisphere,''right'')'};

%Iterator Selection
iterator = 'epocheegnonreferenceanal';

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter, 'eegtetrodepairs',tetfilter,'iterator', iterator);
f = setfilterfunction(f, 'calcripplecoherence', {'ripples','cellinfo','ripple','lowgamma'});
f = runfilter(f);

% save('/data13/mcarr/RipplePaper/ripplecoherence.mat','f')
rip = []; gam = [];
for an = [1 3]
    for d = 1%:length(f(an).output)
        ind = length(f(an).output{d}(1).ripple_corr);
        tmp_rip = []; tmp_gam = [];
        for e = 1:length(f(an).output{d})
            if length(f(an).output{d}(e).ripple_corr) > 10
                if ind == length(f(an).output{d}(e).ripple_corr)
                    tmp_rip = [tmp_rip; f(an).output{d}(e).ripple_corr];
                    tmp_gam = [tmp_gam; f(an).output{d}(e).gamma_corr];
                else
                    if size(tmp_rip,1)>2
                        rip = [rip mean(tmp_rip)];
                        gam = [gam mean(tmp_gam)];
                    end
                    tmp_rip = f(an).output{d}(e).ripple_corr;
                    tmp_gam = f(an).output{d}(e).gamma_corr;
                end
            end
            ind = length(f(an).output{d}(e).ripple_corr);
        end
    end
end

figure
bar(-1:0.05:1,hist(rip,-1:0.05:1)./length(rip),'b')
set(gca,'xlim',[-1 1],'xtick',-1:0.2:1,'ylim',[0 0.11],'ytick',0:0.05:0.15)
xlabel('Correlation between CA1 ripple oscillations across hemispheres')
ylabel('Proportion')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_rippleoscillation_correlation.pdf', m, d, y);
print('-dpdf', savestring)
%mean: 0.0212, %median: 0.02

figure
bar(-1:0.05:1.1,hist(gam,-1:0.05:1.1)./length(gam),'b')
set(gca,'xlim',[-1 1],'xtick',-1:0.2:1,'ylim',[0 0.11],'ytick',0:0.05:0.15)
xlabel('Correlation between CA1 gamma oscillations across hemispheres')
ylabel('Proportion')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_gammaoscillation_correlation.pdf', m, d, y);
print('-dpdf', savestring)