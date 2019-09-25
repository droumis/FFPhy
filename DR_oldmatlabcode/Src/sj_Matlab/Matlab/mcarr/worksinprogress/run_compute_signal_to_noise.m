%Run compute_signal_to_noise
%stats = cell(3,1);

%stats{1} = compute_signal_to_noise('/data13/mcarr/Bon/','bon',3:10);
%stats{2} = compute_signal_to_noise('/data13/mcarr/Fra/','fra',2:12);
%stats{3} = compute_signal_to_noise('/data13/mcarr/Ten/','ten',1:7);

%save('/data13/mcarr/RipplePaper/signal_to_noise_stats.mat','stats');

load '/data13/mcarr/RipplePaper/ripstats.mat'
baseline = [];
for an = 1:length(f)
    for d =1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).ca1_ca1_coherence)
                baseline = [baseline mean(f(an).output{d}(e).ca1_ca3_coherence(1:6,:))];
            end
        end
    end
end
baseline = mean(baseline);
clear an d e f
load '/data13/mcarr/RipplePaper/ripple_phaselocking.mat'
baseline_p = [];
for an =1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).ca1_ca1_phase_std) && f(an).output{d}(e).nrips>5
                baseline_p = [baseline_p nanmean(f(an).output{d}(e).ca1_ca3_phase_std(1:6,2))];
            end
        end
    end
end
baseline_p = mean(baseline_p);
clear an d e f

load '/data13/mcarr/RipplePaper/signal_to_noise_stats.mat'
rip_c = []; norip_c = [];
rip_phase = []; norip_phase = [];
for an = 1:length(stats)
    for d = 1:length(stats{an})
        for e = 2:2:length(stats{an}{d})
            if ~isempty(stats{an}{d}{e}) & isfield(stats{an}{d}{e},'gam_power_rip_norip_pvalue')
                if ~isempty(stats{an}{d}{e}.ripple_ca1_ca3_coherence)
                    [ripval ripind] = max(stats{an}{d}{e}.ripple_ca1_ca3_coherence);
                    rip_c = [rip_c ripval];

                    [noripval noripind] = max(stats{an}{d}{e}.no_ripple_ca1_ca3_coherence);
                    norip_c = [norip_c noripval];
                    if ~isempty(stats{an}{d}{e}.ripple_ca1_ca3_phase)
                        ripval = zeros(size(ripind)); noripval = zeros(size(noripind));
                        for i = 1:length(ripind)
                            ripval(i)= stats{an}{d}{e}.ripple_ca1_ca3_phase(ripind(i),i);
                        end
                        for i = 1:length(noripind)
                            noripval(i)= stats{an}{d}{e}.no_ripple_ca1_ca3_phase(noripind(i),i);
                        end
                        [m r] = anglemean(ripval);
                        rip_phase = [rip_phase; r];

                        [m r] = anglemean(noripval);
                        norip_phase = [norip_phase; r];
                    end
                end
            end
        end
    end
end
clear an d e i m noripind noripval r ripind ripval

means = [mean(rip_c)-baseline mean(norip_c)-baseline mean(rip_phase)-baseline_p mean(norip_phase)-baseline_p];
errors = [std(rip_c)./sqrt(length(rip_c)-1) std(norip_c)./sqrt(length(norip_c)-1) std(rip_phase)./sqrt(length(rip_phase)-1) std(norip_phase)./sqrt(length(norip_phase)-1)];

figure
hold on
bar([1 2 4 5],means,1,'b')
errorbar2([1 2 4 5],means, errors,'k')
set(gca,'xtick',[1 2 4 5],'ylim',[-0.02 0.2],'ytick',0:0.05:1,'xticklabel',[{'SWR coherence'},{'No SWR coherence'},{'SWR phase locking'},{'No SWR phase locking'}])
xlabel('Time since ripple detection (s)')
ylabel('CA1-CA3 coherence')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_signaltonoise_ca1ca3coherence.pdf', m, d, y);
print('-dpdf', savestring)

