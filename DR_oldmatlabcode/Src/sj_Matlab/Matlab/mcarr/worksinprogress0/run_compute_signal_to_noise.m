%Run compute_signal_to_noise
stats = cell(3,1);

stats{1} = compute_signal_to_noise('/data13/mcarr/Bon/','bon',3:10);
stats{2} = compute_signal_to_noise('/data13/mcarr/Fra/','fra',2:12);
stats{3} = compute_signal_to_noise('/data13/mcarr/Ten/','ten',1:7);

save('/data13/mcarr/RipplePaper/signal_to_noise_stats.mat','stats');

time = -0.15:0.01:0.35;

%Plot coherence and phase variance for gamma power matched ripples and non-ripples
rip_std = []; norip_std = [];
rip_c = []; norip_c = [];
for an = [1 3]%:length(stats)
    for d = 1:length(stats{an})
        for e = 1:length(stats{an}{d})
            if ~isempty(stats{an}{d}{e}) & isfield(stats{an}{d}{e},'gam_power_rip_norip_pvalue')
                rip_std = [rip_std std(stats{an}{d}{e}.ripple_ca1_ca3_phase,[],2)];
                norip_std = [norip_std std(stats{an}{d}{e}.no_ripple_ca1_ca3_phase,[],2)];
                rip_c = [rip_c mean(stats{an}{d}{e}.ripple_ca1_ca3_coherence,2)];
                norip_c = [norip_c mean(stats{an}{d}{e}.no_ripple_ca1_ca3_coherence,2)];
            end
        end
    end
end

rip_baseline = repmat(mean(rip_c(1:6,:)),size(rip_c,1),1);
norip_baseline = repmat(mean(norip_c(1:6,:)),size(norip_c,1),1);

x = [time time(end:-1:1)];
mean_r = mean(rip_c-rip_baseline,2); se_r = std(rip_c-rip_baseline,[],2)./sqrt(size(rip_c,2)-1);
mean_nr = mean(norip_c-norip_baseline,2); se_nr = std(norip_c-norip_baseline,[],2)./sqrt(size(norip_c,2)-1);

figure
hold on
plot(time,mean_r,'r',time,mean_nr,'b')
fill(x,[mean_r + se_r; mean_r(end:-1:1) - se_r(end:-1:1)],'r','EdgeColor','none')
fill(x,[mean_nr + se_nr; mean_nr(end:-1:1) - se_nr(end:-1:1)],'b','EdgeColor','none')
set(gca,'xlim',[time(1) time(end)])
xlabel('Time since ripple detection (s)')
ylabel('delta CA1-CA3 coherence')

x = [time time(end:-1:1)];
mean_r = mean(rip_std,2); se_r = std(rip_std,[],2)./sqrt(size(rip_c,2)-1);
mean_nr = mean(norip_std,2); se_nr = std(norip_std,[],2)./sqrt(size(norip_c,2)-1);

figure
hold on
plot(time,mean_r,'r',time,mean_nr,'b')
fill(x,[mean_r + se_r; mean_r(end:-1:1) - se_r(end:-1:1)],'r','EdgeColor','none')
fill(x,[mean_nr + se_nr; mean_nr(end:-1:1) - se_nr(end:-1:1)],'b','EdgeColor','none')
set(gca,'xlim',[time(1) time(end)])
xlabel('Time since ripple detection (s)')
ylabel('Angular variance of phase coherence')