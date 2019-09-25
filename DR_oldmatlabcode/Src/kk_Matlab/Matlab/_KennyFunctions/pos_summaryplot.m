%% plots velocity, epoch-by-epoch, for a given day

%% ingredient: pos


figure
hold on
d=2;     % select day

for e=1:length(pos{d})
    subplot(length(pos{d}),1,e)
    plot(pos{d}{e}.data(:,5));
    ylim([0 60])
end