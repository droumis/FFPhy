function rp = plotmeanse(data1, data2, label1, label2,labely, varargin)
%plotmeanse(data1, data2, label1, label2,labely, varargin)
%data1 and data2 are data to be plotted
%label1 and label2 are names of groups

data1 = data1(~isnan(data1));
data2 = data2(~isnan(data2));

figure(gca)
subplot(1,2,2);
hold on
bar([1 2], [mean(data1) mean(data2)], 0.8)
errorbar2([1 2], [mean(data1) mean(data2)],  [stderr(data1) stderr(data2)] , 0.3, 'k') %'linewidth', 3)
xlim([0.3 2.7])
set(gca, 'fontsize', 24)
if ~isempty(label1) & ~isempty(label2)
    set(gca, 'xtick', [1 2], 'xticklabel', {label1, label2})
end
if ~isempty( labely)
    ylabel(labely)
end
rp = ranksum(data1,data2);


end
