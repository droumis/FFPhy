%Run compute_gamma_ripple_glm

stats_ten = compute_gamma_ripple_glm('/data13/mcarr/Ten/','ten',1:7);
stats_bon = compute_gamma_ripple_glm('/data13/mcarr/Bon/','bon',3:10);
stats_fra = compute_gamma_ripple_glm('/data13/mcarr/Fra/','fra',2:12);

stats(1).animal = 'Ten';
stats(1).ca1_beta = stats_ten.ca1_beta;
stats(1).ca1_pvalue = stats_ten.ca1_pvalue;
stats(1).ca1_prediction = stats_ten.ca1_prediction;
stats(1).ca3_beta = stats_ten.ca3_beta;
stats(1).ca3_pvalue = stats_ten.ca3_pvalue;
stats(1).ca3_prediction = stats_ten.ca3_prediction;

stats(1).epoch = [0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 1 0 0 1 0 1 0 1 0,...
    0 1 0 1 0 1 0 0 1 0 1 0];

stats(2).animal = 'Bond';
stats(2).ca1_beta = stats_bon.ca1_beta;
stats(2).ca1_pvalue = stats_bon.ca1_pvalue;
stats(2).ca1_prediction = stats_bon.ca1_prediction;
stats(2).ca3_beta = stats_bon.ca3_beta;
stats(2).ca3_pvalue = stats_bon.ca3_pvalue;
stats(2).ca3_prediction = stats_bon.ca3_prediction;
stats(2).epoch = [0 1 0 1 0 1 0 0 1 0 1 0 1 0 0 1 0 1 0 1 0 0 1 0 1 0 1,...
    0 0 1 0 1 0 1 0 0 1 0 1 0 1 0 0 1 0 1 0 1 0 0 1 0 1 0 1 0];

stats(3).animal = 'Frank';
stats(3).ca1_beta = stats_fra.ca1_beta;
stats(3).ca1_pvalue = stats_fra.ca1_pvalue;
stats(3).ca1_prediction = stats_fra.ca1_prediction;
stats(3).ca3_beta = stats_fra.ca3_beta;
stats(3).ca3_pvalue = stats_fra.ca3_pvalue;
stats(3).ca3_prediction = stats_fra.ca3_prediction;
stats(3).epoch = [0 1 0 1 0 1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 0 1 0 1 0 1 0,...
    0 1 0 1 0 1 0 0 1 0 1 0 1 0 0 1 0 1 0 1 0 0 1 0 1 0 1 0 0 1 0 1 0 1 0,...
    0 1 0 1 0 1 0 0 1 0 1 0 1 0];


save('/data13/mcarr/RipplePaper/glmstats.mat','stats');

%Load and plot prediction of ripples based on gamma power
load('/data13/mcarr/RipplePaper/glmstats.mat')

y1 = []; y3 = [];
x = [-10:0.5:10 10:-0.5:-10];
for an = 1:length(stats)
    y1 = [y1 stats(an).ca1_prediction];
    y3 = [y3 stats(an).ca3_prediction];
end

%Average together the data greater than 10
tmp = y1;
tmp(41:end,:) = [];
tmp(41,:) = nanmean(y1(41:end,:));
y1 = tmp;

tmp = y3;
tmp(41:end,:) = [];
tmp(41,:) = nanmean(y3(41:end,:));
y3 = tmp; clear tmp

mean_1 = nanmean(y1,2)'; mean_3 = nanmean(y3,2)';
se_1 = nanstd(y1,[],2)'./sqrt(size(y1,2)-1);    
fill_1 = [mean_1+se_1 mean_1(end:-1:1)-se_1(end:-1:1)];
invalid = isnan(fill_1); x_1 = x(~invalid); fill_1(invalid) = [];
se_3 = nanstd(y3,[],2)'./sqrt(size(y3,2)-1);
fill_3 = [mean_3+se_3 mean_3(end:-1:1)-se_3(end:-1:1)];
invalid = isnan(fill_3); x_3 = x(~invalid); fill_3(invalid) = [];


figure
hold on
plot(-10:.5:10,mean_1,'b',-10:0.5:10,mean_3,'r')
legend([{'CA1'},{'CA3'}],'location','NorthWest')
fill(x_1,fill_1,'b','EdgeColor','None')
fill(x_3,fill_3,'r','EdgeColor','None')
set(gca,'xlim',[-2 11],'ylim',[0 1],'xtick',-1:1:10,'ytick',0:0.2:1,'yticklabel',0:20:100)
box off
xlabel('Gamma power (z-score)')
ylabel('Probability of observing an SWR')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_gamma_predicting_SWRs.pdf', m, d, y);
print('-dpdf', savestring)

p1 = []; p3 = [];
for an = 1:length(stats)
    p1 = [p1 stats(an).ca1_pvalue(2,:)];
    p3 = [p3 stats(an).ca3_pvalue(2,:)];
end
% 77% of sessions significant for CA3, 83% of sessions significant for CA1