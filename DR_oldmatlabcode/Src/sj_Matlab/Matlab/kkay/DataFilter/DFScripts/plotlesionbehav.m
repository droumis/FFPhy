function [h] = plotABbehav(bp)
% function [h] = plotABbehav(behavperform)
%
% plots the estimated probability correct and confidence bounds for inbound and
% outbound runs for A W-track. 
%
% behavperform is the structure created by createbehavperform

% figure out which track was presented first
%adayfirst = bp(1).dayepoch(1,1) < bp(2).dayepoch(1,1);
%absameday = bp(1).dayepoch(1,1) = bp(2).dayepoch(1,1);
%aepochfirst = bp(1).dayepoch(1,2) < bp(2).dayepoch(1,2);


figure
subplot(2,1,1);
plot(bp(1).inprobcorrect(2:end, 1), 'b', 'LineWidth', 1.5);
hold on
plot(bp(1).inprobcorrect(2:end, 2), 'r', 'LineWidth', 0.5);
plot(bp(1).inprobcorrect(2:end, 3), 'r', 'LineWidth', 0.5);
for i = 1:size(bp(1).dayintrials,1)
    % plot the trial number at the end of the day
    if (isfinite(bp(1).dayintrials(i,2))) 
	plot([bp(1).dayintrials(i,2) bp(1).dayintrials(i,2)], [0 1], 'k'); 
    end
    day = sprintf('%d %d', bp(1).dayepoch(i,:));
    t = text(bp(1).dayintrials(i,1), .25, day);
    set(t, 'FontSize', 8);
end
plot(get(gca, 'XLim'), [.5 .5], 'k');
title(gca, 'Inbound')

subplot(2,1,2);
plot(bp(1).outprobcorrect(2:end, 1), 'b', 'LineWidth', 1.5);
hold on
plot(bp(1).outprobcorrect(2:end, 2), 'r', 'LineWidth', 0.5);
plot(bp(1).outprobcorrect(2:end, 3), 'r', 'LineWidth', 0.5);
for i = 1:size(bp(1).dayouttrials,1)
    % plot the trial number at the end of the day
    if (isfinite(bp(1).dayouttrials(i,2))) 
	plot([bp(1).dayouttrials(i,2) bp(1).dayouttrials(i,2)], [0 1], 'k'); 
    end
    day = sprintf('%d %d', bp(1).dayepoch(i,:));
    t = text(bp(1).dayouttrials(i,1), .25, day);
    set(t, 'FontSize', 8);
end
plot(get(gca, 'XLim'), [.5 .5], 'k');
title(gca, 'Outbound')
supertitle(sprintf('%s %s', bp(1).subject, bp(1).group));
