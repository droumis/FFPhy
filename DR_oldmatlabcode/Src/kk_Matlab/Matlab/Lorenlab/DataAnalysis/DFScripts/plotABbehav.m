function [h] = plotABbehav(bp)
% function [h] = plotABbehav(behavperform)
%
% plots the estimated probability correct and confidence bounds for inbound and
% outbound runs for Track A and B W-tracks. 
%
% behavperform is the structure created by createbehavperform

% figure out which track was presented first
adayfirst = bp(1).dayepoch(1,1) <= bp(2).dayepoch(1,1);
aepochfirst = bp(1).dayepoch(1,2) < bp(2).dayepoch(1,2);

if (adayfirst) 
    afirst = 1;
elseif (






figure
subplot(2,1,1);
plot(


