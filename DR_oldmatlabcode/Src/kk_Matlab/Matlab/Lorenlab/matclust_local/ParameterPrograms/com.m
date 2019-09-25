function [cm,cms] = com(w)
% function cm = com(w)

L = size(w,1);
cms = ([1:L]*w ./ sum(w,1));
% cm = mean(cms);
sw = sum(w,2);
sw = sw - min(sw);
cm = [1:L]*sw / sum(sw);
