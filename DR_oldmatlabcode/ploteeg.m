function ploteeg
ep = 1
tett = 8

JK_ripple_ms


eegeight = eeg{1,25}{1,ep}{1,tett}.data;
lines = ripples{1, 25}{1,ep}{1, tett}.endind;
figure
plot(eegeight, 'Color', 'c')
for i = 1:length(lines)
line([lines(i) lines(i)],[-500 500], 'Color', 'g')
end
% 
% JK_ripple_ms(ep,9);
% eegnine = eeg{1,25}{1,ep}{1,9}.data;
% linesnine = ripples{1, 25}{1,ep}{1, 9}.endind;
% figure 2
% plot(eegnine, 'Color', 'b')
% for i = 1:length(linesnine)
% line([lines(i) lines(i)],[-1500 1500], 'Color', 'r')
% end
