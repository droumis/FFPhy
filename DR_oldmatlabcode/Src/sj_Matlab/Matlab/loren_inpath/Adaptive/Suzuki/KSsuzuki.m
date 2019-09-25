orient tall
subplot(4,2,1);
getKSPlot(CKSStat);
title('Cortical');
set(gca, 'XTick', [], 'YTick', []);
subplot(4,2,2);
getKSPlot(AKSStat);
title('Hippocampal');
set(gca, 'XTick', [], 'YTick', []);

