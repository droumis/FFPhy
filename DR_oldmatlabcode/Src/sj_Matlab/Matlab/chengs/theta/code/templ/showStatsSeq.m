%s=collectStatsSeq('r',16);

figure(1); clf

subplot(3,1,1);
errorbar(1:16, mean(s.val), std(s.val)/sqrt(size(s.val,1)));
title([s.name ': mean with std dev'])
axis tight

subplot(3,1,2);
errorbar(1:16, mean(s.val), std(s.val));
title([s.name ': mean with std error'])
axis tight

subplot(3,1,3);
%plot(s.val');
plot(s.val','.','MarkerSize',2);
title([s.name ': scatter'])
axis tight

figure(2); clf

subplot(3,1,1);
hist(s.val(:,1))
title([s.name ': at 1'])
axis tight

subplot(3,1,2);
hist(s.val(:,8))
title([s.name ': at 8'])

axis tight
subplot(3,1,3);
hist(s.val(:,16))
title([s.name ': at 16'])
axis tight

