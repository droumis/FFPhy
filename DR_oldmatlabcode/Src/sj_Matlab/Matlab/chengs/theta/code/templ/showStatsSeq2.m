s=[];
av=[];
va=[];
n=[];
for i=1:16
    s{i}=collectStats('r',i);
    n(i)= length(s{i}.val);
    av(i)= mean(s{i}.val);
    va(i)= var(s{i}.val);
end


figure(1); clf

subplot(3,1,1);
errorbar(1:16, av, sqrt(va)./sqrt(n));
title([s{1}.name ': mean with std dev'])
axis tight

subplot(3,1,2);
errorbar(1:16, av, sqrt(va));
title([s{1}.name ': mean with std error'])
axis tight

subplot(3,1,3);
%plot(s.val');
% plot(s.val','.','MarkerSize',2);
% title([s.name ': scatter'])
% axis tight

figure(2); clf

subplot(3,1,1);
hist(s{1}.val)
title([s{1}.name ': at 1'])
axis tight

subplot(3,1,2);
hist(s{8}.val)
title([s{1}.name ': at 8'])

axis tight
subplot(3,1,3);
hist(s{16}.val)
title([s{1}.name ': at 16'])
axis tight

