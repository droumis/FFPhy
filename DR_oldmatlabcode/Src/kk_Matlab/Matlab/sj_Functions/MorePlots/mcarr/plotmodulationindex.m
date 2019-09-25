%% Plot Modulation Index Example.

load '/data13/mcarr/Fiv/Fivmodindex.mat'
g = gaussian2(10,2);

for i = 1:24;
temp = [];
for j = 1:length(modindex.output{1}(i).modulationindex)
    if j==1
        temp = [conv2(abs(modindex.output{1}(i).modulationindex{j}(:,3:17)),g,'same')];
    else
        temp = [temp + conv2(abs(modindex.output{1}(i).modulationindex{j}(:,3:17)),g,'same')];
    end
end

temp = temp./length(modindex.output{1}(i).modulationindex);

imagesc([4:1:18],[5:5:200],temp)
axis xy
colorbar
set(gca,'fontsize',18)
xlabel('Frequency for Phase (Hz)','fontsize',18)
ylabel('Frequency for Amplitude (Hz)','fontsize',18)
title('Modulation Index','fontsize',18)
pause
end

%% Combine across segments and look at theta modulation over days
load '/data13/mcarr/Fiv/Fivmodindex.mat'
g = gaussian2(15,5);
expose1 = [1 2 3 4 5 6 7 10 13 16];

lincohere = [];
epoch = modindex.epochs{1};
for k= 1:length(expose1)
    i= expose1(k);
    temp = [];
    for j = 1:length(modindex.output{1}(i).modulationindex)
        if j==1
            temp = mean(conv2(abs(modindex.output{1}(i).modulationindex{j}(:,3:17)),g,'same'),2);
        else
            temp = [temp mean(conv2(abs(modindex.output{1}(i).modulationindex{j}(:,3:17)),g,'same'),2)];
        end
    end
    %temp = temp./length(modindex.output{1}(i).modulationindex);
    lincohere = [lincohere temp];
end

figure
imagesc(1:length(lincohere),[5:5:200],lincohere)
axis xy
colorbar
set(gca,'fontsize',18)
