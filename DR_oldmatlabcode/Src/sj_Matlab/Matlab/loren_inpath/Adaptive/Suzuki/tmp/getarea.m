ntotaltrials = length(dataset{6}.cobj.cond_no);

for i = 1:ntotaltrials
    dataset{6}.cell{1}.fixrate(i) = getspikerate(CP.x, dataset{6}.cell{1}.trialtimerate(i,:), 0, 300);
    dataset{6}.cell{1}.stimrate(i) = getspikerate(CP.x, dataset{6}.cell{1}.trialtimerate(i,:), 300, 800);
    dataset{6}.cell{1}.delayrate(i) = getspikerate(CP.x, dataset{6}.cell{1}.trialtimerate(i,:), 800, 1500);
    dataset{6}.cell{1}.resprate(i) = getspikerate(CP.x, dataset{6}.cell{1}.trialtimerate(i,:), 1500, 2200);
end

