animprefix = 'dwi'
n = 1
clr = ['br']
CI = 1
taskswitch = 1 %1 to plot both tasks, switching period
switchday = 6 %first day of switch
titlename = [animprefix, ' New Algorithm'];

%----
task = ['Old'; 'New'];
for tk = 1:taskswitch+1
    color = clr(tk)

    eval(['load ', animprefix, 'estprobcorr_outbnd_', task(tk,:),'Tsk.mat'])
    eval(['load ', animprefix, task(tk,:),'estpc.mat'])
    eval(['totalbp = totalbpout', task(tk,:),';']);
    epochs = unique(totalbp(:,1:2), 'rows');
    endday = find(diff(epochs(:,1)));

    enddaytrial = [];
    for i = 1:size(epochs, 1)
        bp{epochs(i,1)}{epochs(i,2)} = totalbp(ismember(totalbp(:,1:2),epochs(i,:), 'rows'), 5);
        bptrial{epochs(i,1)}{epochs(i,2)} =  totalbp(ismember(totalbp(:,1:2),epochs(i,:), 'rows'), 3);
        %         if ismember(i, endday)
        %             enddaytrial = [enddaytrial; bptrial{epochs(i,1)}{epochs(i,2)}(end)];
        %         end
    end

    %     figure(n)
    %     hold on
    %     for i = 1:length(enddaytrial)
    %         plot([enddaytrial(i) enddaytrial(i)], [0 1], 'k', 'linewidth', 2)
    %     end

    for i = 1:size(epochs, 1)
        figure(n)
        hold on
        if ismember(i, endday)
            plot(bptrial{epochs(i,1)}{epochs(i,2)}(end)*[1 1],[ 0 1] , 'k', 'linewidth', 2)
        else
            plot(bptrial{epochs(i,1)}{epochs(i,2)}(end)*[1 1],[ 0 1] , 'k--', 'linewidth', 1)
        end
        if CI == 1
            plot(bptrial{epochs(i,1)}{epochs(i,2)}(:,1), estpc{epochs(i,1)}{epochs(i,2)}(:,2), [color, '--'], 'linewidth', 2)
            plot(bptrial{epochs(i,1)}{epochs(i,2)}(:,1), estpc{epochs(i,1)}{epochs(i,2)}(:,4), [color, '--'], 'linewidth', 2)
        end
        plot(bptrial{epochs(i,1)}{epochs(i,2)}(:,1), estpc{epochs(i,1)}{epochs(i,2)}(:,3), color, 'linewidth', 2)

    end
end

tmp = epochs(epochs(:,1) == switchday, 2);
switchep = tmp(1);
figure(n)
title(titlename)
if taskswitch == 1
    xlim([bptrial{switchday}{switchep}(1)  bptrial{epochs(end,1)}{epochs(end,2)}(end) ])
else
    xlim([ 0 bptrial{switchday}{switchep}(1)])

end
