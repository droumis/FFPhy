task = ['New'] %'Old' or 'New'
animprefix = 'cal'
dir = 'in'

%----
eval(['load ', animprefix, 'estprobcorr_', dir, 'bnd_', task,'Tsk.mat'])
eval(['totalbp = totalbp',dir, task,';']);
epochs = unique(totalbp(:,1:2), 'rows');

for i = 1:size(epochs, 1)
    d = epochs(i,1); e = epochs(i,2);
    bp{d}{e} = totalbp(ismember(totalbp(:,1:2),epochs(i,:), 'rows'), 5);
    bptrial{d}{e} =  totalbp(ismember(totalbp(:,1:2),epochs(i,:), 'rows'), 3);
end

i = 1
behavperf = bp{epochs(i,1)}{epochs(i,2)};
estpc{epochs(i,1)}{epochs(i,2)} = getestprobcorrectnew(behavperf, 0)

for i = 2:size(epochs, 1)
    behavperf = bp{epochs(i,1)}{epochs(i,2)};
    if length(behavperf) > 1
        estpc{epochs(i,1)}{epochs(i,2)} = getestprobcorrectnew(behavperf, 1)
    end
end


eval(['save ',animprefix,task, 'estpc', dir, ' estpc'])