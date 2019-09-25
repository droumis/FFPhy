



% check theta peak of each reference channel

params.Fs = continuous(1).Fs;
params.fpass = [0 20];


for i = 1:numel(continuous)
  subplot(1,i,numel(continuous));


