
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48','S58','S59','S60','S61'};
threshold = {0.006; 0.010}; % 6 ms and 10 ms thresholds for bursts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)

  clear('unit');
  % normalized_power
  files_list = dir(sprintf('%s/%s/spike/*_unit.mat',path_prefix,subjects{s}));
  for i = 1:numel(files_list)
    filename = sprintf('%s/%s/spike/%s', ...
      path_prefix,subjects{s},files_list(i).name);
    disp(filename);
    load(filename);
    for j = 1:numel(unit)
      [isi_pre, isi_post, burstflag{1,1}] = unit_isi(unit(j),threshold{1});
      [isi_pre, isi_post, burstflag{2,1}] = unit_isi(unit(j),threshold{2});
      unit(j).isi_pre = isi_pre;
      unit(j).isi_post = isi_post;
      unit(j).burst = struct('threshold',threshold,'flag',burstflag);
      clear('isi_pre','isi_post','burstflag');
    end
    disp(unit);
    save(filename,'unit');
  end
end
