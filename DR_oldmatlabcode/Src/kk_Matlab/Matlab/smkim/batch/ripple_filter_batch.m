
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48', 'S58', 'S59', 'S60', 'S61'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)

  files_list = dir(sprintf('%s/%s/continuous/*_lfp.mat',path_prefix,subjects{s}));

  for i = 1:numel(files_list)
    load(sprintf('%s/%s/continuous/%s',path_prefix,subjects{s},files_list(i).name));

    ripple = filter_continuous(lfp,'ripple_filter.mat') 
    save(sprintf('%s/%s/continuous/%s',path_prefix,subjects{s}, ...
        [files_list(i).name(1:end-8) '_ripple.mat']),'ripple');
    clear('lfp','ripple');
  end

end

