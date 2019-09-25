
subjects = {'S48','S58','S59','S60','S61'};

for i = 1:numel(subjects)
  subject = subjects{i};

  % load spike_config  
  load(sprintf('/data14/smkim/%s/spike/%s_spike_config.mat',subject,subject));
  
  tetfolders = dir(sprintf('/data14/smkim/%s/spike/tetrode*',subject));
  for j = 1:length(tetfolders)
    tetfolder = tetfolders(j).name;
    spike_files = dir(sprintf('/data14/smkim/%s/spike/%s/*_spike.mat',subject,tetfolder));
    for k = 1:length(spike_files)
      spike_file = sprintf('/data14/smkim/%s/spike/%s/%s',subject,tetfolder, ...
          spike_files(k).name);
      try
        spike = load(spike_file);
      catch
        warning('could not load %s',spike_file);
        continue;
      end
      
      %spike.source = [spike_file(1:(regexp(spike_file,'_spike.mat','once')-1)) '.tt']

      match_idx = find( ...
          ([spike_config(:).day] == spike.day) & ...
          ([spike_config(:).tetrode] == spike.tetrode) );

      spike.sources = { spike.source, spike_config(match_idx).sources{:} };
      spike = rmfield(spike,'source');
      save(spike_file,'-struct','spike');
      clear('spike','match_idx');

    end
  end
end

