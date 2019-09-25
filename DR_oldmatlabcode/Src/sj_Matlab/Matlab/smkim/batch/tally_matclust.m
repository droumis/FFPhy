
clear('tally');

tally = struct( ...
  'uid' , {}, ...
  'subject'             , {}                              , ...
  'day'                 , {}                                  , ...
  'tetrode'             , {}                                   , ...
  'region'              , {}                                    , ...
  'epoch'               , {}                                , ...
  'environment',   {}, ...
  'timerange'           , {}                            , ...
  'overall_mean_rate'   , {} , ...
  'complex_spike_index' , {}        , ...
  'timestamp'           , {}        , ...
  'samples'             , {}        , ...
  'shape'               , {}                                     );

subjects = {'S48','S58','S59','S60','S61'};

for i = 1:numel(subjects)
  subject = subjects{i};
  tetfolders = dir(sprintf('/data14/smkim/%s/spike/tetrode*',subject));
  for j = 1:length(tetfolders)
    tetfolder = tetfolders(j).name;
    matclust_files = dir(sprintf('/data14/smkim/%s/spike/%s/*_matclust.mat',subject,tetfolder));
    for k = 1:length(matclust_files)
      matclust_file = sprintf('/data14/smkim/%s/spike/%s/%s',subject,tetfolder, ...
          matclust_files(k).name);
      disp(matclust_file);
      try
        clust = load(matclust_file);
      catch
        error('could not load file %s',matclust_file);
        continue;
      end
      try
        spike = load(clust.clustattrib.datafile);
      catch
        error('could not load spike waveforms %s',clust.clustattrib.datafile);
      end

      % match sessions to matclust time fitlers
      session = clust.clustdata.customvar.session;
      for l = 1:numel(session)
        % Find all timefilterranges that fall within the session timerange
        m = find((clust.clustdata.timefilterranges(:,1) >= ...
          session(l).timerange(1)) & ...
          (clust.clustdata.timefilterranges(:,2) <= ...
          session(l).timerange(2)));
        % Remove the "all" filter from consideration (this can happen if the data
        % were acquired in a single recording session, in which case the "all"
        % interval is degenerate with a time filter)
        m(m == 1) = []; 
        % Sort the rows of timefilterranges that fall within the interval defined by
        % session(l).timerange
        t = sortrows(clust.clustdata.timefilterranges(m,:));
        % Check for agreement with the session time interval
        if ~all([t(1,1) t(end,end)] == session(l).timerange)
          error('matclust time filters do not agree with session time boundaries');
        end
        % If more than one timefilter falls within the interval, check that they
        % completely cover the interval and coincide at their boundaries
        if (size(t,1) > 1) && ~all(t(2:end,1) == t(1:end-1,2))
          error('matclust time filters do not agree with session time boundaries');
        end
        % Check the names of the matching time filters, and warn if they don't match
        % the name of the session
        if ~all(strcmp(session(l).epoch, ...
            regexp(clust.clustdata.timefilternames(m),'^\w+','match','once')))
          warning('clustdata.timefilternames do not match epoch name');
        end
        % Record these time filter indices for future lookup
        session(l).timefilters = m;
      end
      % Sanity check: make sure that no two sessions correspond to the same time
      % filters
      for l1 = 1:numel(session)
        for l2 = l1:numel(session)
          if (l1 ~= l2) && ~isempty(intersect( ...
              session(l1).timefilters,session(l2).timefilters))
            error('bug?! different sessions map onto the same time filters');
          end
        end
      end

      for l = 1:length(clust.clustattrib.clusters)
        if isempty(clust.clustattrib.clusters{l})
          continue;
        else
          % Iterate through defineaxes(:,3) to find time filters for this cluster
          for m = 1:size(clust.clustattrib.clusters{l}.defineaxes,1)
            t = clust.clustattrib.filterindex{ ...
                clust.clustattrib.clusters{l}.defineaxes(m,3)}(1);
            if ~ismember(t,unique(clust.clustdata.timefiltermemmap( ...
                clust.clustdata.timefiltermemmap ~= 0)))
              error('clustattrib.clusters{%d}.defineaxes is invalid',i);
            else
              % Keep record of which time filters this cluster is defined over (i.e.,
              % in which time filters is a cluster box is drawn?)
              if ~isfield(clust.clustattrib.clusters{l},'timefilters')
                clust.clustattrib.clusters{l}.timefilters = t;
              else
                clust.clustattrib.clusters{l}.timefilters = unique( ...
                    [clust.clustattrib.clusters{l}.timefilters; t]);
              end
            end
          end
        end
      end

      for l = 1:length(clust.clustattrib.clusters)
        if isempty(clust.clustattrib.clusters{l}) || isempty(clust.clustattrib.clusters{l}.index)
          continue;
        end
        for m = 1:numel(session)
          if ismember(1,clust.clustattrib.clusters{l}.timefilters)
            % cluster is defined for all times, including the current session
          else
            overlap = intersect(session(m).timefilters, ...
                clust.clustattrib.clusters{l}.timefilters);
            if ~isempty(overlap)
              if isempty(setdiff(session(m).timefilters,overlap))
                % cluster is defined for all time filters that correspond to this
                % session
                disp(sprintf('cluster %d in %s is defined for all of %s',l, ...
                    matclust_file,session(m).epoch));
              else
                % cluster is defined for some, but not all time filters that
                % correspond to this session (partial overlap)
                continue; % very important to skip to next iteration, because unit is not defined in this session!
              end
            else
              continue; % very important to skip to next iteration, because unit is not defined in this session!
            end
          end
          % Indices of the threshold-crossing events during this session that
          % belong % to this cluster
          event_idx = clust.clustattrib.clusters{l}.index( ...
              uint32(clust.clustdata.params(clust.clustattrib.clusters{l}.index,1) >= ...
              session(m).timerange(1)) & ...
              uint32(clust.clustdata.params(clust.clustattrib.clusters{l}.index,1) < ...
              session(m).timerange(2)));
          [junk, max_energy_channel] = ...
              max(median(sum(spike.samples(:,:,event_idx).^2,1),3));
          % complex spike index: proportion of first-order timestamp diffs that
          % are less than 15 milliseconds, such that the amplitude of the
          % second spike is smaller than the amplitude of the first spike in
          % the pair
          timestamp_diffs = diff(uint32(clust.clustdata.params(event_idx,1)));
          csi = nnz((timestamp_diffs < 150) & ...
              (clust.clustdata.params(event_idx(2:end),1+max_energy_channel) < ...
              clust.clustdata.params(event_idx(1:end-1),1+max_energy_channel))) / ...
              numel(timestamp_diffs);
          uid = sprintf('%s_day%d_tetrode%02d_cluster%d',session(m).subject,session(m).day, ...
              spike.tetrode,l);
          try
            waveforms = load(clust.clustdata.datafile);
          catch
            error('could not load waveforms from file %s', clust.clustdata.datafile);
          end
          waveforms.timestamp = waveforms.timestamp(event_idx);
          waveforms.samples = waveforms.samples(:,:,event_idx);
          shape = measure_spike_shape(waveforms,10,1e8);
          tally(end+1,1) = struct( ...
              'uid'                 , {uid}                                             , ...
              'subject'             , {session(m).subject}                              , ...
              'day'                 , {session(m).day}                                  , ...
              'tetrode'             , {spike.tetrode}                                   , ...
              'region'              , {spike.region}                                    , ...
              'epoch'               , {session(m).epoch}                                , ...
              'environment'               , {session(m).environment}                                , ...
              'timerange'           , {session(m).timerange}                            , ...
              'overall_mean_rate'   , {numel(event_idx)/double(diff(session(m).timerange))*1e4} , ...
              'complex_spike_index' , {csi}                                             , ...
              'timestamp'           , {timestamp}                                       , ...
              'samples'             , {samples}                                         , ...
              'shape'               , {shape}                                           );
        end
      end   
    end
  end

end

save('/home/smkim/scratch/tally.mat','tally');

