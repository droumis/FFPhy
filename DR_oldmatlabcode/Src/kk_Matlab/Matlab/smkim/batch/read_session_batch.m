
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)
  clear('session');
  this_subject = subjects{s};

  txt_filename = sprintf('%s/%s/%s_sessions.txt',path_prefix,this_subject,...
      this_subject);
  if (exist(txt_filename) ~= 2)
    error('file %s does not exist',txt_filename);
  end
  try
    session = read_session(txt_filename);
  catch
    error('read_session("%s") failed',txt_filename);
  end
  
  for i = 1:length(session)
    % try to read the videosync file to check whether session start/end times
    % are correct
    flist = dir(sprintf('%s/%s/behavior/%s_day%d.videosync',path_prefix, ...
        this_subject,this_subject,session(i).day));
    if (numel(flist) ~= 1)
      error('i made a mistake');
    else
      videosync_filename = sprintf('%s/%s/behavior/%s',path_prefix, ...
          this_subject,flist(1).name);
      try
        sync_signal = read_continuous_mex(videosync_filename, ...
            session(i).timerange(1),session(i).timerange(2));
        assert(length(sync_signal) == 1);
        clear('sync_signal');
      catch
        disp(session(i));
        error('session includes a gap in data acquisition');
      end
    end
  end

  struct_filename = sprintf('%s/%s/%s_session.mat',path_prefix, ...
      this_subject,this_subject);
  try
    save(struct_filename,'session');
  catch
    error('could not save session as %s',struct_filename);
  end

end
clear;

