
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48','S58','S59','S60','S61'};
executable = '/home/smkim/bin/nspike_fixpos';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)
  this_subject = subjects{s};
  if (exist([path_prefix '/' this_subject]) ~= 7);  
    error('could not find folder %s/%s',path_prefix,this_subject);
  end 

  % load session info
  try
    load(sprintf('%s/%s/%s_session.mat',path_prefix, ...
        this_subject,this_subject));
    assert(exist('session') == 1);
  catch
    error('could not load session info for %s',this_subject);
  end
  
  % run nspike_fixpos for each session
  for i = 1:numel(session)
    if ~strcmp(this_subject,session(i).subject)
      error('subject of session info (%s) does not match expected subject %s', ...
          session(i).subject,this_subject);
    end
    day = session(i).day;
    epoch = session(i).epoch;
    timerange = session(i).timerange;

    if ~strncmp('run',epoch,3)
      % skip if it isn't a run
      continue;
    end

    mpegfilename = sprintf('%s/%s/behavior/%s_day%d.mpeg', ...
        path_prefix,this_subject,this_subject,day);
    offsetfilename = sprintf('%s/%s/behavior/%s_day%d.mpegoffset', ...
        path_prefix,this_subject,this_subject,day);
    timestampsfilename = sprintf('%s/%s/behavior/%s_day%d.videotimestamps', ...
        path_prefix,this_subject,this_subject,day);
    outfilename = sprintf('%s/%s/behavior/%s_day%d_%s.pos',path_prefix, ...
        this_subject,this_subject,day,epoch);
    % start and end with extra 1.5 second buffer at edges of the epoch
    tstart = ts2str(timerange(1) - 1.5e4);
    tstart = tstart(2:9);
    tend = ts2str(timerange(end) + 1.5e4);
    tend = tend(2:9);
    if (system(sprintf(['%s -p %s -f64 %s -t %s -o %s -thresh 200 ' ...
        '-skip 0 -tstart %s -tend %s -speed 30 -lights 2'], ...
        executable,mpegfilename,offsetfilename,timestampsfilename, ...
        outfilename,tstart,tend)) < 0)
      error(['could not run nspike_fixpos for %s day %d %s', ...
          'this_subject,day,epoch']);
    end

  end

end
clear;

