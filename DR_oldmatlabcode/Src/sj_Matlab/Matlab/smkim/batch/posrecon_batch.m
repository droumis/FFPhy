
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
%subjects = {'S48','S58','S59','S60','S61'};
subjects = {'S61'};
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
    indexfilename = sprintf('%s/%s/behavior/%s_day%d.mpegoffset', ...
        path_prefix,this_subject,this_subject,day);
    timestampsfilename = sprintf('%s/%s/behavior/%s_day%d.videotimestamps', ...
        path_prefix,this_subject,this_subject,day);

    clear('cdata','rawpos');
    obj = posrecon.controller( ...
        mpegfilename,indexfilename,timestampsfilename,session(i).timerange);
    disp(session(i));
    % wait for obj to be deleted. the user is expected to export rawpos and
    % cdata
    while isvalid(obj)
      uiwait;
    end

    % selectively inherit fieldnames
    inherited_fields = {'subject','day','epoch','environment'};
    for f = 1:numel(inherited_fields)
      rawpos.(inherited_fields{f}) = session(i).(inherited_fields{f});
    end
    rawpos.units = 'pixels';
    rawpos.video_frame = cdata;
    rawpos.sources = { ...
        mpegfilename, indexfilename, timestampsfilename, session(i).source};

    if ~is_rawpos(rawpos)
      error('something went wrong');
    end
    
    try
      save(sprintf('%s/%s/behavior/%s_day%d_%s_rawpos.mat', ...
          path_prefix,this_subject,this_subject,day,epoch),'rawpos');
    catch
      error('could not save');
    end

  end
end
clear;

