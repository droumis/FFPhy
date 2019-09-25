
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S61'}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)
  clear('timestamps');
  this_subject = subjects{s};
  try
    load(sprintf('%s/%s/%s_session.mat',path_prefix,this_subject,this_subject));
    assert(exist('session') == 1);
  catch
    error('could not load session struct array for %s',this_subject);
  end
  unique_days = unique([session(:).day]);
  for d = 1:numel(unique_days)
    this_day = unique_days(d);
    videosync_filename = sprintf('%s/%s/behavior/%s_day%d.videosync', ...
        path_prefix,this_subject,this_subject,this_day);
    cpudsptimecheck_filename = sprintf( ...
        '%s/%s/behavior/%s_day%d.cpudsptimecheck', ...
        path_prefix,this_subject,this_subject,this_day);
    cpupostimestamp_filename = sprintf( ...
        '%s/%s/behavior/%s_day%d.cpupostimestamp', ...
        path_prefix,this_subject,this_subject,this_day);
    disp(cpudsptimecheck_filename)
    try
      timestamps = reconcile_video_timestamps(cpudsptimecheck_filename, ...
          cpupostimestamp_filename,videosync_filename);
    catch
      error('could not execute reconcile_video_timestamps(%s,%s,%s)', ...
          cpudsptimecheck_filename,cpupostimestamp_filename, ...
          videosync_filename);
    end
    save(sprintf('%s/%s/behavior/%s_day%d_videotimestamps.mat', ...
        path_prefix,this_subject,this_subject,this_day),'timestamps');
    binary_filename = sprintf('%s/%s/behavior/%s_day%d.videotimestamps', ...
        path_prefix,this_subject,this_subject,this_day);
    fid = fopen(binary_filename,'w');
    if (fid == -1)
      error('could not open write file %s for writing',binary_filename);
    end
    header = { ...
        '%BEGINHEADER', ...
        '% File type:	Binary', ...
        '% Extraction type:	video frame timestamps from DSP clock', ... 
        '% Fields:	  timestamp (uint32_t)', ...
        '%%ENDHEADER' };
    count = fprintf(fid,'%s\n',header{:});
    if (count ~= numel(horzcat(header{:})) + numel(header))
      error('could not write header to file %s',binary_filename);
    end        
    count = fwrite(fid,timestamps,'uint32');
    if (count ~= numel(timestamps))
      error('could not write data to file %s',binary_filename);
    end
    if (fclose(fid) ~= 0)
      error('could not close output file %s',binary_filename);
    end
    clear('timestamps');
  end
end
clear;
