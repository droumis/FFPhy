
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
%subjects = {'S48','S55','S58','S59','S60','S61'};
subjects = {'S48'};
dio_hostname = 'drizzle';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subjects)
  clear('dio');
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
    dio_filename = sprintf('%s/%s/behavior/%s_day%d.dio',path_prefix, ...
        this_subject,this_subject,this_day);
    if (exist(dio_filename) ~= 2)
      error('dio file %s does not exist',dio_filename);
    end
    config_filename = sprintf('%s/%s/%s_day%d.%s.config',path_prefix, ...
        this_subject,this_subject,this_day,dio_hostname);
    if (exist(config_filename) ~= 2)
      error('config file %s does not exist',config_filename);
    end
    try
      tmp_dio = read_dio(dio_filename,config_filename, ...
          session(find([session(:).day] == this_day)));
    catch
      error('read_dio failed: dio filename %s, config filename %s', ...
          dio_filename,config_filename);
    end
    if (exist('dio') ~= 1)
      dio = tmp_dio;
    else
      dio = [dio; tmp_dio];
    end
  end
  % save combined dio struct array across all days for each subject
  try
    save(sprintf('%s/%s/behavior/%s_dio.mat',path_prefix,this_subject, ...
      this_subject),'dio');
  catch
    error('could not save dio data for subject %s',this_subject);
  end
end
clear;
