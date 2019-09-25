
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set these values by hand
path_prefix = '/data14/smkim';
subjects = {'S48','S58','S59','S60','S61'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read_rawpos(posfilename,mpegfilename,mpegoffsetfilename,timestampfilename,session)

for s = 1:numel(subjects)
  this_subject = subjects{s};
  try
    load(sprintf('%s/%s/%s_session.mat',path_prefix, ...
        this_subject,this_subject));
    assert(exist('session') == 1);
  catch
    error('could not load session info for %s',this_subject);
  end

  for i = 1:numel(session)
    day = session(i).day;
    epoch = session(i).epoch;
    timerange = session(i).timerange;
    posfilename = sprintf('%s/%s/behavior/%s_day%d_%s.pos', ...
        path_prefix,this_subject,this_subject,day,epoch);
    mpegfilename = sprintf('%s/%s/behavior/%s_day%d.mpeg', ...
        path_prefix,this_subject,this_subject,day);
    mpegoffsetfilename = sprintf('%s/%s/behavior/%s_day%d.mpegoffset', ...
        path_prefix,this_subject,this_subject,day);
    timestampfilename = sprintf('%s/%s/behavior/%s_day%d.videotimestamps', ...
        path_prefix,this_subject,this_subject,day);
    if (exist(posfilename) == 2)
      rawpos = read_rawpos(posfilename,mpegfilename, ...
          mpegoffsetfilename,timestampfilename,session(i));
    else
      warning('expected pos file %s was not found',posfilename);
      continue; % skip to next iteration of the for loop
    end
    if (timerange(1) <= rawpos.timestamp(1))
      warning('%s first timestamp %d >= session start %d', ...
          posfilename,rawpos.timestamp(1),timerange(1));
    end
    if (timerange(end) >= rawpos.timestamp(end))
      warning('%s last timestamp %d <= session end %d = %d', ...
          posfilename,rawpos.timestamp(end),timerange(end));
    end
    save(sprintf('%s/%s/behavior/%s_day%d_%s_rawpos.mat', ...
        path_prefix,this_subject,this_subject,day,epoch),'rawpos');
    clear('rawpos');
  end

end
%clear;

%{
% append one more sample to the pos struct (duplicate the last one)
posrecstruct = struct( ...
    'name',{'timestamp','xfront','yfront','xback','yback'}, ...
    'type',{'uint32','int16','int16','int16','int16'}, ...
    'count',{1,1,1,1,1});
pos = read_binary_records(posfilename,'%%ENDHEADER\n',posrecstruct,Inf);
pos.timestamp = [pos.timestamp; timestamps.timestamp(idx+1)];
pos.xfront = [pos.xfront; pos.xfront(end)];
pos.yfront = [pos.yfront; pos.yfront(end)];
pos.xback = [pos.xback; pos.xback(end)];
pos.yback = [pos.yback; pos.yback(end)];
% append this extra sample to end of file
fid = fopen(posfilename,'a');
fwrite(fid,pos.timestamp(end),'uint32');
fwrite(fid,int16(pos.xfront(end)),'int16');
fwrite(fid,int16(pos.yfront(end)),'int16');
fwrite(fid,int16(pos.xback(end)),'int16');
fwrite(fid,int16(pos.yback(end)),'int16');
fclose(fid);
%}

