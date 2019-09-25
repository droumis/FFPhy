function rawpos = run(mpegfilename,indexfilename,timestampsfilename,timerange)
% function rawpos = run(mpegfilename,indexfilename,timestampsfilename,timerange)
%RUN Wrapper function for using posrecon package.
%
%   RAWPOS = POSRECON.RUN(MPEGFILENAME,INDEXFILENAME,TIMESTAMPSFILENAME,TIMERANGE)
%   returns a RAWPOS struct. TIMERANGE must be a 2-element row vector of uint32
%   timestamps in increasing order.
%
%Depends on:
%   IS_SESSION (written by smk)
%   TS2STR (written by smk)
%   IS_RAWPOS (written by smk)
%   VIDEO_READER class (written by smk)
%   KMEANS (MATLAB Statistics Toolbox)
%
%Written by SMK, 2009 October 20.
%

  if (exist('ts2str') ~= 2)
    error('posrecon package depends on m-file TS2STR (written by smk)');
  end
  if (exist('video_reader') ~= 2)
    error('posrecon package depends on VIDEO_READER class (written by smk)');
  end
  if (exist('kmeans') ~= 2)
    error(['posrecon package depends on m-file KMEANS ' ...
        '(MATLAB Statistics Toolbox)']);
  end

  if ~isa(timerange,'uint32') || ~isvector(timerange) || ...
      (size(timerange,1) ~= 1) || (size(timerange,2) ~= 2)
    error(['timerange must be a 2-element row vector of uint32 timestamps ' ...
        '[start_timestamp end_timestamp]']);
  end
    
  obj = posrecon.controller( ...
      mpegfilename,indexfilename,timestampsfilename,timerange);
  % wait for obj to be deleted. the user is expected to export rawpos and
  % cdata
  while isvalid(obj)
    uiwait;
  end

  if (evalin('base','exist(''cdata'')') == 1)
    evalin('base','rawpos.video_frame = cdata;');
  end

end

