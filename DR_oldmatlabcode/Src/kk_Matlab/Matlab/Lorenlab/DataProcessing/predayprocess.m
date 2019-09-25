function predayprocess(daydirect, animdirect, fileprefix, daynum, varargin)
% function predayprocess(daydirect, animdirect, fileprefix, daynum, varargin)
%
% processstimdio = 1;
% generateTask = 1;
% generateTetinfo = 1;

processstimdio = 1;
generateTask = 1;
generateTetinfo = 1;
varprefix = '';

[otherOptions] = procOptions(varargin);

if (exist(fullfile(pwd,daydirect),'dir') == 7) % is day direct here now?
  daydirect = fullfile(pwd,daydirect);
end

[daypath,dayname] = fileparts(daydirect);

% timesdir = dir(fullfile(animdirec,'times.mat'));
% if isempty(timesdir)
  % error('Could not find times.mat file showing epoch boundaries.');
% end
% times = load(fullfile(daydirect,'times.mat'));
% NEpochs = size(times.ranges,1) - 1;

if (processstimdio)
  stimdiofiles = dir(fullfile(daydirect,'*stimdio*.mat'));
  if length(stimdiofiles) > 1
    error('[dayprocess: processstimdio] Multiple stimdio mat files found.');
  end
  tmp_stimdio = load(fullfile(daydirect,stimdiofiles.name),'stimdio');
  stimdio{daynum} = tmp_stimdio.stimdio;

  eval([varprefix,'stimdio = stimdio;']);
  save(fullfile(animdirect,sprintf('%sstimdio%02d.mat',fileprefix,daynum)),[varprefix,'stimdio']);
end

if (generateTask)
  taskfiles = dir(fullfile(daydirect,'*task*.mat'));
  if length(taskfiles) == 1
    tmp_task = load(fullfile(daydirect,taskfiles.name),'task');
    task{daynum} = tmp_task.task;
    eval([varprefix,'task = task;']);
    save(fullfile(animdirect,sprintf('%stask%02d.mat',fileprefix,daynum)),[varprefix,'task']);
  end
end

if (generateTetinfo)
  tetfiles = dir(fullfile(daydirect,'*tetinfo*.mat'));
  if length(tetfiles) == 1
    tmp_tet = load(fullfile(daydirect,tetfiles.name),'tetinfo');
    % tetinfo is a single file for all days and epochs
    % check to see if it already exists
    fname = fullfile(animdirect,sprintf('%stetinfo.mat',fileprefix));
    load(fname);
    tetinfo{daynum} = tmp_tet.tetinfo;
    eval([varprefix,'tetinfo = tetinfo;']);
    save(fname,[varprefix,'tetinfo']);
  end
end

