function tetinfo = generateTetInfoFromConfig(animalDir,varargin)

masterMachine = '';

otherArgs = procOptions(varargin);

configFilenames = dir(fullfile(animalDir,'*.dat.config'));

for i = 1:length(configFilenames)
  tok = regexp(configFilenames(i).name,'.*\.(.*)\.dat\.config','tokens');
  if ~isempty(tok)
    machines{i} = tok{1}{1};
  end
end

if length(configFilenames)==0
  fprintf('No config files found.\n');
  return;
end

umachines = unique(machines);
fprintf('Machine names detected:');
for m = 1:length(umachines)
  fprintf(' [%d]%s',m,umachines{m});
end
m = input('. Please enter machine number: ');

machine = umachines{m};

configFilename = configFilenames(find(strcmp(machines,machine),1)).name;
configFilename = fullfile(animalDir,configFilename);

confFile = fopen(configFilename,'rt');

% go through config file until we get to the channel defs
s = fgets(confFile);
t = regexp(s,'\s*channel\s*(\d+)','tokens');
while isempty(t)
  s = fgets(confFile);
  t = regexp(s,'\s*channel\s*(\d+)','tokens');
end
% ok, now, a channel at a time, read in important info
ch = str2num(t{1}{1});
s = fgets(confFile);
t = regexp(s,'^\s*(\w+)\s*(\d+)\s*(\d*)','tokens');
while ischar(s) & ~isempty(t)
  switch lower(t{1}{1})
  case 'number'
    tetdata(ch+1).number = str2num(t{1}{2});
  case 'depth'
    tetdata(ch+1).depth = str2num(t{1}{2});
  case 'electchan'
    tetdata(ch+1).chan = str2num(t{1}{2});
  case 'refelect'
    tetdata(ch+1).refelect = str2num(t{1}{2});
  case 'refchan'
    tetdata(ch+1).refchan = str2num(t{1}{2});
  case 'filter'
    tetdata(ch+1).filter = [str2num(t{1}{2}) str2num(t{1}{3})];
  case 'dspchan'
    tetdata(ch+1).dspchan = str2num(t{1}{2});
  case 'channel'
    ch = str2num(t{1}{2});
  end
  s = fgets(confFile);
  t = regexp(s,'^\s*(\w+)\s*(\d+)\s*(\d*)','tokens');
end


% Get rid of EEG video channel:
tetdata([tetdata.dspchan] == 63) = [];


% Figure out number of epochs (tetinfo is {day}{epoch})
timesdir = dir(fullfile(animalDir,'times.mat'));
if isempty(timesdir)
  nepochs = input(strcat('Could not find [times.eeg] file. How many epochs?'));
else
  load(fullfile(animalDir,'times.mat'));
  nepochs = size(ranges,1)-1;
end

for n = 1:nepochs
  for i = 1:length(tetdata)
    tetinfo{n}{tetdata(i).number}.tetrode = tetdata(i).number;
    tetinfo{n}{tetdata(i).number}.chan = tetdata(i).chan;
    tetinfo{n}{tetdata(i).number}.depth = tetdata(i).depth;
    tetinfo{n}{tetdata(i).number}.refelect = tetdata(i).refelect;
    tetinfo{n}{tetdata(i).number}.refchan = tetdata(i).refchan;
    tetinfo{n}{tetdata(i).number}.filters = tetdata(i).filter;
  end
end


