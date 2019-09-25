function f = setstimfilter(f, filterString)
% function f = setdiofilter(f, filterString)
% For each epoch in the filter F, this function finds the time
% indices to triggered pulsetimes.

TIMESTAMPRATE = 10000;

if isfield(f,'times')
  f = rmfield(f,'times');
end
if isfield(f,'stiminfo')
  f = rmfield(f,'stiminfo');
end

for an = 1:length(f)
  if isempty(f(an).animal)
    error(['You must define an animal for the filter before filtering the tetrodes'])
  end
  if isempty(f(an).epochs)
    %error(['You must define the desired epochs for the filter before filtering the tetrodes'])
    continue;
  end
  datadir = f(an).animal{2};
  animalprefix = f(an).animal{3};
  stimdio = loaddatastruct(datadir,animalprefix,'stimdio');

  if ~iscell(filterString)
    filterString = {filterString};
  end

  for i = 1:length(f(an).epochs)
    f(an).times{i} = [];
    f(an).pulses{i} = [];
    if isempty(f(an).epochs{i})
      continue;
    end
    for j = 1:size(f(an).epochs{i},1)
      if length(stimdio) < f(an).epochs{i}(j,1) || ...
        isempty(stimdio{f(an).epochs{i}(j,1)}) || ...
        length(stimdio{f(an).epochs{i}(j,1)}) < f(an).epochs{i}(j,2) || ...
        isempty(stimdio{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}) 
        pulses = [];
      else
        pulses = stimdio{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)};
      end

      for k = 1:length(filterString) % allowing for the possibility of different DIO filters in one main filter
        if isempty(pulses) || isempty(pulses.pulsetimes)
          f(an).times{i}{j}{k} = [];
          f(an).pulses{i}{j}{k} = [];
          continue;
        end

        if ~isempty(filterString{k})
          filtresult = evaluatefilter2({pulses}, filterString{k}, 'activefield','pulseind');
          idx = find(filtresult{1}(:,2));
          f(an).times{i}{j}{k} = pulses.pulsetimes(idx) / TIMESTAMPRATE;
          f(an).pulses{i}{j}{k} = pulses.pulseStruct(idx);
        else
          f(an).times{i}{j}{k} = pulses.pulsetimes(:,1) / TIMESTAMPRATE;
          f(an).pulses{i}{j}{k} = pulses.pulseStruct;
        end
      end

      if (length(filterString) == 1)
        f(an).times{i}{j} = f(an).times{i}{j}{1};
        f(an).pulses{i}{j} = f(an).pulses{i}{j}{1};
      end
    end
    if ((size(f(an).epochs{i},1) == 1) & (length(filterString) == 1))
      f(an).times{i} = f(an).times{i}{1};
      f(an).pulses{i} = f(an).pulses{i}{1};
    end
  end
end

