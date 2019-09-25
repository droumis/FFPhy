function f = setdiofilter(f, filterString)
% function f = setdiofilter(f, filterString)
% For each epoch in the filter F, this function finds the time
% indices to triggered pulsetimes.

% pin = 15;
pin = 48; % default
TIMESTAMPRATE = 10000;

if isfield(f,'times')
  f = rmfield(f,'times');
  f = rmfield(f,'dioinfo');
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
  dio = loaddatastruct(datadir,animalprefix,'DIO');

  if isfield(f(an),'diopin') && ~isempty(f(an).diopin)
    pin = f(an).diopin;
  else
    f(an).diopin = pin;
    fprintf('setdiofilter: Setting default DIO pin 48\n');
  end

  if ~iscell(filterString)
    filterString = {filterString};
  end

  for i = 1:length(f(an).epochs)
    if isempty(f(an).epochs{i})
      f(an).times{i} = [];
      f(an).dioinfo{i} = [];
      continue;
    end
    for j = 1:size(f(an).epochs{i},1)
      pulses = dio{f(an).epochs{i}(j,1)}{f(an).epochs{i}(j,2)}{pin};
      if isempty(pulses.pulsetimes)
        f(an).times{i}{j} = [];
        f(an).dioinfo{i}{j} = [];
      end
      for k = 1:length(filterString) % allowing for the possibility of different DIO filters in one main filter
        if ~isempty(filterString{k})
          filtresult = evaluatefilter2({pulses}, filterString{k}, 'activefield','pulseind');
          idx = find(filtresult{1}(:,2));
          f(an).times{i}{j}{k} = pulses.pulsetimes(idx) / TIMESTAMPRATE;
          f(an).dioinfo{i}{j}{k} = [idx ...
            pulses.timesincelast(idx,1) ...
            pulses.pulselength(idx,1) ...
            pulses.pulseind(idx,1) ...
            pulses.frequency(idx,1) ...
            pulses.index_in_sequence(idx,1) ...
            pulses.frequency2(idx,1)];
        else
          f(an).times{i}{j}{k} = pulses.pulsetimes(:,1) / TIMESTAMPRATE;
          f(an).dioinfo{i}{j}{k} = [[1:length(pulses.pulsetimes)]' ...
            pulses.timesincelast(:,1) ...
            pulses.pulselength(:,1) ...
            pulses.pulseind(:,1) ...
            pulses.frequency(:,1) ...
            pulses.index_in_sequence(:,1) ...
            pulses.frequency2(:,1)];
        end
      end
      if (length(filterString) == 1)
        f(an).times{i}{j} = f(an).times{i}{j}{1};
        f(an).dioinfo{i}{j} = f(an).dioinfo{i}{j}{1};
      end
    end
    if ((size(f(an).epochs{i},1) == 1) & (length(filterString) == 1))
      f(an).times{i} = f(an).times{i}{1};
        f(an).dioinfo{i} = f(an).dioinfo{i}{1};
    end
  end
end
