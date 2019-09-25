function stimdio = createstimstruct(varargin)
% function stimdio = createstimstruct(varargin)
% daydir = pwd;
% pin = [];  --- ZERO based pins
% pindescriptions = {};  --- cell array of strings for labels of DIOs

daydir = pwd;
pin = [];
pindescriptions = {};

[otherArgs] = procOptions(varargin);

if exist(fullfile(daydir,'times.mat'))
  times = gettimes(fullfile(daydir,'times.mat'));
else
  times = [];
end

% Read in DIO file
files = dir(fullfile(daydir,'*.dio'));
if isempty(files)
   error('No dio file found in %s.',daydirect);
end
if length(files) > 1
   error('Multipled dio files found in %s.',daydirect);
end
fid = fopen(fullfile(daydir,files(1).name),'r');
dio = textscan(fid,'%d %s');
fclose(fid);
diotimes = dio{1};

for i = 1:length(dio{2})
  dlen(i) = length(dio{2}{i});
end
if std(dlen) ~= 0
  error('File is corrupted.');
end

diovals = double(str2mat(dio{2})) - str2mat('0');

% correct a bug in the DIO text output - reverse each 16 bit frame
diovals = diovals(:,[[16:-1:1] [32:-1:17] [48:-1:33] [64:-1:49]]);
% clear fid files dio;

if isempty(diovals)
  return;
end

if isempty(pin)
  biphasicStim = -1;
  while ~ismember(upper(biphasicStim),['','M','B'])
    biphasicStim = input('(B)iphasic or (M)onophasic stimulation?  ([M])','s');
  end

  if (isempty(biphasicStim) | upper(biphasicStim)=='M')
    biphasicStim = 0;
  else
    biphasicStim = 1;
  end

  changing_pins = find(sum(diovals)) - 1;
  pin = input(strcat('Use which pin as input [ ',sprintf('%d ',changing_pins),']?'));
end

pin = pin + 1; % Assume inputs are zero-based


for i = 1:length(pin)

  [diopulses, pulseStruct] = extractStimDataFromFile(pin(i), diotimes, diovals, otherArgs{:});

  if ~isempty(times)
    epochs = findEpoch(times,diopulses.pulsetimes(:,1));
  else
    epochs = ones(size(diopulses.pulsetimes(:,1)));
  end
  eps = unique(epochs);

  if (sum(eps==0) > 0)
    warning(sprintf('Stim extraction yields %d pulses not within defined epochs.',...
      sum(epochs==0)));
    eps(eps==0) = [];
  end
  for ep = eps(:)'
    idx = epochs==ep;
    stimdio{ep}{i}.pin = pin(i) - 1;
    if ~isempty(pindescriptions)
      stimdio{ep}{i}.type = pindescriptions{i};
    end
    stimdio{ep}{i}.pulsetimes = diopulses.pulsetimes(idx,:);
    stimdio{ep}{i}.timesincelast = diopulses.timesincelast(idx);
    stimdio{ep}{i}.timeuntilnext = diopulses.timeuntilnext(idx);
    stimdio{ep}{i}.pulselength = diopulses.pulselength(idx,:);
    stimdio{ep}{i}.pulseind = diopulses.pulseind(idx,:);
    stimdio{ep}{i}.frequency = diopulses.frequency(idx,:);
    stimdio{ep}{i}.frequency2 = diopulses.frequency2(idx,:);
    stimdio{ep}{i}.index_in_sequence = diopulses.index_in_sequence(idx,:);
    stimdio{ep}{i}.pulseStruct = pulseStruct(idx);
  end
end

if (nargout==0)
  save stimdio stimdio;
end

function epoch = findEpoch(times, timestamp)
epoch = zeros(length(timestamp),1);
for i = 1:length(times)
  if isempty(times(i).range)
    continue;
  end
  epoch((timestamp >= times(i).range(1)) & (timestamp <= times(i).range(2))) = i;
end

function [diopulses, pulseStruct] = extractStimDataFromFile(pin, diotimes, diovals, varargin)

TIMESTAMPRATE = 10000;
pot_threshold = 10000; % 1 seconds
[otherArgs] = procOptions(varargin);

diovector = diovals(:,pin);
ind = find(diovector,1);

if isempty(ind)
  error('Pin not found in DIO');
end

if ~isempty(ind)
  % Note - I've noticed both rising and falling edges missed;
  % For simplicity I assume that only falling edges are ever
  % neglected.
  i = 1;
  pulse(i).tstart = diotimes(ind);
  done = 0;
  while ~done
     found = 0;
     k = ind + 1;
     while (found == 0)
        if (k > length(diovector))
           found = 1;
           pulse(i).tend = Inf;
           pulse(i).length = NaN;
           done = 1;
        elseif (diovector(k) == 0)
           found = 1;
           pulse(i).tend = diotimes(k);
           pulse(i).length = double(pulse(i).tend - pulse(i).tstart);
           ind = k + 1;
           % find next pulse start
           if (ind <= length(diovector))
              while ((diovector(ind) == 0) & (ind < length(diovector)))
                 ind = ind + 1;
              end
              if (diovector(ind) == 1)
                 i = i + 1;
                 pulse(i).tstart = diotimes(ind);
              else
                 done = 1;
              end
           else
              done = 1;
           end
        elseif (sum(diovals(k,:) ~= diovals(ind,:)) > 0)
           % Did a different pin change, causing DIO to
           % record?
           k = k + 1;
        else  % This accounts for missed falling edge
           % and only occurs if the same DIO string has
           % occured twice with different time values.
           found = 1;
           pulse(i).tend = diotimes(k);
           pulse(i).length = NaN;
           i = i + 1;
           pulse(i).tstart = diotimes(k);
           ind = k;
        end
     end
  end
end
for i = 1:length(pulse)
  pulse(i).ind = i;
end
diopulses.pulsetimes = double([cat(1,pulse.tstart) cat(1,pulse.tend)]);
diopulses.timesincelast = [inf;diff(diopulses.pulsetimes(:,1))];
if length(diopulses.timesincelast) > 1
  diopulses.timeuntilnext = [diopulses.timesincelast(2:end); inf];
else
  diopulses.timeuntilnext = inf;
end
diopulses.pulselength = cat(1,pulse.length);
diopulses.pulseind = cat(1,pulse.ind);


frequency = 10000./diopulses.timesincelast;
frequency2 = frequency;
if length(frequency) > 1
   frequency = [frequency(2:end); nan];
end
NP = length(diopulses.timesincelast);
index_in_sequence = -ones(NP,1);
first_pulse = -1;
for k = 1:NP
   if (diopulses.timesincelast(k) < pot_threshold)
      if (first_pulse < 0)
         first_pulse = k;
      else
         frequency2(k) = frequency(first_pulse-1);
      end
      index_in_sequence(k) = index_in_sequence(k-1) + 1;
   else
      first_pulse = -1;
   end
end
pot_inds = find(index_in_sequence >= 0);
index_in_sequence(pot_inds) = index_in_sequence(pot_inds) + 1;
pot_starts = find(index_in_sequence == 1);
index_in_sequence(pot_starts - 1) = 0;

diopulses.index_in_sequence = index_in_sequence;
diopulses.frequency2 = frequency2;
diopulses.frequency = frequency;


for i = 1:length(pulse)
  pulseStruct(i).pulsetimes = diopulses.pulsetimes(i,:);
  pulseStruct(i).timesincelast = diopulses.timesincelast(i);
  pulseStruct(i).timeuntilnext = diopulses.timeuntilnext(i);
  pulseStruct(i).pulselength = diopulses.pulselength(i,:);
  pulseStruct(i).pulseind = diopulses.pulseind(i,:);
  pulseStruct(i).frequency = diopulses.frequency(i,:);
  pulseStruct(i).frequency2 = diopulses.frequency2(i,:);
  pulseStruct(i).index_in_sequence = diopulses.index_in_sequence(i,:);
end
