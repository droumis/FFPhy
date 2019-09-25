
windowLength = 0.3;
winOffset = 0.1;

t = [0:30000*windowLength-1]/30000 - winOffset;

d = dir('*.cont');
if isempty(d)
    error('No CONT files found.');
end

% generateTimesFromCont
% stimdio = createstimstruct; save stimdio1 stimdio

%%% NOTE- pin "32" appears to correspond to "A" = CA3 = pin 16
%%% NOTE- pin "31" appears to correspond to "B" = EC = pin 17

load times
ne = size(ranges,1)-1;
for i = 1:length(d)
  for e = 1:ne
    filename = d(i).name;
    load stimdio1;
    if (e <= length(stimdio)) && (~isempty(stimdio{e}))
      ts = stimdio{e}.pulsetimes(:,1);
      [eeg1{e}{i}, times, samplingRate] = cont_window_c(filename, ts/10000 - winOffset, windowLength);
    end
    load stimdio2;
    if (e <= length(stimdio)) && (~isempty(stimdio{e}))
      ts = stimdio{e}.pulsetimes(:,1);
      [eeg2{e}{i}, times, samplingRate] = cont_window_c(filename, ts/10000 - winOffset, windowLength);
    end
  end
end

for e = 1:ne
  if (e <= length(eeg1)) && ~isempty(eeg1{e})
    mEEG1{e} = cat(3,eeg1{e}{:});
    tmp = mEEG1{e}(:,:,1:2:end);
    mEEG1{e}(:,:,1:2:end) = mEEG1{e}(:,:,2:2:end);
    mEEG1{e}(:,:,2:2:end) = tmp;
  end

  if (e <= length(eeg2)) && ~isempty(eeg2{e})
    mEEG2{e} = cat(3,eeg2{e}{:});
    tmp = mEEG2{e}(:,:,1:2:end);
    mEEG2{e}(:,:,1:2:end) = mEEG2{e}(:,:,2:2:end);
    mEEG2{e}(:,:,2:2:end) = tmp;
  end
end


