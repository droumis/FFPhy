function [epochs, nsamps, samplingrate] = cont_file_info(filename)

ff = fopen(filename);

instring = '';
while ~strncmp(instring,'%%ENDHEADER',11)
  instring=fgets(ff);
  if strncmp(instring, '% nsamplesperbuf',16)
    nsamps =  sscanf(instring,'%% nsamplesperbuf %d');
    fprintf('N samples per buf is: %d\n', nsamps);
  end
  if strncmp(instring, '% samplingrate', 14)
    samplingrate = sscanf(instring,'%% samplingrate %d');
    fprintf('Samplingrate is: %d\n', samplingrate);
  end
end

timestamps = fread(ff,inf,'uint32',2*nsamps); % read in timestamps, skipping data

dt = nsamps * 10000 / samplingrate;

max_dt = ceil(dt);
min_dt = floor(dt);

epochs(1,1) = timestamps(1);
n_epochs = 1;
k = 2;
while k < length(timestamps);
  if (timestamps(k) - timestamps(k-1)) < min_dt
    keyboard
    error('Timestamp difference too small!');
  elseif (timestamps(k) - timestamps(k-1) > max_dt)
    epochs(n_epochs,2) = timestamps(k-1);
    n_epochs = n_epochs + 1;
    epochs(n_epochs,1) = timestamps(k);
  end
  k = k + 1;
end
epochs(n_epochs,2) = timestamps(end);

fclose(ff);
