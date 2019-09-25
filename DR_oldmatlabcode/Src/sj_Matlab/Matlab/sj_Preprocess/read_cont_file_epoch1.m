function [eegstruct,timestamps, samples, nsamps, samplingrate] = read_cont_file_epoch1(filename, starttime, endtime)
% Returns eegstruct for given timeperiod
%%
ff = fopen(filename);

instring = '';
while ~strncmp(instring,'%%ENDHEADER',11)
  instring=fgets(ff);
  if strncmp(instring, '% nsamplesperbuf',16)
    nsamps =  sscanf(instring,'%% nsamplesperbuf %d');
    %fprintf('N samples per buf is: %d\n', nsamps);
  end
  if strncmp(instring, '% samplingrate', 14)
    samplingrate = sscanf(instring,'%% samplingrate %d');
    %fprintf('Samplingrate is: %d\n', samplingrate);
  end
end

datastart_position = ftell(ff); %Important - this is where binary data starts


if nsamps ==1
    timestamps = fread(ff,inf,'uint32',2*nsamps); % read in timestamps, skipping data (2 bytes)
    fseek(ff,datastart_position+4, 'bof');
    samples = fread(ff,inf,'int16',4*nsamps); % read in data, skipping timestamps (4 bytes)
    
else
    
    fprintf('Its Complicated!\n');
    keyboard;
    error('nsampplesperbuffer is not equal to 1; Its complicated!');
    
    
%     %LONG WAY
%     timestamps = fread(ff,inf,'uint32',2*nsamps); % read in timestamps, skipping data
%     NTIMEPTS = length(timestamps);
%     NTOTALBYTES = NTIMEPTS*4 + NTIMEPTS*nsamps*2; % (Timepts*4 bytes + Nsamps perbuffer*Nbuffers*2 bytes)
%     fseek(ff, datastart_position+4, 'bof'); % Go to first data point skipping the first timestamp
%     cnt=1;
%     position = ftell(fid);
%     while position<=NTOTALBYTES-2*nsamps
%         samples(cnt)=fread(ff,2*nsamps,'int16'); % read in nsamps*2bytes chunk of data until next timestamp - NO Skipping
%         cnt=cnt+1;
%         position = ftell(fid);
%     end
%     
%     %SHORT WAY
% %     data = fread(ff,inf,'short'); %Read in all data as short, which is 16 bits/2 bytes (timestamps is 2 shorts, and smaples is 1 short) 
% %     data = reshape(data,2,length(data)/2); % Will be a multiple of 2
%     
% %  TO COMPLETE    
    
end   
fclose(ff);

% Now Cut Data to only include epoch
[outt] = timetrans({starttime endtime}, 10000, 2); % Convert to 1/10000 s, the units of timestamps
startt=outt(1); endt=outt(2);
start_ind = lookup(startt, timestamps); % Find timestamps idx closest to the starttime and endtime
end_ind = lookup(endt, timestamps);

% Make eegstruct to return
%text_time = cell2mat(timetrans(starttime,10000,1));

text_time{1}=starttime; text_time{2}=endtime;
eegstruct.descript = ['eeg data in file ' filename ':' text_time{1} 'to' text_time{2}];
eegstruct.fields = ['eegamplitude']; 
eegstruct.starttime = timestamps(start_ind)./10000; % Timestamps in 0.1ms units
eegstruct.samprate = samplingrate; 
eegstruct.nsampsperbuffer = nsamps; 
eegstruct.data = samples(start_ind:end_ind); % This only works if lth(samples)==length(timestamps)!! 




%% TO get epochs based on gaps in cont file
% dt = nsamps * 10000 / samplingrate;
% 
% max_dt = ceil(dt);
% min_dt = floor(dt);
% 
% epochs(1,1) = timestamps(1);
% n_epochs = 1;
% k = 2;
% while k < length(timestamps);
%   if (timestamps(k) - timestamps(k-1)) < min_dt
%     keyboard
%     error('Timestamp difference too small!');
%   elseif (timestamps(k) - timestamps(k-1) > max_dt)
%     epochs(n_epochs,2) = timestamps(k-1);
%     n_epochs = n_epochs + 1;
%     epochs(n_epochs,1) = timestamps(k);
%   end
%   k = k + 1;
% end
% epochs(n_epochs,2) = timestamps(end);
% fclose(ff);
