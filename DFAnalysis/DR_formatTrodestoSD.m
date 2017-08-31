
dataz = out.channelData;

%[Time1Channel1, Time1Channel2, .. Time1ChannelN, Time2Channel1, ...]
combtetdata = reshape(dataz', [numel(dataz), 1]);

%write the file flattened to a signed 2byte (int16)
fid = fopen('D09_d06_e02_t19.dat','w'); %open file for writing
fwrite(fid, combtetdata, 'int16');
fclose(fid)

%% check the file.. looks good!
fid = fopen('D09_d06_e02_t19.dat','r');
datacheck = fread(fid,Inf,'int16');
fclose(fid)