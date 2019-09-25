function [timestamps, waves] = readtt_rethresh(filename)
%[timestamps, waves] = readtt_rethresh(filename)
%filename-- a string containing the name of the tt file
%timestamps-- a vector containing the timestamps of the spikes (in uint32 format)
%waves-- a 40 by 4 by n matrix containing 40 points for the four channels for each spike (in int16 format)

fid = fopen(filename,'r');
junk = fread(fid,200,'char');
headersize = strfind(junk','ENDHEADER')+9;
frewind(fid);
junk = fread(fid,headersize,'char');

timestamps = fread(fid,inf,'uint32=>uint32',320);
frewind(fid);
junk = fread(fid,headersize,'char');
junk = fread(fid,1,'uint32');

n=100000;
waves = [];
a = 0;
indextotal = [];
while 1 %~feof(fid) %go until get to end of file
    waveforms = fread(fid,[160,n],'160*int16=>int16',4); %read in 1 million bits

    m = max(waveforms);
    index = find(m > 50 & m < 3000); %find entries greater than 50 and less than 2500mV
    %index = unique(J); %select all columns or wwaveforms that are within threshold
    
    for i = 1:4  %restructure into 3dimensional matrix
        wavestemp(:,:,i) = waveforms(i:4:end,index)';
    end
    wavestemp2 = shiftdim(wavestemp,1);
    waves = cat(3,waves, wavestemp2); %add next section of file
    
    clear wavestemp;
    clear wavestemp2;
indextotal = [indextotal index+a];
 a = a+n; 
 if feof(fid)
     break
 end
end
%keyboard
fclose(fid); %close file
timestamps = timestamps(indextotal);
clear waveforms;



