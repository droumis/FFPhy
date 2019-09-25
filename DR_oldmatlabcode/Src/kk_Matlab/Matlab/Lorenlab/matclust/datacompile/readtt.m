function [timestamps, waves] = readtt(filename)  
%[timestamps, waves] = readtt(filename)
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
waveforms = fread(fid,[160,inf],'160*int16=>int16',4); 


fclose(fid);


for i = 1:4
    waves(:,:,i) = waveforms(i:4:end,:)';
end
waves = shiftdim(waves,1);


clear waveforms;



