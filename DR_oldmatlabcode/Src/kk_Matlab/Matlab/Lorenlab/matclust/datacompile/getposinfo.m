function [times, pos] = getposinfo(filename)

fid = fopen(filename,'r');
junk = fread(fid,200,'char');

headersize = strfind(junk','ENDHEADER')+9;

 frewind(fid);
 junk = fread(fid,headersize,'char');
 
 times = fread(fid,inf,'uint32=>uint32',8);
 frewind(fid);
 junk = fread(fid,headersize,'char');
 junk = fread(fid,1,'uint32');
 pos = fread(fid,[4,inf],'4*int16=>int16',4); 
 
% for i = 1:100 
%  times(i) = fread(fid,1,'uint32');
%  pos(1:4,i) = fread(fid,4,'int16');
% end
 fclose(fid);
pos = pos';

if (length(times) == (size(pos,1)+1))
    %sometimes, the last timestamp is read without a position
    times = times(1:end-1);
end
%[times, pos] = loadposfile(filename,headersize);
