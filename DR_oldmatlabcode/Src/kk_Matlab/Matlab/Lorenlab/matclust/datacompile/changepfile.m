function changepfile(filename,offset)

fid = fopen(filename,'r+');
junk = fread(fid,200,'char');

headersize = strfind(junk','ENDHEADER')+9;

 frewind(fid);
 junk = fread(fid,headersize,'char');
 
 times = fread(fid,inf,'uint32=>uint32',8);

 
 times = times+(uint32(offset)*10000);
 
 frewind(fid);
 junk = fread(fid,headersize,'char');
 
 fwrite(fid,times(1),'uint32');
 fwrite(fid,times(2:end),'uint32',8);
 

 %pos = fread(fid,[4,inf],'4*int16=>int16',4); 
 
% for i = 1:100 
%  times(i) = fread(fid,1,'uint32');
%  pos(1:4,i) = fread(fid,4,'int16');
% end
fclose(fid);
