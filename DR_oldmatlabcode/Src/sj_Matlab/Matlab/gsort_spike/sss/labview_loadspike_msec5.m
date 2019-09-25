function data = labview_loadspike4(fn,nch,Fs,t1,t2)

data =[];
ind1=floor(t1*Fs); %% Start Sample No to read for all channels / Scan No
ind2=floor(t2*Fs); %% End Sample No to Read for all channels / Scan No
%if ind1==0, ind1=1; end

ind1=ind1*nch;
ind2=ind2*nch;

fid=fopen(fn,'rb','b');

fseek_status = fseek(fid, ind1*4, -1)
% if t1==0,
%     fseek_status = fseek(fid, ind1*4, -1) 
% else
%     fseek_status = fseek(fid, ind1*4+4, -1) 
%     x=1
% end
data=fread(fid,[nch,floor((ind2-ind1)/nch)],'float'); 

%data=fscanf(fid,'%g',[8 Inf]);

%s=length(y);
%nn=floor(s/nch);
%y=y(1:nn*nch);
%for ii=1:nch
%    data(1:nn,ii) = y(ii:nch:end);
%end
size(data)
