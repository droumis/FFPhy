function data = labview_loadspike3(fn,nch,Fs,t1,t2)

data =[];
ind1=t1*8*Fs; ind2=t2*8*Fs; 
%if ind1==0, ind1=1; end

fid=fopen(fn,'rb','b'); s=fseek(fid, ind1*8, -1); 
y=fread(fid,ind2-ind1,'float'); s=length(y);
nn=floor(s/nch);
y=y(1:nn*nch);
for ii=1:nch
    data(1:nn,ii) = y(ii:nch:end);
end
size(data)
