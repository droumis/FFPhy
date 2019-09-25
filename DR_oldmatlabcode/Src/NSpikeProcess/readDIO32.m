function [DIOout]=readDIO32(filename)

% change to that directory and saves the figure with file name
% animal_day_epoch

%cd(fileoutdir);
%filename = strcat(animalname,'_',dayt,'e',epocht,'.mat');
%save(filename,'Epochdata');


%cd(currentdir);
fid=fopen(filename);
DIO=fread(fid,'uint32');

i=1;
DIOout=[];

while i<size(DIO,1);
    curr(1,1)=DIO(i,1);
    curr(1,2)=DIO(i+1,1);
    curr(1,3)=DIO(1+2,1);
    DIOout=[DIOout;curr];
    
    i=i+3;
end
DIOtimes=DIOout(:,1);
DIOout=dec2bin(DIOout(:,2));

DIOout=fliplr(double(double(str2mat(DIOout)) - str2mat('0')));

DIOout=[DIOtimes,DIOout];

end




    
    






