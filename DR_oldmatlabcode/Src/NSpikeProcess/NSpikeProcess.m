function NSpikeProcess(filedir,filename,filepref,day)
% Does NSpike preprocessing for the specified directory
% Does full NSpike extract
% filedir - full path of directory where the NSpike file is located
% filename - name of NSpike file to be extracted
% filepref - animal name, same as folder
% day - day to be processed

currentdir = pwd;

cd(filedir);
file=strcat(filedir,filename);
dsz = '';
   if (day < 10)
      dsz = '0';
   end
dayt = num2str(day);

directparam = regexp(filedir(2:end-1),'\/','split');

extractcont=['/opt/hippo/bin/nspike_extract',' -cont ', file];
% extracts dio data as uint32, therefore -dio instead of -diotext is called
extractdio=['/opt/hippo/bin/nspike_extract',' -dio ', file];
extractpos=['/opt/hippo/bin/nspike_extract',' -pos64 ', file];
extractspike=['/opt/hippo/bin/nspike_extract',' -spike ', file];
%timestamp=['/opt/hippo/bin/nspike_postimestamp',' -cd ', filepref,dsz,dayt,'.cpudsptimecheck –cp ',filepref,dsz,dayt,'.cpupostimestamp -ps possynctimes -o ',filepref,dsz,dayt,'.postimestamp'];
timestamp=sprintf('/opt/hippo/bin/nspike_postimestamp -cd %s%s%s.cpudsptimecheck -cp %s%s%s.cpupostimestamp -ps possynctimes -o %s%s%s.postimestamp',...
    filepref,dsz,dayt, filepref,dsz,dayt, filepref,dsz,dayt);



[a,b]=system(extractpos,'-echo');


%[a,b]=system(extractpos,'-echo');

[a,b]=system(extractcont,'-echo');

[a,b]=system(extractdio,'-echo');

[a,b]=system(extractspike,'-echo');

[a,b]=system(timestamp,'-echo');




% Change to animal directory
% 
cd(fullfile('/',directparam{1,1:3}));


animdirect=fullfile('/',directparam{1,1:3});

makedayparms([filepref,dayt]);

cd(filedir);

generateTimesFileFromCont;
% 
% system('chmod 777 *','-echo');
% 
daydirect=filedir;

animdirect=[animdirect,'_'];

if (exist(animdirect))~=7;
            %a folder needs to be created
            !mkdir animdirect 
        end
% 
% % modified from sj_diodayprocess since DIO data is read from uint32 instead
% % of text version
jy_diodayprocess(daydirect,animdirect,filepref,day);
% 
% system('chmod 777 *','-echo');
% 
% end
