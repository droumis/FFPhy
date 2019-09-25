function batchDIOextract(animaldir,day,configfile,wells,odour)
% Extracts all epochs for each day
% animaldir - name of directory containing processed data
% day - day of experiment
% configfile - NSpike configfile for the day, stored in configdir currently
% /data14/jai/Config/
% well - vector indicating odour used eg [1 3]
% odour - odours used for the wells, specify in same order {'orange'
% 'melon'}
% ---currently configured for return home if error strategy--- 


% get the current directory to go back to after the function is done
currentdir = pwd;

% ---- File loading ----
% See if day number needs 0
dsz = '';
   if (day < 10)
      dsz = '0';
   end
   
% Specify data directory and load the file 
animalname = animaldir(1:end-1);
datadir = '/data14/jai/';

% Converts the day and epoch into text
dayt = num2str(day);

sfilename = strcat(datadir,animaldir,'/',animalname,'DIO',dsz,dayt,'.mat');   
 
% Load the mat file
eval(['load ', sfilename]); 

% Get number of epochs

epoch = size(DIO{1,day},2);

% Process all epochs

Data={};
for i=day;
    for j=1:epoch;
        % calls the ReturnHome version of DIO process
    
        Epochdata=DIOprocessReturnHome(animaldir,i,j,configfile,wells,odour);
    
    % calls find correct well version
    %Epochdata=DIOprocess(animaldir,i,j,configfile,wells,odour);
    Data{1,day}{1,j}=Epochdata;
    j=j+1;
    end
end


% change to that directory and saves the figure with file name
% animal_day_epoch

fileoutdir=strcat(datadir,animaldir,'/');

cd(fileoutdir);
filename = strcat(animalname,'_',dayt,'.mat');
save(filename,'Data');
