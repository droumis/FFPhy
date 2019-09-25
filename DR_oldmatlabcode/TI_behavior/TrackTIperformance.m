% dio key.. testdio.dio
%bits 16-23 tracklight 1 - home leds
%bits 24-30 trackrew 1 - home rew
%bits 0-7 trackrew1 - homerew

%missing something here...
%rew in are 46 - 54. home sensor is 54... yes
%rew out are 40 - 46. home rew out is 46.. yes
%leds are 32 - 38. home leds are 38, 39..
% ...missing 51, 47,

%if there's an output, not at the home rew, than it must have been a correct trial..
% so i can go through and count the number of rewards given out as the total num of correct.
% To get the total number of trials, i need to get the total number of non-home rew sensors that immediately follow a home output...

%home rew in. home rew out. A



% abimp = str2double( regexp( fscanf( fopen('T03_AB.dio.txt') ,'%1c'), '[0-1]', 'match'));
clear all
fid = fopen('/data19/droumis/TransInf/T03_AB/T03_AB.dio.txt','r'); %open the file and give it an int handle

%get num of rows in text

numrows = 0;
tline = fgetl(fid);
while ischar(tline)
    tline = fgetl(fid);
    numrows = numrows+1;
end
fclose(fid);
fid = fopen('/data19/droumis/TransInf/T03_AB/T03_AB.dio.txt','r'); %need to do this again to start at the first line

% diobits = nan(numrows,4); % prealocate diobits for different number of bits at each time point.. it's never more than 3 right?
for i = 1:numrows;
    
    line = fgetl(fid); %this will get the next line in the txt file..
    diostr = textscan(line, '%s'); %this will separate the time stamp from the dio bits
    timestruct = str2num(diostr{1}{1});%put the time stamp in the first col
    diovec = double(diostr{1}{2} == '1'); %separate the dio bits as a vector 
    diostruct(i,:) = [timestruct diovec]; %append to the time stamp
    bits = find(diovec == 1);
    diobits{i,1} = [timestruct bits]; %create structure with time stamp and bits on... 
    bitsstruct{i,1} = bits;
        
   %translate the diovec to bits corresponding
    
    
end