function Plot_spectrogram(animaldir,anim,tet,day,epoch,tstart,tend)

% get the current directory to go back to after the function is done
currentdir = pwd

% ---- File loading ----
% See if day number needs 0
dsz = '';
   if (day < 10)
      dsz = '0';
   end
   
% Specify data directory and load the file 
animalname = animaldir;
datadir = '/data14/jai/';

% Converts the day and epoch into text
dayt = num2str(day);
epocht = num2str(epoch);

tsz = '';
    if (tet < 10)
        tsz = '0';
    end
    
tett = num2str(tet);

% Loads the eeg data
eegfilename = strcat(datadir,animaldir,'_/','EEG/',anim,'eeg',dsz,dayt,'-',epocht,'-',tsz,tett,'.mat');
eval(['load ', eegfilename]);

 
params = {};
params.Fs = 1500;
params.tapers = [3 5];

% set options
params.fpass = [0 500];
win = [1 0.05];

%search for time

tst=tstart;
tnt=tend;

tstart=timetrans({tstart},10000,2);
tend=timetrans({tend},10000,2);
 
% assign a temporary variable for eeg
e1 = eeg{1,day}{1,epoch}{1,tet}.data;
timeoffset=eeg{1,day}{1,epoch}{1,tet}.starttime(1,1)*10000;
samplecount1=ceil((tstart-timeoffset)*params.Fs/10000);
samplecount2=floor((tend-timeoffset)*params.Fs/10000);
e1 = eeg{1,day}{1,epoch}{1,tet}.data(samplecount1:samplecount2,1);
e1=e1';
  
   
% compute full spectrum
[S,t,f]= mtspecgramc(e1,[win(1) win(2)],params);

figure;
plot_matrix(S,t,f); 
title(['Spectrogram ',animalname,': Day ',dayt,' tetrode ', tett,' epoch ',epocht, ' (', tst, ' - ',tnt,')' ]);
ylabel('Hz');
xlabel('s');

cd('/data26/lynne/eeg/');

figurename = strcat(anim,'_d',dayt,'_e',epocht,'_t',tett,'_',tst,'-',tnt);
saveas(gcf, figurename, 'pdf');

cd(currentdir);
