function [params] = sj_HPexpt_baselinespecgram_gen(prefix, days, epochs, tets, varargin)
% Shantanu - Nov 2012. Generic version of get baseline spectrgroam and save
% This is for getting baseline values for given eeg tets for normalization when needed. Save These. 

% sj_HPexpt_baselinespecgram('HPa', 1, [1:5], [1:7], 0);



if nargin<1,
    keyboard
    error('Please enter Expt Prefix and Day No!');
end
if nargin<2,
    keyboard
    error('Please enter Day No!');
end
if nargin<3,
    epoch=1; %% Epochs 
end
if nargin<4,
    tets=1; %
end


% Define Chronux params
% -------------------------------------------
movingwin = [1000 100]/1000; %movingwin = [100 10]/1000;                
params.Fs = 1500;
params.fpass = [0 40]; % params.fpass = [0 400];
params.tapers = [3 5];
params.err = [2 0.05];

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'movingwin'
            movingwin = varargin{option+1};
        case 'fpass'
    	    params.fpass = varargin{option+1};
    end
end

savetag = '';


% SET DATA
% -------------------------------------------

switch prefix
    case 'HPa'
        rawdir = '/data25/sjadhav/HPExpt/HPa/';
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct/';
end
dir2=directoryname;


% Get EEGs and spectrogram it
% ------------------------------
cd([directoryname,'/EEG/']);
savedir = [directoryname,'/EEGSpec/'];

for d=1:length(days)
    
    day = days(d);
    if (day<10)
        daystring = ['0',num2str(day)];
    else
        daystring = num2str(day);
    end
    
    
    for t=1:length(tets)
        tet=tets(t);
        
        if (tet<10)
            tetstring = ['0',num2str(tet)];
        else
            tetstring = num2str(tet);
        end
        eegspec = []; dummy=[]; % Save 1 file for each tet in a day 
              
        for ep=1:length(epochs)
            eeg=[]; eeggnd=[];
            epoch=epochs(ep);
            disp(['Doing Day ',num2str(day) ', Ep ',num2str(epoch),', Tet ',num2str(tet)]);
            curreegfile = [dir2,'/EEG/',prefix,'eeg',daystring,'-',num2str(epoch),'-',tetstring];
            load(curreegfile);
            %eegspec{day}{epoch}{tet}.specgram = mtspecgramc(eeg{day}{epoch}{tet}.data,movingwin,params);
            specgram = mtspecgramc(eeg{day}{epoch}{tet}.data,movingwin,params);
            eegspec{day}{epoch}{tet}.meanspec = mean(specgram,1);
            eegspec{day}{epoch}{tet}.stdspec = std(specgram,1);
            dummy=[dummy; specgram]; % For combining across epochs            
            
            % Dont save specgram in file - its too big
            %eegspec{day}{epoch}{tet}.specgram  = [];
            
        end % end epochs
        % Also save mean and std for whole day. Save in fields for 1st epoch - [FIXED This has become last epoch now - by accident. Fix later]
        eegspec{day}{1}{tet}.meandayspec=mean(dummy,1);
        eegspec{day}{1}{tet}.stddayspec=std(dummy,1);
        % Save File for current day and tet
        savefile = [savedir,prefix,'eegspec',savetag,daystring,'-Tet',tetstring]
        save(savefile,'eegspec');
       
    end % end tets
    
    
end % end days


















i=1;
%set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.03 Screen(3)*0.2 Screen(4)*0.4])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



