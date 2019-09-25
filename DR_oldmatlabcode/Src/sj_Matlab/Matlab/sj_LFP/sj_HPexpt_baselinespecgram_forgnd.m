function [params] = sj_HPexpt_baselinespecgram_forgnd(prefix, days, epochs, tets, varargin)
% Shantanu - Nov 2012. 
% To avoid errors, make a separate file for wrt gnd spectrogram and ref spectrogram
% Analogous code in sj_HPexpt_baselinespecgram.m and sj_HPexpt_baselinespecgram_forref.m respectively

% From Kenny and Maggies - event_spectrograms .m and calcriptriggerredspectrograms.m respectively
% This is for getting baseline values for given eeg tets for normalization when needed. Save These. 

% sj_HPexpt_baselinespecgram_gnd('HPa', 8, [1:5], [1:16]);
% sj_HPexpt_baselinespecgram_gnd('HPa', 2, [1:5], [1:16]);

% sj_HPexpt_baselinespecgram_gnd('HPb', 1, [1:7], 1:16, 1)
% sj_HPexpt_baselinespecgram_gnd('HPa', 2, [1:5], [14,15], 1);
% sj_HPexpt_baselinespecgram_gnd('HPb', 2, [1:5], [4,9,16], 1);
% sj_HPexpt_baselinespecgram_gnd('HPb', 6, [1:5], [4,9,16], 1);



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
% if nargin<5,
%     do_wrtgnd=1; % Whether to also do with respect to ground
% end


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

savetag = 'tmp';

if params.fpass(2) == 400
    savetag = '';
    %movingwin = [100 10]/1000; 
end
if params.fpass(2) == 100
    savetag = 'mid';
    %movingwin = [400 40]/1000;
end
if params.fpass(2) == 40
    savetag = 'low';
    %movingwin = [1000 100]/1000;
end
if params.fpass(2) <= 10
    savetag = 'floor';
    %movingwin = [8000 800]/1000; 
end



% SET DATA
% -------------------------------------------

switch prefix
    case 'HPa'
        rawdir = '/data25/sjadhav/HPExpt/HPa/';
        directoryname = '/data25/sjadhav/HPExpt/HPa_direct/';
    case 'HPb'
        rawdir = '/data25/sjadhav/HPExpt/HPb/';
        directoryname = '/data25/sjadhav/HPExpt/HPb_direct/'; 
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
        eeggndspec = []; dummy_gnd = []; % Save 1 file for each tet in a day 
         
        flaggnd=0;
        for ep=1:length(epochs)
            eeggnd=[];
            epoch=epochs(ep);
            
            % Do the gnd channel         
            curreeggndfile = [dir2,'/EEG/',prefix,'eeggnd',daystring,'-',num2str(epoch),'-',tetstring];
            if (exist([curreeggndfile,'.mat'],'file'))==2
                flaggnd=1;
                load(curreeggndfile);
                specgram_gnd = mtspecgramc(eeggnd{day}{epoch}{tet}.data,movingwin,params);
                % Dont save specgram in file - its too big
                % --------------------------------------------
                eeggndspec{day}{epoch}{tet}.meanspec = mean(specgram_gnd,1);
                eeggndspec{day}{epoch}{tet}.stdspec = std(specgram_gnd,1);
                dummy_gnd=[dummy_gnd;specgram_gnd];
                % Dont save specgram in file - its too big
                %eeggndspec{day}{epoch}{tet}.specgram  = specgram_gnd;
                eeggndspec{day}{1}{tet}.meandayspec=mean(dummy_gnd,1);
                eeggndspec{day}{1}{tet}.stddayspec=std(dummy_gnd,1);
                savefile = [savedir,prefix,'eeggndspec',savetag,daystring,'-Tet',tetstring]
                save(savefile,'eeggndspec');
            else
                disp(['EEG wrt Gnd file seems to not exist for Tet ',tetstring]);
                disp(['Tet ',tetstring, ' must be a Ref. Use normal eeg file spectrogram']);
            end
            
        end % end epochs
    end % end tets
  
end % end days





i=1;
%set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.03 Screen(3)*0.2 Screen(4)*0.4])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



