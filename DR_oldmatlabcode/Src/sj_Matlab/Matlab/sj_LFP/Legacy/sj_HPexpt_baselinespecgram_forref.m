function [params] = sj_HPexpt_baselinespecgram_forref(prefix, days, epochs, tets, reftag)
% Shantanu - Nov 2012. 
% To avoid errors, make a separate file for wrt gnd spectrogram and ref spectrogram
% For this file, have to specify a Tag which tells waht the Reference is
% Analogous code in sj_HPexpt_baselinespecgram.m and sj_HPexpt_baselinespecgram_forgnd.m respectively

% From Kenny and Maggies - event_spectrograms .m and calcriptriggerredspectrograms.m respectively
% This is for getting baseline values for given eeg tets for normalization when needed. Save These. 

% sj_HPexpt_baselinespecgram_forref('HPa', 8, [1:5], [15:20],'Ref19Pfc');
% sj_HPexpt_baselinespecgram_forref('HPa', 2, [1:5], [15:20],'Ref19Pfc');

% sj_HPexpt_baselinespecgram_forref('HPb', 1, [1:7], 8:14,'Ref7Hp');
% sj_HPexpt_baselinespecgram_forref('HPb', 2, [1:5], 8:14, 'Ref7Hp');
% sj_HPexpt_baselinespecgram_forref('HPb', 5, [1:5], 8:14, 'Ref7Hp');
% sj_HPexpt_baselinespecgram_forref('HPb', 6, [1:5], 8:14, 'Ref7Hp');



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
if nargin<5,
    reftag='Ref7Hp'; % Whether to also do with respect to ground
end


% Define Chronux params
% -------------------------------------------
movingwin = [100 10]/1000;          
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];

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
            curreegreffile = [dir2,'/EEG/',prefix,'eeg',reftag,daystring,'-',num2str(epoch),'-',tetstring];
            if (exist([curreegreffile,'.mat'],'file'))==2
                flaggnd=1;
                load(curreegreffile);
                specgram_gnd = mtspecgramc(eegref{day}{epoch}{tet}.data,movingwin,params);
                eegrefspec{day}{epoch}{tet}.meanspec = mean(specgram_gnd,1);
                eegrefspec{day}{epoch}{tet}.stdspec = std(specgram_gnd,1);
                dummy_gnd=[dummy_gnd;specgram_gnd];
                % Dont save specgram in file - its too big
                %eeggndspec{day}{epoch}{tet}.specgram  = specgram_gnd;
                eegrefspec{day}{1}{tet}.meandayspec=mean(dummy_gnd,1);
                eegrefspec{day}{1}{tet}.stddayspec=std(dummy_gnd,1);
                savefile = [savedir,prefix,'eeg',reftag,'spec',daystring,'-Tet',tetstring]
                save(savefile,'eegrefspec');
            else
                disp(['EEG wrt Given Ref file seems to not exist for Tet ',tetstring]);
            end
            
        end % end epochs
    end % end tets
  
end % end days





i=1;
%set(gcf,'Position',[Screen(3)*0.205 Screen(4)*0.03 Screen(3)*0.2 Screen(4)*0.4])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



