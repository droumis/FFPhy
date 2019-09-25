function sj_nmakeparamsNoPos_pc(datafile,tetnum, waves, timestamps, THRESH, MAXALLOWEDAMP, dopc, Samprange, donoise, doaddnspikefile)
%
% Usually called from sj_makedayparams_pc
% From sj_nmakeparmas_add. Getting rid of additional parameters, and only keeping PCs
% Shantanu - May 2012
%
%
% Change nmakeparamsNoPos to calculate PCs and also make .fet file for
% Klusatakiwk and .spk file for Klusters
% Shantanu, 04/21/2010

% nmakeparams(datafile,tetnum, waves, timestamps)
%nmakeparams(datafile,tetnum, waves, timestamps, THRESH, MAXALLOWEDAMP)
%USE WITH NEW NSPIKE RIGS
% this program makes the parameter file that matclust can open.
% datafile is the name of the .mat file with the waveform data (do not include
% the '.mat' in the name). Tetnum is the tetrode number.  
% Waves is an 3 dimensional matrix containing the spike waveforms. It is 
% spikelength-by-numchannels-by-N. Timestamps is a vector with length N containing
%the time of each spike. Thrsh and maxallowedamp are optional.  They describe the spike threshhold
%(in at least one channel) and the maximum allowed amplitude (accross all channels)
%for a particular spike to be included in the parameters. Defaults are 0 and 2500 (microvolts). This should be run from 
% the directory that has the waveform data.  The program saves a parameter file in the same
% directory titled datafile_params.mat, and the waves and timestamps in datafile.  The program also computes position
% as a parameter, which means that the .p files for all epochs for that day
% must be stored in the directory above the current directory.
% example:  makeparams('rat02-12',12) -> output file:  rat02-12_params.mat

if (nargin < 6)
    THRESH = 0;  %threshhold: at least one spike must excede this amplitude (in micro volts)
end
if (nargin < 5)
    MAXALLOWEDAMP = 2500;  %maximum allowed amplitude in micro volts.  This filters out noise and makes the graphing funtions faster 
end

if (nargin < 7)
    dopc = 0;
end

if (nargin < 8)
    Samprange = [4:32];
end

if (nargin < 9)
    donoise = 0;
end

if (nargin < 10)
    doaddnspikefile = 0;  %Save in Kkwik and gspikefile format or not
end

dogspikefile=0; 


UnitsPerSec = 10000;  %for the .tt timestamps, every 10000 equals one second
SamplingRate = 31200; %sampling rate within each waveform
TimePerStep = (1/SamplingRate)*UnitsPerSec;
TimePerWindow = TimePerStep*40;
%MAXAMP = 5;     %max voltage (in milivolts)

if doaddnspikefile ~=0
    % Get Name for Feature File for Kkwik
    [day,elec] = strtok(datafile,'-');
    elec = elec(2:end);
    if elec(1)=='0',
        elec=elec(2:end);
    end
    if ~(exist('Kkwik','dir')),
        mkdir Kkwik;
    end
    if ~(exist('Gsort','dir')),
        mkdir Gsort;
    end
    fetfile = [day '.fet.' elec];
    spkfile = [day '.spk.' elec];
    gspikefile = [day '_tet' elec '_gspike'];
end


%load(datafile);

currdir = pwd;
cd ..
%read the .dat.config files in the directory above to get the gain and
%threshhold info
tetrodeinfo = getTetrodeConfig2;

cd(currdir);
%gains = tetrodeinfo(tetnum).gain;
thresh = tetrodeinfo(tetnum).thresh;
%thresh = round(((thresh/1000)*32767)/MAXAMP);
%gains(find(gains < 5)) = inf;  %if the gain is zero, then change it to infinity so the the division below won't cause events to equal inf


%timestamps = [];
%waves = [];
tmp = [];
filedata.params = zeros(length(tmp),5);
%loadload(datafile)

progress = 0;

%triggers = max((squeeze(waves(9,:,:)))')
triggers = thresh;
%waves = waves(:,:,1);

% Noise by Amplitude Plot Thresh if asked for - Shantanu
if donoise ==1
    [waves, timestamps, noiseidx] = sss_windowcut_waves_hist(waves, timestamps, nch);
end


% parmcalc calculates the waveform paramters and returns a matrix with 6
% columns for peak1, peak2, peak3, peak4, maximumwidth, maxheight change
% (from peak to trough)

filedata.params = parmcalc(timestamps,waves,triggers)';

% change the peak amplitudes into microvolts 
%filedata.params(:,1:4) = filedata.params(:,1:4) * diag(repmat((5/32767),1,4)./gains);
%filedata.params(:,1:4) = filedata.params(:,1:4)*1000000;


goodspikes = find( ((filedata.params(:,1)<MAXALLOWEDAMP)&(filedata.params(:,2)<MAXALLOWEDAMP)&(filedata.params(:,3)<MAXALLOWEDAMP)&(filedata.params(:,4)<MAXALLOWEDAMP)) & ...
                   ((filedata.params(:,1)>THRESH)|(filedata.params(:,2)>THRESH)|(filedata.params(:,3)>THRESH)|(filedata.params(:,4)>THRESH)) );

% GET ADDITIONAL PARAMETERS
disp('      Addn Parameter Calc');
% Additional Parameters:  if dopc=1, 3xNch PCs, usually 12  
[addn_params] = sj_addpc(waves,triggers,dopc);

disp('      Making File Data');


if ~isempty(goodspikes)
   filedata.params = filedata.params(goodspikes,:);
   timestamps = timestamps(goodspikes);
   addn_params = addn_params(goodspikes,:);

   % add some jitter to the width points (becuase of poor sampling, this
   % decreases point overlap
   addnoise = 2 * rand(size(filedata.params,1),1)-1;
   filedata.params(:,5) = filedata.params(:,5)+addnoise;
   %also, add noise to the amplitude data
   addnoise = rand(size(filedata.params,1),4)-.5;
   filedata.params(:,1:4) = filedata.params(:,1:4)+addnoise;
   
   addnoise = [];
   
    % The final parameter matrix:
    % Save One with all Params and One for Matclust
    if dopc==1
        
        % All PCs File      
        % Time + 6 base params + 12 PCs
        filedata.params = [double(timestamps) filedata.params(:,1:6) addn_params];
        
        % we also need to give the names of these parameters
        filedata.paramnames = {'Time',
            'Channel 1 Max',
            'Channel 2 Max',
            'Channel 3 Max',
            'Channel 4 Max',
            'Max width',
            'Max height change'
            'Wave1PC1',
            'Wave1PC2',
            'Wave1PC3',
            'Wave2PC1',
            'Wave2PC2',
            'Wave2PC3',
            'Wave3PC1',
            'Wave3PC2',
            'Wave3PC3',
            'Wave4PC1',
            'Wave4PC2',
            'Wave4PC3'};
        
        % All Param File
        filedata.filename = datafile;
        eval(['save ',datafile,'_params_allpcs filedata -V6']);
        
        %----------------------------------------------------
        
        % Param File for Matclust
        % {Time + 6 base params} already in + 4PCs 
        filedata.params = [filedata.params(:,1:7) ...
            addn_params(:,1) addn_params(:,4) addn_params(:,7) addn_params(:,10)];
        
        % we also need to give the names of these parameters
        filedata.paramnames = {'Time',
            'Channel 1 Max',
            'Channel 2 Max',
            'Channel 3 Max',
            'Channel 4 Max',
            'Max width',
            'Max height change'
            'Wave1PC1',
            'Wave2PC1',
            'Wave3PC1',
            'Wave4PC1'};
        
        % Waveform Data File matclust also needs the name of the waveform data file
        filedata.filename = datafile;
        eval(['save ',datafile,'_params_pc filedata -V6']);     
        
    else  % No PC
  
        % --------------------------------------------------------
        
         % Param File for Matclust
        % {Time + 6 base params + 2pos} already in 
        filedata.params = [filedata.params(:,1:7)];
        
        % we also need to give the names of these parameters
         filedata.paramnames = {'Time',
            'Channel 1 Max',
            'Channel 2 Max',
            'Channel 3 Max',
            'Channel 4 Max',
            'Max width',
            'Max height change'};
        
        % Waveform Data File matclust also needs the name of the waveform data file
        filedata.filename = datafile;
        eval(['save ',datafile,'_params filedata -V6']);     
        
    
    end
   
   disp('      Saving Waves');
   
   waves = waves(:,:,goodspikes);
   save(datafile,'waves','timestamps');
   
  if doaddnspikefile~=0
        % Save binary file for Klusters
        % Save as single column, Samp1Elec1Spk1, Samp1,Elec2Spk1... Samp40Elec4Spk1
        wavesli = permute(waves,[2 1 3]);
        wavesli = wavesli(:);
        cd Kkwik
        fid = fopen(spkfile, 'w'); % Do I need big endian or native or ... Check later
        fwrite(fid, wavesli, 'int16'); % Need 16 or 32 bit binary, not default uint8
        fclose(fid);
        cd ..
    end
   
end % end goodspikes

%toc



