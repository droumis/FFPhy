function sj_nmakeparamsNoPos_withpc(datafile,tetnum, waves, timestamps, THRESH, MAXALLOWEDAMP, dogspikefile, Samprange)
%
% Usually called from sj_makedayparams_withpc
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
    dogspikefile = 0;  %Save in gspikefile format or not
end

if (nargin < 8)
    Samprange = [4:32];
end


UnitsPerSec = 10000;  %for the .tt timestamps, every 10000 equals one second
SamplingRate = 31200; %sampling rate within each waveform
TimePerStep = (1/SamplingRate)*UnitsPerSec;
TimePerWindow = TimePerStep*40;
%MAXAMP = 5;     %max voltage (in milivolts)

% Get Name for Feature File for Kkwik
[day,elec] = strtok(datafile,'-');
elec = elec(2:end);
if elec(1)=='0',
    elec=elec(2:end);
end
fetfile = [day '.fet.' elec];
spkfile = [day '.spk.' elec];
gspikefile = [day '_tet' elec '_gspike'];

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

% parmcalc calculates the waveform paramters and returns a matrix with 6
% columns for peak1, peak2, peak3, peak4, maximumwidth, maxheight change
% (from peak to trough)

filedata.params = parmcalc(timestamps,waves,triggers)';

% change the peak amplitudes into microvolts 
%filedata.params(:,1:4) = filedata.params(:,1:4) * diag(repmat((5/32767),1,4)./gains);
%filedata.params(:,1:4) = filedata.params(:,1:4)*1000000;


goodspikes = find( ((filedata.params(:,1)<MAXALLOWEDAMP)&(filedata.params(:,2)<MAXALLOWEDAMP)&(filedata.params(:,3)<MAXALLOWEDAMP)&(filedata.params(:,4)<MAXALLOWEDAMP)) & ...
                   ((filedata.params(:,1)>THRESH)|(filedata.params(:,2)>THRESH)|(filedata.params(:,3)>THRESH)|(filedata.params(:,4)>THRESH)) );

% GET PCs 
spc = sj_pcasvd_matspikes(waves,size(waves,2),Samprange);
%filedata.params(:,10:21)=spc; %Dont Put In Yet

if ~isempty(goodspikes)
   filedata.params = filedata.params(goodspikes,:);
   timestamps = timestamps(goodspikes);
   % add some jitter to the width points (becuase of poor sampling, this
   % decreases point overlap
   addnoise = 2 * rand(size(filedata.params,1),1)-1;
   filedata.params(:,5) = filedata.params(:,5)+addnoise;
   %also, add noise to the amplitude data
   addnoise = rand(size(filedata.params,1),4)-.5;
   filedata.params(:,1:4) = filedata.params(:,1:4)+addnoise;
   
   addnoise = [];
   

   
   %the final parameter matrix has 7 columns + 12 PCs + Time in Kkwik format (Normalized to start at 0)
   % Time in Kwik/Klusters format in multiples of sampling interval (Fs*Time in Secs)
   
   filedata.params = [double(timestamps) filedata.params spc (Fs.*(double(timestamps-timestamps(1))./UnitsPerSec))];
   
   
   % we also need to give the names of these parameters
   filedata.paramnames = {'Time', 
               'Channel 1 Max',
               'Channel 2 Max',
               'Channel 3 Max',
               'Channel 4 Max',
               'Max width',
               'Max height change',
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
               'Wave4PC3', 
               'KkwikTime'};
               
               
   % matclust also needs the name of the waveform data file
   filedata.filename = datafile;
   
   eval(['save ',datafile,'_params filedata -V6']);
   
   waves = waves(:,:,goodspikes);
   save(datafile,'waves','timestamps');
   
   %% Save binary file for Klusters
   % Save as single column, Samp1Elec1Spk1, Samp1,Elec2Spk1... Samp40Elec4Spk1 
   wavesli = permute(waves,[2 1 3]); 
   wavesli = wavesli(:); 
   fid = fopen(spkfile, 'w'); % Do I need big endian or native or ... Check later
   fwrite(fid, wavesli, 'int16'); % Need 16 or 32 bit binary, not default uint8
   fclose(fid);
   
end

clear wavesli
clear spc;

% Make a feature file (.fet) for Klustakwik and Klusters
% Order as 12 PCs, 4 amps, 1 Kkwik time
if ~isempty(goodspikes)
    
    idxs = [8 9 10 11 12 13 14 15 16 17 18 19 2 3 4 5 20];
    nfets = length(idxs); % its 12+4+1=17
    fet_matrix = zeros(length(timestamps),length(idxs));
    fet_matrix = filedata.params(:,idxs);
    
    % Write to ascii file  
%    Meth 1 - use Matlab save 
%    save(fetfile,'nfets','fet_matrix','-ascii');
    
    % Meth 2 - use fprintf
     fid = fopen(fetfile, 'wt');
     fprintf(fid, '%2.0f\n', nfets);
     fprintf(fid, '%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.2f %8.2f %8.2f %8.2f %8.2f\n', fet_matrix');
     fclose(fid);
     
        
     clear filedata;
     
      %% Save in gspikefile format
   if dogspikefile==1
       spikes.waveforms_ch1=squeeze(waves(:,1,:))';
       spikes.waveforms_ch2=squeeze(waves(:,2,:))';
       spikes.waveforms_ch3=squeeze(waves(:,3,:))';
       spikes.waveforms_ch4=squeeze(waves(:,4,:))';
       spikes.Fs=Fs;
       spikes.swtimes=double(timestamps)./UnitsPerSec; % in sec
       spikes.fstimes=spikes.swtimes*1000;    % in ms
       spikes.ftimes=int16(spikes.fstimes);
       spikes.spiketimes=spikes.swtimes;
       spikes.threshT=9;
       spikes.threshV=[-Inf min(thresh)];
       spikes.waveforms=[spikes.waveforms_ch1, spikes.waveforms_ch2, spikes.waveforms_ch3, spikes.waveforms_ch4];
   
       save(gspikefile,'spikes');
   end
    
    
end





