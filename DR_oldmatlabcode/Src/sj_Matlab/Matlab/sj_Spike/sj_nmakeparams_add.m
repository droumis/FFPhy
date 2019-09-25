function sj_nmakeparams_add(datafile,tetnum, waves, timestamps, THRESH, MAXALLOWEDAMP, dopc, Samprange, donoise, doaddnspikefile)
%
% Usually called from sj_makedayparams_dio
%
% Changing sj_nmakeparams_withpc - Adding params like Width, Peak, Valley, and changing PC / Removing PC
% Have to dejitter = peak align for this
% sj_nmakeparams_withpc - Change nmakeparams to calculate PCs and also make .fet file for
% Klusatakiwk and .spk file for Klusters, and gspikefile for gsort
% Shantanu, 04/21/2010

%nmakeparams(datafile,tetnum, waves, timestamps)
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

if (nargin < 5)
    THRESH = 0;  %threshhold: at least one spike must excede this amplitude (in micro volts)
end
if (nargin < 6)
    MAXALLOWEDAMP = 1000;  %maximum allowed amplitude in micro volts.  This filters out noise and makes the graphing funtions faster
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
Fs = SamplingRate/1000; % in kHz
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

% Tetrode Info
%read the .dat.config files in the directory above to get the gain and
%threshhold info
tetrodeinfo = getTetrodeConfig2;
posinfo = [];
times = [];
pnames = dir('*.p');

% Position - load all the position data from the directory above (uses all .p files)
for i = 1:length(pnames)
    [tmptimes, tmppos] = getposinfo(pnames(i).name);
    times = [times; tmptimes];
    posinfo = [posinfo; tmppos];
end
[times, I] = sort(times);
posinfo = posinfo(I,:);

% Set Up
cd(currdir);
%gains = tetrodeinfo(tetnum).gain;
thresh = tetrodeinfo(tetnum).thresh;
%thresh = round(((thresh/1000)*32767)/MAXAMP);
%gains(find(gains < 5)) = inf;  %if the gain is zero, then change it to
%infinity so the the division below won't cause events to equal inf
%timestamps = [];
%waves = [];
tmp = [];
filedata.params = zeros(length(tmp),5);
%loadload(datafile)
progress = 0;
%triggers = max((squeeze(waves(9,:,:)))')
triggers = thresh;
%waves = waves(:,:,1);

% Noise by Amplitude Plot Thresh if asked for
if donoise ==1
    [waves, timestamps, noiseidx] = sss_windowcut_waves_hist(waves, timestamps, nch);
end

% Calculate Base Parameters

% parmcalc calculates the waveform paramters and returns a matrix with 6
% columns for peak1, peak2, peak3, peak4, maximumwidth, maxheight change
% (from peak to trough)

nch = size(waves,2);

disp('      Base Parameter Calc');
filedata.params = parmcalc(timestamps,waves,triggers)';
% change the peak amplitudes into microvolts
%filedata.params(:,1:4) = filedata.params(:,1:4) * diag(repmat((5/32767),1,4)./gains);
%filedata.params(:,1:4) = filedata.params(:,1:4)*1000000;


% Remove Spikes Not In Range / Not Good
goodspikes = find( ((filedata.params(:,1)<MAXALLOWEDAMP)&(filedata.params(:,2)<MAXALLOWEDAMP)&(filedata.params(:,3)<MAXALLOWEDAMP)&(filedata.params(:,4)<MAXALLOWEDAMP)) & ...
    ((filedata.params(:,1)>THRESH)|(filedata.params(:,2)>THRESH)|(filedata.params(:,3)>THRESH)|(filedata.params(:,4)>THRESH)) );

% GET ADDITIONAL PARAMETERS
disp('      Addn Parameter Calc');
% Additional Parameters:  1)Peak-Val Ratio, 2)Peak-Val Latency, 3)Width-FWHM 4)Valley, 5)Valley a, 6) if dopc=1, 3xNch PCs, usually 12  
[addn_params] = sj_addparams(waves,triggers,dopc);

disp('      Making File Data');

%tic
% Make filedata
if ~isempty(goodspikes)
    filedata.params = filedata.params(goodspikes,:);
    timestamps = timestamps(goodspikes);
    addn_params = addn_params(goodspikes,:);
    
    % add some jitter to the width points (because of poor sampling), this decreases point overlap
    addnoise = 2 * rand(size(filedata.params,1),1)-1; % Range -1 to 1
    filedata.params(:,5) = filedata.params(:,5)+addnoise;
    %also, add noise to the amplitude data
    addnoise = rand(size(filedata.params,1),4)-.5;
    filedata.params(:,1:4) = filedata.params(:,1:4)+addnoise; % Range -0.5 to 0.5
    
    % Addn Params - Add Noise to all except peak-valley ratio and PC, ie 4
    % params of 2)Peak-Val Latency, 3)Width-FWHM, 4)Valley, 5)Valley a
    
    % Add Noise to Peak-Valley Ratio
%     addnoise = rand(size(filedata.params,1),4)-.5; % Range -0.5 to 0.5
%     addn_params(:,1:4) = addn_params(:,1:4) + addnoise;
    
    % Add Noise to Peak-Valley Latency and Width
    addnoise = 2 * rand(size(addn_params,1),8)-1;
    addn_params(:,5:12) = addn_params(:,5:12) + addnoise; % Range -1 to 1
    
    % Add Noise to Valley and valleya
    addnoise = rand(size(filedata.params,1),8)-.5; % Range -0.5 to 0.5
    addn_params(:,13:20) = addn_params(:,13:20) + addnoise; 

    % Add Noise to PCs
    
    addnoise = [];
    
    spikepos = double(posinfo(lookup(double(timestamps),double(times)),:));
    %just in case some points were bad, we set them to 0
    spikepos(find(spikepos < 0)) = 0;
    spikepos(find(spikepos > 1000)) = 0;
    
    % we also add some jitter to the position info to decrease point overlap
    noise = .002;  %percentage of noise to add to the position info
    ma = max(spikepos(:,1:2));
    mi = min(spikepos(:,1:2));
    addnoise = randn(size(spikepos,1),2);   % Standard Normal Distr with Mean 0
    addnoise(:,1) = addnoise(:,1)*((ma(1)-mi(1))*noise);
    

    % The final parameter matrix:
    % Save One with all Params and One for Matclust
    if dopc==1
        
        % All Params File      
        % Time + 6 base params + 2pos + 20 addnparams + 12 PCs
        filedata.params = [double(timestamps) filedata.params(:,1:6) spikepos(:,1:2)+addnoise addn_params];
        
        % we also need to give the names of these parameters
        filedata.paramnames = {'Time',
            'Channel 1 Max',
            'Channel 2 Max',
            'Channel 3 Max',
            'Channel 4 Max',
            'Max width',
            'Max height change'
            'X position',
            'Y position',
            'PeakValRatio1',
            'PeakValRatio2',
            'PeakValRatio3',
            'PeakValRatio4',
            'PeakValLat1',
            'PeakValLat2',
            'PeakValLat3',
            'PeakValLat4',   
            'Width1',
            'Width2',
            'Width3',
            'Width4',
            'Valley1',
            'Valley2',
            'Valley3',
            'Valley4',
            'Valleya1',
            'Valleya2',
            'Valleya3',
            'Valleya4',
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
        eval(['save ',datafile,'_allparams filedata -V6']);
        
        %----------------------------------------------------
        
        % Param File for Matclust
        % {Time + 6 base params + 2pos} already in + 4 addnparams + 4PCs + 4PCs
        filedata.params = [filedata.params(:,1:9) addn_params(:,13:16) ...
            addn_params(:,21) addn_params(:,24) addn_params(:,27) addn_params(:,30)];
        
        % we also need to give the names of these parameters
        filedata.paramnames = {'Time',
            'Channel 1 Max',
            'Channel 2 Max',
            'Channel 3 Max',
            'Channel 4 Max',
            'Max width',
            'Max height change'
            'X position',
            'Y position',
            'Valley1',
            'Valley2',
            'Valley3',
            'Valley4',
            'Wave1PC1',
            'Wave2PC1',
            'Wave3PC1',
            'Wave4PC1'};
        
        % Waveform Data File matclust also needs the name of the waveform data file
        filedata.filename = datafile;
        eval(['save ',datafile,'_params filedata -V6']);     
        
    else  % No PC
        
        % All Params
        % Time + 6 base params + 2pos + 20 addnparams 
            filedata.params = [double(timestamps) filedata.params(:,1:4) spikepos(:,1:2)+addnoise addn_params];
                    
        % we also need to give the names of these parameters
       filedata.paramnames = {'Time',
            'Channel 1 Max',
            'Channel 2 Max',
            'Channel 3 Max',
            'Channel 4 Max',
            'Max width',
            'Max height change'
            'X position',
            'Y position',
            'PeakValRatio1',
            'PeakValRatio2',
            'PeakValRatio3',
            'PeakValRatio4',
            'PeakValLat1',
            'PeakValLat2',
            'PeakValLat3',
            'PeakValLat4',   
            'Width1',
            'Width2',
            'Width3',
            'Width4',
            'Valley1',
            'Valley2',
            'Valley3',
            'Valley4',
            'Valleya1',
            'Valleya2',
            'Valleya3',
            'Valleya4'};
        
        
          % All Param File
        filedata.filename = datafile;
        eval(['save ',datafile,'_allparams filedata -V6']);
        
        % --------------------------------------------------------
        
         % Param File for Matclust
        % {Time + 6 base params + 2pos} already in + 4 addnparams 
        filedata.params = [filedata.params(:,1:9) addn_params(:,13:16)];
        
        % we also need to give the names of these parameters
         filedata.paramnames = {'Time',
            'Channel 1 Max',
            'Channel 2 Max',
            'Channel 3 Max',
            'Channel 4 Max',
            'Max width',
            'Max height change'
            'X position',
            'Y position',
            'Valley1',
            'Valley2',
            'Valley3',
            'Valley4'};
        
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
 
% if doaddnspikefile~=0
%     clear wavesli
%     clear spikepos;
%     % Make a feature file (.fet) for Klustakwik and Klusters
%     % Order as 12 PCs, 4 amps, 1 Kkwik time
%     if ~isempty(goodspikes)
%         
%         idxs = [10 11 12 13 14 15 16 17 18 19 20 21 2 3 4 5 22];
%         nfets = length(idxs); % its always 12+4+1=17
%         fet_matrix = zeros(length(timestamps),length(idxs));
%         fet_matrix = filedata.params(:,idxs);
%         
%         % Write to ascii file
%         %    Meth 1 - use Matlab save
%         %    save(fetfile,'nfets','fet_matrix','-ascii');
%         
%         % Meth 2 - use fprintf
%         cd Kkwik
%         fid = fopen(fetfile, 'wt');
%         fprintf(fid, '%2.0f\n', nfets);
%         fprintf(fid, '%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.2f %8.2f %8.2f %8.2f %8.2f\n', fet_matrix');
%         fclose(fid);
%         cd ..
%         
%         %% Save in gspikefile format
%         if dogspikefile==1
%             spikes.waveforms=[];
%             for i=1:nch
%                 cmd=sprintf('spikes.waveforms_ch%d=double(squeeze(waves(:,%d,:))\'');',i,i); eval(cmd);
%                 cmd=sprintf('spikes.waveforms = [spikes.waveforms, spikes.waveforms_ch%d];',i); eval(cmd);
%             end
%             spikes.Fs=Fs;
%             spikes.swtimes=double(timestamps)./UnitsPerSec; % in sec
%             spikes.fstimes=spikes.swtimes*1000;    % in ms
%             spikes.ftimes=int16(spikes.fstimes);
%             spikes.spiketimes=spikes.swtimes;
%             spikes.threshT=9;
%             spikes.threshV(1:nch,1)=-Inf;
%             spikes.threshV(1:nch,2)=[thresh'];
%             spikes.noiseidx = noiseidx;
%             spikes.matclustparams = filedata.params(:,1:9);
%             for i=1:9, spikes.matclustparamnames{i} = filedata.paramnames{i}; end
%             
%             cd Gsort
%             save(gspikefile,'spikes');
%             cd ..
%         end
%         
%         clear filedata;
%     end
% end



