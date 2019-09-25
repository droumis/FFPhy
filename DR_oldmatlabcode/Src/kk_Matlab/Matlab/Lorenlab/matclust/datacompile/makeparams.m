function makeparams(datafile,tetnum)

% this program makes the parameter file that matclust can open.
% datafile is the name of the .mat file with the waveform data (do not include
% the '.mat' in the name). Tetnum is the tetrode number.  This should be run from 
% the directory that has the waveform data.  The program saves a parameter file in the same
% directory titled datafile_params.mat.  The program also computes position
% as a parameter, which means that the .p files for all epochs for that day
% must be stored in the directory above the current directory.
% example:  makeparams('rat02-12',12) -> output file:  rat02-12_params.mat

UnitsPerSec = 10000;  %for the .tt timestamps, every 10000 equals one second
SamplingRate = 31200; %sampling rate within each waveform
TimePerStep = (1/SamplingRate)*UnitsPerSec;
TimePerWindow = TimePerStep*40;


load(datafile);

currdir = pwd;
cd ..
%read the .dat.config files in the directory above to get the gain and
%threshhold info
tetrodeinfo = getTetrodeConfig;
posinfo = [];
times = [];
pnames = dir('*.p');
%load all the position data from the directory above (uses all .p files)
for i = 1:length(pnames)
   [tmptimes, tmppos] = getposinfo(pnames(i).name);
   times = [times; tmptimes];
   posinfo = [posinfo; tmppos];
end
[times, I] = sort(times);
posinfo = posinfo(I,:);


cd(currdir);
gains = tetrodeinfo(tetnum).gain;
thresh = tetrodeinfo(tetnum).thresh;
gains(find(gains < 5)) = inf;  %if the gain is zero, then change it to infinity so the the division below won't cause events to equal inf


timestamps = [];
waves = [];
filedata.params = zeros(length(timestamps),5);
load(datafile)

progress = 0;

%triggers = max((squeeze(waves(9,:,:)))')
triggers = thresh;

%waves = waves(:,:,1);

% parmcalc calculates the waveform paramters and returns a matrix with 6
% columns for peak1, peak2, peak3, peak4, maximumwidth, maxheight change
% (from peak to trough)
filedata.params = parmcalc(timestamps,waves,triggers)';

% change the peak amplitudes into microvolts 
filedata.params(:,1:4) = filedata.params(:,1:4) * diag(repmat((5/32767),1,4)./gains);
filedata.params(:,1:4) = filedata.params(:,1:4)*1000000;

% add some jitter to the width points (becuase of poor sampling, this
% decreases point overlap
addnoise = 2 * rand(size(filedata.params,1),1)-1;
filedata.params(:,5) = filedata.params(:,5)+addnoise;


addnoise = [];

spikepos = double(posinfo(lookup(double(timestamps),double(times)),:));
%just in case some points were bad, we set them to 0
spikepos(find(spikepos < 0)) = 0; 
spikepos(find(spikepos > 1000)) = 0; 



% we also add some jitter to the position info to decrease point overlap
noise = .002;  %percentage of noise to add to the position info
ma = max(spikepos(:,1:2));
mi = min(spikepos(:,1:2));
addnoise = randn(size(spikepos,1),2);
addnoise(:,1) = addnoise(:,1)*((ma(1)-mi(1))*noise);

%the final parameter matrix has 9 columns
filedata.params = [double(timestamps) filedata.params spikepos(:,1:2)+addnoise];


% we also need to give the names of these parameters
filedata.paramnames = {'Time', 
              'Channel 1 Max',
              'Channel 2 Max',
              'Channel 3 Max',
              'Channel 4 Max',
              'Max width',
              'Max height change'
              'X position',
              'Y position'};
              
% matclust also needs the name of the waveform data file
filedata.filename = datafile;

eval(['save ',datafile,'_params filedata -V6']);



clear spikepos;
clear filedata;



