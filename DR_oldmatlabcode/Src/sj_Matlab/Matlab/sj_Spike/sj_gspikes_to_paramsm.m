function sj_gspikes_to_paramsm(gspikefile, nch, elec, Nspike, dopos, domatclust)
% sj_gspikes_to_paramsm('08_030610_tet7_gspike_tetNF2_sortf4', 4, '7', 1,1, 1);
% Use the spikes structure in gsort and generate parameters for Matclust: 
% Save in spikes structure as well as  for KKwik and
% also as matclust .mat file and matclust _params.mat file if asked
% for

%%%%%%%%%%%%%%%
if (nargin < 2)
    nch = 4;
end

if (nargin < 3)
    elec = '1';
end

if (nargin < 4)
    Nspike = 1; % when data comes from Nspike
end
if nargin<5,
    dopos=1;
end

if nargin<6,
    domatclust=0;
end
 


%%%%%%%%%%%%%%
load(gspikefile);
Samprange=[4:32];

load generic_wave34;  % in /home/Src/Matlab

%%%%%% Time stamps in xx.y ms - resolution of 0.1ms
Fs = spikes.Fs; % in kHz
SamplingRate = Fs*1000; %sampling rate within each waveform
%TimePerStep = (1/SamplingRate)*UnitsPerSec;
%TimePerWindow = TimePerStep*40;
%MAXAMP = 5;     %max voltage (in milivolts)

%%%%%%%%% Get File Names 
new_gspikefile = [gspikefile '_params'];
matclust_file = [gspikefile '_matc'];
matclust_paramfile = [gspikefile '_matc_params'];


%%%%%%%%% Calc Params

currdir = pwd;
cd ../..
    
tetrodeinfo = getTetrodeConfig2;
tetnum=str2num(elec);
thresh = tetrodeinfo(tetnum).thresh;
triggers = thresh;
waves=int16(zeros(size(spikes.waveforms_ch1,2),nch,size(spikes.waveforms_ch1,1)));
for i=1:nch
    if triggers(i)<200 % if channel is valid
        cmd=sprintf('waves(:,%d,:)=int16(spikes.waveforms_ch%d\'');',i,i); eval(cmd);
    else
        waves(:,i,:) = repmat(generic,1,size(spikes.waveforms_ch1,1));
    end
end
timestamps=uint32(spikes.fstimes*10);
filedata.params = parmcalc(timestamps,waves,triggers)';
% add noise to the amplitude data
addnoise = rand(size(filedata.params,1),4)-.5;
filedata.params(:,1:4) = filedata.params(:,1:4)+addnoise;
% add some jitter to the width points (becuase of poor sampling, this
% decreases point overlap
addnoise = 2 * rand(size(filedata.params,1),1)-1;
filedata.params(:,5) = filedata.params(:,5)+addnoise;
filedata.params = [double(timestamps) filedata.params];



% Position Data

if dopos==1,  

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
    
    filedata.params = [filedata.params spikepos(:,1:2)+addnoise];
end


cd(currdir);

if ~isfield(spikes,'matclustparams'),
    spikes.matclustparams = filedata.params;
    spikes.matclustparamnames = {'Time',
        'Channel 1 Max',
        'Channel 2 Max',
        'Channel 3 Max',
        'Channel 4 Max',
        'Max width',
        'Max height change';
        'X position';
        'Y position'};
else
    disp('Matclust Params already exist - Removing field to recalculate');
    rmfield(spikes,'matclustparams');
    spikes.matclustparams = filedata.params;
    spikes.matclustparamnames = {'Time',
        'Channel 1 Max',
        'Channel 2 Max',
        'Channel 3 Max',
        'Channel 4 Max',
        'Max width',
        'Max height change';
        'X position';
        'Y position'};
end

save(new_gspikefile,'spikes');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matclust format
if domatclust==1,
    
    % 1st file
    if ~exist('timestamps','var'),
        timestamps=uint32(spikes.fstimes*10);
    end
    save(matclust_file,'timestamps','waves');
    
    % 2nd file
    
    filedata.paramnames = {'Time',
        'Channel 1 Max',
        'Channel 2 Max',
        'Channel 3 Max',
        'Channel 4 Max',
        'Max width',
        'Max height change',
        'X position',
        'Y position'};
    
    filedata.filename = gspikefile;
    save(matclust_paramfile,'filedata');
 
end







