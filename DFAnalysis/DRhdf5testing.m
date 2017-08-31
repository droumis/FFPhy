


%% step 1: run TrodestoMatlab to get raw continuous
%out = readTrodesFileContinuous(filename,channels,skipTime,varargin)
%%make sure to use version  with the new field 'skiptime'.. set to 1 for
%%D05 data

%% step 2: run gapcheck to make sure no gaps in timestamps


%% step 3: run convertTrodestoSD to get .dat

%[Time1Channel1, Time1Channel2, .. Time1ChannelN, Time2Channel1, ...]
combtetdata = reshape(dataz', [numel(dataz), 1]);

%write the file flattened to a signed 2byte (int16)
fid = fopen('D09_d06_e02_t19.dat','w'); %open file for writing
fwrite(fid, combtetdata, 'int16');
fclose(fid)

%% step 4: run phy detect myfile.prm to get Kwik and Kwx files.. make prm,prb files manually
%in sys terminal
%source activate phy
%phy detect file.prm

%% step 5: grab the spiketimes from kwik and the wavs from kwx

%%waveforms from the KWX file
wavsall = h5read('D09_d06_e02_t19_test1.kwx','/channel_groups/0/waveforms_filtered');

%%spiketimes from the Kwik file
timesamples = h5read('D09_d06_e02_t19_test1.kwik','/channel_groups/0/spikes/time_samples');

%% step 6: get the peak amplitudes from wavs

peakamps = max(wavsall,[],2);

%% step 7: pos reconstruct in trodes cam mod to get xy position.. load in XY
%currently i need to run cameramodule from qt.. i dunno why is doesnt
%execute from terminal of file explorer
%fuck there's something wrong with this.... camera module isn't saving
%tracking file..
posfields = readCameraModulePositionFile('D09_d06_e02_12-12-2015(11_05_45).videoPositionTracking');
postimestamps = readCameraModuleTimeStamps('D09_d06_e02_12-12-2015(11_05_45).videoTimeStamps');
%% step 8: put the spiketimes, wavs, peaks, timestamps, posxy together 










%% scrap

info = h5info('test_hybrid_120sec.kwik')
myinfo = h5info('D09_d06_e02_t19_test1.kwik');
%%
h5disp('test_hybrid_120sec.kwik')
%%

data = h5read('test_hybrid_120sec.kwik',datasetname)

%get the time samples
timesamples = h5read('test_hybrid_120sec.kwik','/channel_groups/0/spikes/time_samples')
timesamples = h5read('D09_d06_e02_t19_test1.kwik','/channel_groups/0/spikes/time_samples')


amplitudes = h5read('test_hybrid_120sec.kwik','/channel_groups/0/clusters/main/10/quality_measures/amplitude')


clustwhat = h5read('test_hybrid_120sec.kwik','/channel_groups/0/spikes/clusters/main')

%%waveforms from the KWX file
wavsall = h5read('test_hybrid_120sec.kwx','/channel_groups/0/waveforms_filtered');
close all; wav1 = wavsall(22,:,505); plot(wav1)

