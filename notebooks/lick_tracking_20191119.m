

clockrate=30000;
likethresh = .999;
animal = 'YM3';
%% load data
% extract script:
% home/droumis/Src/dataEngineering/extract_rec/DR_extract_script_LickTracking.py
andef = animaldef('ym3');
date = '20191107';
file = dir([andef{end} '*.csv']);
fprintf('importing::: %s\n',file.name)
DLC_tracking = csvread([file.folder '/' file.name],3);

% Get timestamps from .rec file
timedir = '20191016_Lotus_03_lineartrack.time/';
basenmrec='20191016_Lotus_03';


timefile=dir([andef{5} date '/*.time/*continuoustime.dat']);
ptp_ctime_filename = [timefile.folder '/' timefile.name];
ptp_ctime=readTrodesExtractedDataFile(ptp_ctime_filename);
tt = double(ptp_ctime.fields(1).data);
% st = double(ptp_ctime.fields(2).data);

%%
% Get camera timestamps
camHWlog_file = dir([andef{4} date '/*.cameraHWSync']);
camHWlog_filename = [camHWlog_file.folder '/' camHWlog_file.name];
camHWlog = readTrodesExtractedDataFile(camHWlog_filename);

t_camHW_ns = double(camHWlog.fields(3).data);
n_camHW_framecount = camHWlog.fields(2).data;
n_dropped_frames = sum(diff(n_camHW_framecount)-1);

t_camHW_s = t_camHW_ns/1e9;
t_camHW_ms = t_camHW_ns/1e6;

fprintf('# of CamHW timestamps: %d \n', length(t_camHW_ms));
fprintf('# of dropped frames: %d \n', n_dropped_frames);
est_framerate = median(1./diff(t_camHW_ns./1e9));
fprintf('Estimated frame rate from camHW: %0.3f \n', est_framerate);
fprintf('First/Last record: %0.3f %0.3f (%0.3f elapsed) \n', t_camHW_s(1), t_camHW_s(end), diff(t_camHW_s([1 end])));

%% regress cam ptp time to trodes time
XX = tt;
Y = st;
bls = polyfit(XX,Y,1);

% 'reconstructed' sysClock times of each Trodes packet
st_fit = tt * bls(1) + bls(2);
cam_tt_fit = (t_camHW_ns -bls(2)) ./ bls(1);
cam_rt_fit=cam_tt_fit./clockrate;
median(diff(cam_rt_fit(1:1000))) %= 8 ms
