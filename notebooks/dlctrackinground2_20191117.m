
%{
process ym2, ym3 to get ILPC

%}


clockrate=30000;
likethresh = .999;
animal = 'YM2';
%% load data
data_dir='/stelmo/demetris/lick_tracking/rawr/20191016/20191016_Lotus_03_lineartrack/';
filename = '20191016_Lotus_03_lineartrackDeepCut_resnet50_lineartrackOct17shuffle1_1030000.csv';
DLC_tracking = csvread([data_dir filename], 3);

% Get timestamps from .rec file
timedir = '20191016_Lotus_03_lineartrack.time/';
basenmrec='20191016_Lotus_03';
ptp_ctime_filename=[data_dir timedir basenmrec '_lineartrack.continuoustime.dat'];
ptp_ctime=readTrodesExtractedDataFile(ptp_ctime_filename);
tt = double(ptp_ctime.fields(1).data);
st = double(ptp_ctime.fields(2).data);
% tt_diff = diff(tt);
% max(tt_diff) % should be 1 % it's 570


