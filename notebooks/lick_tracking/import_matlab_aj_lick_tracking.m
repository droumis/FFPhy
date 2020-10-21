%% import_matlab_aj_lick_tracking:
%this is a quick and dirty way to load single epoch .rec and DLC data 

clear all;
close all;
clc;

%% !!! only need to do once !!!
%%go to the folder with all the raw files. using trodes to matlab functions to extract all the data 

filename = '20191016_Lotus_03_lineartrack'
% extractTimeBinaryFile(filename);
% extractSpikeBinaryFiles(filename);
% extractLFPBinaryFiles(filename);
% extractDioBinaryFiles(filename);
%all files extracted succesfully in the same folder

%%change directory
data_dir='/opt/data15/animals/Lotus/20191016'
cd(data_dir)

clockrate=30000;
ref_tetrode_num=(13)

basenm='20191016_Lotus_';
rec=['03']; 
well1_Dio=[18]; % add the correct DIO 
well2_Dio=[18]; % add the correct DIO 

%% 1-Get timestamps from .rec file
ptp_time_filename=[basenm, num2str(rec), '_','lineartrack.time.dat'];
ptp_time=readTrodesExtractedDataFile(ptp_time_filename);

ptp_ctime_filename=[basenm, num2str(rec), '_','lineartrack.continuoustime.dat']
ptp_ctime=readTrodesExtractedDataFile(ptp_ctime_filename);

% Trodes time (30kHz, ticks elapsed since MCU was restarted)
tt = double(ptp_ctime.fields(1).data);
% System time (PC clock time when each sample arrived from the MCU)
% CHANGed! from milliseconds to nanoseconds in 1.8.2
st = double(ptp_ctime.fields(2).data);

% check for dropped packets
tt_diff = diff(tt);
max(tt_diff); %good if this value is 1

%% 2-Get camera timestamps from CamHWlog
camHWlog_filename=[basenm, num2str(rec), '_','lineartrack','.1.videoTimeStamps.cameraHWSync'];
camHWlog = readTrodesExtractedDataFile(camHWlog_filename);

%camHWlog.fields.name;
t_camHW_ns = double(camHWlog.fields(3).data);
n_camHW_framecount = camHWlog.fields(2).data;
n_dropped_frames = sum(diff(n_camHW_framecount)-1);

t_camHW_s = t_camHW_ns/1e9;
t_camHW_ms = t_camHW_ns/1e6;

camHWlog_n_records = size(t_camHW_ns,1);
fprintf('\n');
fprintf('# of CamHW timestamps: %d \n', length(t_camHW_ms))
fprintf('# of dropped frames: %d \n', n_dropped_frames)

est_framerate = median(1./diff(t_camHW_ns./1e9));
fprintf('Estimated frame rate from camHW: %0.3f \n', est_framerate);
fprintf('First/Last record: %0.3f %0.3f (%0.3f elapsed) \n', t_camHW_s(1), t_camHW_s(end), diff(t_camHW_s([1 end])));

%% 3-Get corresponding Trodes times for camHWlog times
% A simple lookup (below), adds jitter from the MCU->PC communication
% so we'll perform a fit of those timestamps instead
% tt_index = lookup(t_camHW_ms, st);
% cam_tt = tt(tt_index);

% There is some jitter in the arrival times of packets from the MCU (as
% reflected in the sysclock records in the .rec file. Let's assume that
% Trodes clock is actually regular, and recover the clock correspondence
% with a regression
XX = tt; %camHWclock; 
Y = st; %camlogPC; 
bls = polyfit(XX,Y,1);

% 'reconstructed' sysClock times of each Trodes packet
st_fit = tt * bls(1) + bls(2);

% Trodes times corresponding to each camera frame. Yay!!
cam_tt_fit = (t_camHW_ms -bls(2)) ./ bls(1);
cam_rt_fit=cam_tt_fit./clockrate;
% median(diff(cam_rt_fit))=8.001242280 ms

%% 4-Load Deeplabcut results 
fprintf('\n')
dlc_filename=[basenm, num2str(rec), '_','lineartrackDeepCut_resnet50_lineartrackshuffle1_1030000.csv'];
dlc_results=readtable(dlc_filename);

%modify the table headers
dlc_results.Properties.VariableNames={'camfram' 'nose_x' 'nose_y' 'nose_likelihood' 'eye_x' 'eye_y' 'eye_likelihood' 'tongue_x' 'tongue_y' 'tongue_likelihood' 'jaw_x' 'jaw_y' 'jaw_likelihood' 'ear_x' 'ear_y' 'ear_likelihood'};

%load nose
dlc_nose_x=str2double((dlc_results.nose_x(3:end)));
dlc_nose_y=str2double(dlc_results.nose_y(3:end));
dlc_nose_fit=str2double(dlc_results.nose_likelihood(3:end));

%load eye
dlc_eye_x = str2double(dlc_results.eye_x(3:end));
dlc_eye_Y = str2double(dlc_results.eye_y(3:end));
dlc_eye_fit= str2double(dlc_results.eye_likelihood(3:end));

%jaw
dlc_jaw_x= str2double(dlc_results.jaw_x(3:end));
dlc_jaw_y= str2double(dlc_results.jaw_y(3:end));
dlc_jaw_fit= str2double(dlc_results.jaw_likelihood(3:end));

%tongue
dlc_tongue_x= str2double(dlc_results.tongue_x(3:end));
dlc_tongue_y=str2double(dlc_results.tongue_y(3:end));
dlc_tongue_fit=str2double(dlc_results.tongue_likelihood(3:end));

%ear
dlc_ear_x= str2double(dlc_results.ear_x(3:end));
dlc_ear_y= str2double(dlc_results.ear_y(3:end));
dlc_ear_fit= str2double(dlc_results.ear_likelihood(3:end));

dlc_n_records = size(dlc_ear_x,1);
fprintf('# of DeepLabCut timestamps: %d \n', dlc_n_records)

fprintf('Time estimate = %g s \n', dlc_n_records/est_framerate) %=1005.1

plot([dlc_nose_x dlc_eye_x dlc_jaw_x dlc_tongue_x dlc_ear_x], '.')
legend('nose', 'eye', 'jaw', 'tongue', 'ear')
clear VarName*;

%% 5-Get DIOs directly from the .rec file

dout_stim_filename=[basenm,num2str(rec),'_lineartrack.dio_Dout',num2str(laserDio),'.dat'];
dout_stim_dat=readTrodesExtractedDataFile(dout_stim_filename);

% Dout_UP_16_TT are those times when the DIO was==1 
dout_stim_tt = double(dout_stim_dat.fields(1).data(dout_stim_dat.fields(2).data==1));

%dout16_rt is the real time in the recoridng file
dout_stim_rt=double(dout_stim_tt)./clockrate;

% %% 6-Load LFP data
% 
% LFP_filename=[basenm, num2str(rec), '_','lineartrack','.LFP_nt', num2str(ref_tetrode_num),'ch1.dat'];
% LFP=readTrodesExtractedDataFile(LFP_filename);
% LFP_data=LFP.fields.data;
% 
% lfp_time_filename =[basenm,num2str(rec),'_', 'lineartrack','.timestamps.dat'];
% LFP_time=readTrodesExtractedDataFile(lfp_time_filename);
% LFP_time_dat=LFP_time.fields.data;
% LFP_rt_dat=double(LFP_time_dat)./clockrate;
% hist(double(diff(LFP_time_dat)), 100);
% % all values are 20 points away
% 
% % estimation of recording time accoridng to length of LFP_time_dat
% (length(LFP_time_dat)/1500)/60 %=29.89minutes
% 
% %%plot results
% load('thetafilter.mat');
% theta_filtered_lfp=filtfilt(thetafilter.tf.num,thetafilter.tf.den,double(LFP_data));
% 
% %calculating hilbert transform phase for theta 
% HT = hilbert(theta_filtered_lfp);
% amplitude = sqrt(real(HT).^2 + imag(HT).^2);
% phase = angle(HT);
% 
% %% 7-plot the data
% close all; figure(1); hold on
% plot(LFP_rt_dat,LFP_data, '-k')
% plot([dout_stim_rt dout_stim_rt],[-1000 1000], '-m')
% %ylim([-4000 4000]);
% %xlim([655 660]);
% plot(LFP_rt_dat,theta_filtered_lfp, '-r', 'LineWidth', 0.5);
% %plot(LFP_rt_dat(:), (phase.*100)+2500, 'g','LineWidth', 1)
% legend({'raw lfp', 'stim times','theta filtered lfp'})
% xlabel('time')
% ylabel('amplitude uV')
% title('hippocampal lfp and laser stim')
% 
% %% 8-use cont utilities to plot triggered averegaes etc 
% cdat_lfp=imcont('data', LFP_data,'timestamp', LFP_rt_dat);
% % power spectrum using welch method 
% [p1,f1]=contpsd(cdat_lfp, 'window_t', 1);
% figure(2)
% plot(f1,p1);hold on;
% xlim([0 20])
% 
% stim_period=[2570 2582;2630 2642;2691 2703;2751 2763;2810 2822;2870 2882;2931 2943;2991 3004;3051 3063;3111 3123]; 
% %656 659;666 669;676 680;685 690; 705 708; 752 755;775 785; 945 955;
% [p2,f2]=contpsd(cdat_lfp,'window_t',1,'segs',stim_period);
% plot(f2,p2);
% 
% %%an easy thing to evaluate would be the psd inside and outside stim
% %%periods 
% 
% %% caculate the phase of theta duirng stims 
% %calculate the phase of theta during swing
% 
% phase_stim = nan(length(dout_stim_rt),1);
% for s = 1:length(dout_stim_rt)
% ind = find(LFP_rt_dat>= dout_stim_rt(s));
% phase_stim(s) = phase(ind(1));
% end
% 
% circ_mean(phase_stim)
% circ_r(phase_stim)
% circ_rtest(phase_stim)
% 
% %%should only collect first spike in a burst 
% %%should calculate the start and end times of the different stim
% %%contogneties 
% %%should take into account the start and end of each 