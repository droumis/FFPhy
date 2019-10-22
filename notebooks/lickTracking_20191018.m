







%% load  data
loaddir = '/home/jojo/DeepLabCut-master/videos/lick_tracking/';
filename = '20191016_Lotus_03_lineartrackDeepCut_resnet50_lineartrackOct17shuffle1_1030000.csv';
DLC_tracking = csvread([loaddir filename], 3);

%% load DIO data
data_dir='/stelmo/demetris/Lotus/AJ_processed/20191016/20191016_Lotus_03_lineartrack/';
timedir = '20191016_Lotus_03_lineartrack.time/';
basenm='20191016_Lotus_';
cd(data_dir)
%% 1-Get timestamps from .rec file
clockrate=30000;
rec=['03']; 
% ptp_time_filename=[timedir, basenm, num2str(rec), '_','lineartrack.time.dat'];
% ptp_time=readTrodesExtractedDataFile(ptp_time_filename);
ptp_ctime_filename=[timedir, basenm, num2str(rec), '_','lineartrack.continuoustime.dat'];
ptp_ctime=readTrodesExtractedDataFile(ptp_ctime_filename);
tt = double(ptp_ctime.fields(1).data);
st = double(ptp_ctime.fields(2).data);
tt_diff = diff(tt);
max(tt_diff) %good if this value is 1 % it's 570

%% 2-Get camera timestamps from CamHWlog

camHWlog_filename=[basenm,num2str(rec),'_','lineartrack','.1.videoTimeStamps.cameraHWSync'];
camHWlog = readTrodesExtractedDataFile(camHWlog_filename);

%camHWlog.fields.name;
t_camHW_ns = double(camHWlog.fields(3).data);
n_camHW_framecount = camHWlog.fields(2).data;
n_dropped_frames = sum(diff(n_camHW_framecount)-1);

t_camHW_s = t_camHW_ns/1e9;
t_camHW_ms = t_camHW_ns/1e6;

% camHWlog_n_records = size(t_camHW_ns,1);
fprintf('\n');
fprintf('# of CamHW timestamps: %d \n', length(t_camHW_ms));
fprintf('# of dropped frames: %d \n', n_dropped_frames);

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
% not yay. why is t_camHW_ms less than bls(2)? e15 vs e18
% i changed t_camHW_ms to t_camHW_ns bc the trodes data is in ns.. i think
% that's the fix
cam_tt_fit = (t_camHW_ns -bls(2)) ./ bls(1);
cam_rt_fit=cam_tt_fit./clockrate;
median(diff(cam_rt_fit(1:1000)))%= 8 ms
% ok i confirmed with AJ that things are probably in ns now and that the
% way i understand things is correct.. 


% %% 4-Load Deeplabcut results 
% fprintf('\n')
% dlc_filename=[basenm, num2str(rec), '_','lineartrackDeepCut_resnet50_lineartrackshuffle1_1030000.csv'];
% dlc_results=readtable(dlc_filename);
% 
% %modify the table headers
% dlc_results.Properties.VariableNames={'camfram' 'nose_x' 'nose_y' 'nose_likelihood' 'eye_x' 'eye_y' 'eye_likelihood' 'tongue_x' 'tongue_y' 'tongue_likelihood' 'jaw_x' 'jaw_y' 'jaw_likelihood' 'ear_x' 'ear_y' 'ear_likelihood'};
% 
% %load nose
% dlc_nose_x=str2double((dlc_results.nose_x(3:end)));
% dlc_nose_y=str2double(dlc_results.nose_y(3:end));
% dlc_nose_fit=str2double(dlc_results.nose_likelihood(3:end));
% 
% %load eye
% dlc_eye_x = str2double(dlc_results.eye_x(3:end));
% dlc_eye_Y = str2double(dlc_results.eye_y(3:end));
% dlc_eye_fit= str2double(dlc_results.eye_likelihood(3:end));
% 
% %jaw
% dlc_jaw_x= str2double(dlc_results.jaw_x(3:end));
% dlc_jaw_y= str2double(dlc_results.jaw_y(3:end));
% dlc_jaw_fit= str2double(dlc_results.jaw_likelihood(3:end));
% 
% %tongue
% dlc_tongue_x= str2double(dlc_results.tongue_x(3:end));
% dlc_tongue_y=str2double(dlc_results.tongue_y(3:end));
% dlc_tongue_fit=str2double(dlc_results.tongue_likelihood(3:end));
% 
% %ear
% dlc_ear_x= str2double(dlc_results.ear_x(3:end));
% dlc_ear_y= str2double(dlc_results.ear_y(3:end));
% dlc_ear_fit= str2double(dlc_results.ear_likelihood(3:end));
% 
% dlc_n_records = size(dlc_ear_x,1);
% fprintf('# of DeepLabCut timestamps: %d \n', dlc_n_records)
% 
% fprintf('Time estimate = %g s \n', dlc_n_records/est_framerate) %=1005.1
% %% AJ plot 
% figure
% plot([dlc_nose_x dlc_eye_x dlc_jaw_x dlc_tongue_x dlc_ear_x], '.')
% legend('nose', 'eye', 'jaw', 'tongue', 'ear')
% clear VarName*;
%% Demetris..
diodir = '/stelmo/demetris/Lotus/AJ_processed/20191016/20191016_Lotus_03_lineartrack/20191016_Lotus_03_lineartrack.DIO/';
dio = struct;
dct = {'in', 'out'};
c = 0;
for d = 1:length(dct)
    for i = 1:32
        %     dio(i).label = sprintf('Din%d', i);
        di = readTrodesExtractedDataFile([diodir sprintf(...
            '20191016_Lotus_03_lineartrack.dio_D%s%d.dat',dct{d}, i)]);
        di.times = di.fields(1).data;
        di.numtimes = numel(di.times);
        di.state = di.fields(2).data;
        try
            c = c+1;
            dio(c) = di;
        catch
            dio = di;
        end
    end
end

%% gather the dio
numtimes = cell2mat({dio.numtimes}');
validDIO = find(numtimes > 1);

lickInDIO = 8; % 8
lickTimes = double(dio(lickInDIO).times(dio(lickInDIO).state==1));
lickTimes = double(lickTimes)./clockrate;

rewOutDIO = 37; % 37
rewTimes = double(dio(rewOutDIO).times(dio(rewOutDIO).state==1));
rewTimes = double(rewTimes)./clockrate;

%% get lick bout intervals
% bounds where at least 5 seconds ILI
% only include lick bout intervals that are at least 1 second long
maxAllowedBoutGap = 1; % seconds
minBoutLicks = 10; % licks

boutStart = lickTimes(find(diff([diff(lickTimes) < maxAllowedBoutGap]) == 1)+1);
boutEnd = lickTimes(find(diff([diff(lickTimes) > maxAllowedBoutGap]) == 1)+1);
while boutStart(1) > boutEnd(1)
    boutEnd(1) = [];
    if boutEnd(end)<boutStart(end)
        boutStart(end) = [];
    end
end
% filter out bouts with less than boutNum licks
% try
licksInbout = logical(isExcluded(lickTimes, [boutStart boutEnd]));
% catch
%     fprintf('error defining lick bouts for %d %d\n', day, epoch)
%     continue
% end
N = histcounts(lickTimes(licksInbout), sort([boutStart; boutEnd]));
inclItv = N(1:2:end) > minBoutLicks;
boutStart = boutStart(inclItv);
boutEnd = boutEnd(inclItv);
lickBIntvl = [boutStart boutEnd];

% lickBoutBol = logical(isExcluded(st_fit, lickBIntvl));
%     line([dtimes';dtimes'], ylim)
%     pause
%     dioID = cellfun(@(x) str2double(regexp(x.original_id,'\d*','Match')), ...
%         dio(usedios(di)), 'un', 1);
inclTimes = logical(isExcluded(cam_rt_fit, lickBIntvl));
%% gather pos data
% time = [1:size(DLC_tracking,1)] * 1/125;
L = 3;
nose = DLC_tracking(:,2:4);
noseValid = all([[nose(:,L) > Pp.likethresh] inclTimes],2);

eye = DLC_tracking(:,5:7);
eyeValid = all([[eye(:,L) > Pp.likethresh] inclTimes],2);

tongue = DLC_tracking(:,8:10);
tonValid = all([[tongue(:,L) > Pp.likethresh] inclTimes],2);

jaw = DLC_tracking(:,11:13);
jawValid = all([[jaw(:,L) > Pp.likethresh] inclTimes],2);

ear = DLC_tracking(:,14:16);
earValid = all([[ear(:,L) > Pp.likethresh] inclTimes],2);

% xpos = DLC_tracking(:,2:3:end);
% ypos = DLC_tracking(:,3:3:end);
% lkl = DLC_tracking(:,4:3:end);

%% plot each X Y lkl
figname = 'licktracking';
pausefigs = 1;
savefigs = 0;

Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
            'savefigs', savefigs);
        
sf1 = subaxis(5,1,1,Pp.posparams{:});
scatter(cam_rt_fit(noseValid)', nose(noseValid,1)-mean(nose(noseValid,1)), Pp.Msz, 'k.')
% hold on
% scatter(cam_rt_fit(noseValid)', nose(noseValid,2)-mean(nose(noseValid,2)), Pp.Msz, 'm.')
% hold off
axis tight
ylabel('nose')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1.5, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1);

sf2 = subaxis(5,1,2,Pp.posparams{:});
scatter(cam_rt_fit(eyeValid)', eye(eyeValid,1)-mean(eye(eyeValid,1)), Pp.Msz, 'k.')
% hold on
% scatter(cam_rt_fit(eyeValid)', eye(eyeValid,2)-mean(eye(eyeValid,2)), Pp.Msz, 'm.')
% hold off
axis tight
ylabel('eye')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1.5, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1);

sf3 = subaxis(5,1,3,Pp.posparams{:});
% scatter(cam_rt_fit(tonValid)', tongue(tonValid,1)-mean(tongue(tonValid,1)), Pp.Msz, 'k.')
% hold on 
scatter(cam_rt_fit(tonValid)', tongue(tonValid,2)-mean(tongue(tonValid,2)), Pp.Msz, 'k.')
% hold off
axis tight
ylabel('tongue x')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1.5, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1);

sf4 = subaxis(5,1,4,Pp.posparams{:});
scatter(cam_rt_fit(jawValid)', jaw(jawValid,1)-mean(jaw(jawValid,1)), Pp.Msz, 'k.')
% hold on
% scatter(cam_rt_fit(jawValid)', jaw(jawValid,2)-mean(jaw(jawValid,2)), Pp.Msz, 'm.')
% hold off
axis tight
ylabel('jaw')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1.5, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1);

sf5 = subaxis(5,1,5,Pp.posparams{:});
scatter(cam_rt_fit(earValid)', ear(earValid,1)-mean(ear(earValid,1)), Pp.Msz, 'k.')
% hold on
% scatter(cam_rt_fit(earValid)', ear(earValid,2)-mean(ear(earValid,2)), Pp.Msz, 'm.')
% hold off
axis tight
ylabel('ear')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1.5, 'linestyle', ':');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1);

h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p=pan; p.Motion='horizontal';

allAxesInFigure = findall(gcf,'type','axes');
linkaxes(allAxesInFigure, 'x');

%% filter data, get phase of dio, polar plot, report Raleigh stat
% make cmw centered at 7 Hz. fft.
wp = getWaveParams('Lick7Hz'); 

% creating the morlet wavelet by combining the complex sine wave and the gaussian
timewv = [1:length(cam_rt_fit)]./wp.srate;
sine_wave = exp(2*1i*pi*wp.freqLog.*timewv'); %make a complex sine wave at freq
waveStd = wp.nWaveCyc/(2*pi*wp.freqLog); % std of wavelet gaussian. depends on freq, #cycles
gausWin = exp(-cam_rt_fit.^2./(2*waveStd^2)); % gaussian
wavelet = sine_wave .* gausWin; 
wavefft = fft(wavelet); %,nConv2pow); % take the fft of the wavelet pad with zeros to the next power of 2 for speed
% normalize wavelet to a maximum of 1 to ensure convolution units are same as data
waveletFFT(:,1) = (wavefft ./ max(wavefft))';

% fft data, conv with wavelet, then unfft
dataFFT = fft(nose(:,1));
dconv = dataFFT.*waveletFFT;
as = ifft(astmp);

phase = angle(as);
%%
Fs = 125;  % Sampling Frequency

Fstop1 = 4;                % First Stopband Frequency
Fpass1 = 6;                % First Passband Frequency
Fpass2 = 12;               % Second Passband Frequency
Fstop2 = 14;               % Second Stopband Frequency
Dstop1 = 0.0031622776602;  % First Stopband Attenuation
Dpass  = 0.057501127785;   % Passband Ripple
Dstop2 = 0.0031622776602;  % Second Stopband Attenuation
dens   = 20;               % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                          0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

%% aj's method
data = {jaw, tongue, eye, nose, ear};
figure
for s = 1:5
    dataFilt=filtfilt(b,1,double(data{s}(:,1)));
    HT = hilbert(dataFilt);
    amplitude = sqrt(real(HT).^2 + imag(HT).^2);
    dataPhase = angle(HT);
    
    lickCamIdx = knnsearch(cam_rt_fit, lickTimes);
    lickDataPhase = dataPhase(lickCamIdx);
    
    [Rp, z] = circ_rtest(lickDataPhase);
    fprintf('Rp:%.04f z:%.04f\n', Rp, z)
    % PLOT
    polarhistogram(lickDataPhase, 32, 'Normalization', 'pdf', 'edgealpha', .1); %,...
    %     'facecolor', [.5 .5 .5]);
    hold on
end
legend({'jaw', 'tongue', 'eye', 'nose','ear' })
hold off
%% 
figname = 'licktracking';
pausefigs = 1;
savefigs = 0;

Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
            'savefigs', savefigs);

sf1 = subaxis(5,1,1,Pp.posparams{:});
scatter(cam_rt_fit(noseValid)', nose(noseValid,1)-mean(nose(noseValid,1)), Pp.Msz, 'k.')
hold on;
plot(cam_rt_fit(noseValid)', dataFilt(noseValid))
axis tight
ylabel('nose')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1.5, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1);

xlim([763 768]);
%%
sf2 = subaxis(5,1,2,Pp.posparams{:});
scatter(cam_rt_fit(eyeValid)', eye(eyeValid,1)-mean(eye(eyeValid,1)), Pp.Msz, 'k.')
% hold on
% scatter(cam_rt_fit(eyeValid)', eye(eyeValid,2)-mean(eye(eyeValid,2)), Pp.Msz, 'm.')
% hold off
axis tight
ylabel('eye')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1.5, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1);

sf3 = subaxis(5,1,3,Pp.posparams{:});
% scatter(cam_rt_fit(tonValid)', tongue(tonValid,1)-mean(tongue(tonValid,1)), Pp.Msz, 'k.')
% hold on 
scatter(cam_rt_fit(tonValid)', tongue(tonValid,2)-mean(tongue(tonValid,2)), Pp.Msz, 'k.')
% hold off
axis tight
ylabel('tongue x')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1.5, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1);

sf4 = subaxis(5,1,4,Pp.posparams{:});
scatter(cam_rt_fit(jawValid)', jaw(jawValid,1)-mean(jaw(jawValid,1)), Pp.Msz, 'k.')
% hold on
% scatter(cam_rt_fit(jawValid)', jaw(jawValid,2)-mean(jaw(jawValid,2)), Pp.Msz, 'm.')
% hold off
axis tight
ylabel('jaw')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1.5, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1);

sf5 = subaxis(5,1,5,Pp.posparams{:});
scatter(cam_rt_fit(earValid)', ear(earValid,1)-mean(ear(earValid,1)), Pp.Msz, 'k.')
% hold on
% scatter(cam_rt_fit(earValid)', ear(earValid,2)-mean(ear(earValid,2)), Pp.Msz, 'm.')
% hold off
axis tight
ylabel('ear')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1.5, 'linestyle', ':');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1);

h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p=pan; p.Motion='horizontal';

allAxesInFigure = findall(gcf,'type','axes');
linkaxes(allAxesInFigure, 'x');

%%
% %% or using hilbert
% Fs = 125;                                               % Sampling Frequency
% Fn = Fs/2;                                              % Nyquist Frequency
% Wp = [6  10]/Fn;                                       % Theta Passband
% Ws = Wp .* [0.9  1.15];                                 % Bandstop Frequencies
% Rp = 10;                                                % Passband Ripple
% Rs = 40;                                                % Stopband Ripple
% [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                         % Determine Optimal Order
% [b,a] = cheby2(n,Rs,Ws);                                % Transfer Function Coefficients
% [sos_theta,g_theta] = tf2sos(b,a);                      % Second-Order-Section For Stability
% figure(1)
% freqz(sos_theta, 4096, Fs);                             % Filter Bode Plot
%%
% y = filtfilt(d,x);
% filtnose = filtfilt(nose(:,1));
% HT = hilbert(filtnose);
% phase = angle(HT);






