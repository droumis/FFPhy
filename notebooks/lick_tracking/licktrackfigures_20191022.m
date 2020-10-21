
%{
initial work in 
/home/droumis/Src/Matlab/filterframework_dr/notebooks/lickTracking_20191018.m
%}

clockrate=30000;
likethresh = .999;
animal = 'Lotus';
%% load data
% /stelmo/demetris/lick_tracking/DLC/projects/lickTrack_20191017
data_dir='/stelmo/demetris/Lotus/AJ_processed/20191016/20191016_Lotus_03_lineartrack/';
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

%%
% Get camera timestamps
camHWlog_filename = [data_dir basenmrec '_lineartrack.1.videoTimeStamps.cameraHWSync'];
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

%% Get lick and reward DIO..
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
% gather the dio
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

%% gather pos tracking data

L = 3;
nose = DLC_tracking(:,2:4);
noseValid = all([[nose(:,L) > likethresh] inclTimes],2);

eye = DLC_tracking(:,5:7);
eyeValid = all([[eye(:,L) > likethresh] inclTimes],2);

tongue = DLC_tracking(:,8:10);
tonValid = all([[tongue(:,L) > likethresh] inclTimes],2);

jaw = DLC_tracking(:,11:13);
jawValid = all([[jaw(:,L) > likethresh] inclTimes],2);

ear = DLC_tracking(:,14:16);
earValid = all([[ear(:,L) > likethresh] inclTimes],2);

%% FIRPM Parks-McClellan optimal equiripple FIR filter design
Fs = 125;  % Sampling Frequency

Fstop1 = 4;                % First Stopband Frequency
Fpass1 = 7;                % First Passband Frequency
Fpass2 = 9;               % Second Passband Frequency
Fstop2 = 12;               % Second Stopband Frequency
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

jawFilt = filtfilt(b,1,double(jaw(:,1)));
tonFilt = filtfilt(b,1,double(tongue(:,1)));
eyeFilt = filtfilt(b,1,double(eye(:,1)));
noseFilt = filtfilt(b,1,double(nose(:,1)));
earFilt = filtfilt(b,1,double(ear(:,1)));

jawPhase = angle(hilbert(jawFilt));
tonPhase = angle(hilbert(tonFilt));
eyePhase = angle(hilbert(eyeFilt));
nosePhase = angle(hilbert(noseFilt));
earPhase = angle(hilbert(earFilt));
%% Feature phase at DIn times
lickCamIdx = knnsearch(cam_rt_fit, lickTimes);
jawDInPhase = jawPhase(lickCamIdx);
tonDInPhase = tonPhase(lickCamIdx);
eyeDInPhase = eyePhase(lickCamIdx);
noseDInPhase = nosePhase(lickCamIdx);
earDInPhase = earPhase(lickCamIdx);

%% plot timeseries tracking of lickbouts with lick DIOs 
figname = 'FeatureTracking';
savefigs = 0;

Pp=load_plotting_params({'defaults',figname}, 'savefigs', savefigs);
        
sf1 = subaxis(5,1,1,Pp.posparams{:});
% plot(cam_rt_fit(noseValid)', noseFilt(noseValid)-mean(noseFilt(noseValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
% plot(cam_rt_fit(noseValid)', nosePhase(noseValid)-mean(nosePhase(noseValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
scatter(cam_rt_fit(noseValid)', nose(noseValid,1)-mean(nose(noseValid,1)), Pp.Msz, ...
    '.', 'MarkerEdgeColor', [0 0 1])
axis tight
ylabel('nose')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1.5, 'linestyle', '--');
xticks([])
ylim([-20 20])

sf2 = subaxis(5,1,2,Pp.posparams{:});
% plot(cam_rt_fit(eyeValid)', eyeFilt(eyeValid)-mean(eyeFilt(jawValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
% plot(cam_rt_fit(eyeValid)', eyePhase(eyeValid)-mean(eyePhase(jawValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
scatter(cam_rt_fit(eyeValid)', eye(eyeValid,1)-mean(eye(eyeValid,1)), Pp.Msz, '.', ...
    'MarkerEdgeColor', [.1, .9, 1])
axis tight
ylabel('eye')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1.5, 'linestyle', '--');
xticks([])
ylim([-20 20])

sf3 = subaxis(5,1,3,Pp.posparams{:});
% plot(cam_rt_fit(earValid)', earFilt(earValid)-mean(earFilt(earValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
% plot(cam_rt_fit(earValid)', earPhase(earValid)-mean(earPhase(earValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
scatter(cam_rt_fit(earValid)', ear(earValid,1)-mean(ear(earValid,1)), Pp.Msz, '.', ...
    'MarkerEdgeColor', [1, 0, 0])
axis tight
ylabel('ear')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1.5, 'linestyle', '--');
ylim([-20 20])
xticks([])

sf4 = subaxis(5,1,4,Pp.posparams{:});
% plot(cam_rt_fit(jawValid)', jawFilt(jawValid)-mean(jawFilt(jawValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
% plot(cam_rt_fit(jawValid)', jawPhase(jawValid)-mean(jawPhase(jawValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
scatter(cam_rt_fit(jawValid)', jaw(jawValid,1)-mean(jaw(jawValid,1)), Pp.Msz, '.', ...
    'MarkerEdgeColor', [.8, .4, .05]); 
axis tight
ylabel('jaw')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1.5, 'linestyle', '--');
xticks([])
ylim([-20 20])

sf5 = subaxis(5,1,5,Pp.posparams{:});
% plot(cam_rt_fit(tonValid)', tonFilt(tonValid)-mean(tonFilt(tonValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
% plot(cam_rt_fit(tonValid)', tonPhase(tonValid)-mean(tonPhase(tonValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
scatter(cam_rt_fit(tonValid)', tongue(tonValid,2)-mean(tongue(tonValid,2)), Pp.Msz, '.',...
    'MarkerEdgeColor', [0 .6 0])

axis tight
ylabel('tongue x')
line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1, 'linestyle', '-');
line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1.5, 'linestyle', '--');
xlabel('time s')

h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p=pan; p.Motion='horizontal';
allAxesInFigure = findall(gcf,'type','axes');
linkaxes(allAxesInFigure, 'x');
xlim([696 700.2])

setSuperAxTitle('Feature Tracking with DIO');

if savefigs
    save_figure('/stelmo/demetris/analysis/figures/DLCtracking/', [figname animal])
end

%% PLOT polar histogram
figname = 'DInPhase';
savefigs = 1;

Pp=load_plotting_params({'defaults',figname}, 'savefigs', savefigs);

polarhistogram(jawDInPhase, 32, 'Normalization', 'pdf', 'edgealpha', .1, 'FaceColor', [.8, .4, .05]); 
hold on
% polarhistogram(tonDInPhase, 32, 'Normalization', 'pdf', 'edgealpha', .1);
polarhistogram(eyeDInPhase, 32, 'Normalization', 'pdf', 'edgealpha', .1, 'FaceColor', [.1, .9, 1]);
% polarhistogram(noseDInPhase, 32, 'Normalization', 'pdf', 'edgealpha', .1);
% polarhistogram(earDInPhase, 32, 'Normalization', 'pdf', 'edgealpha', .1);
legend({'jaw', 'eye'})
hold off

% [jawRp, jawZ] = circ_rtest(jawLickPhase);
% fprintf('Rp:%.04f z:%.04f\n', Rp, z)
title('Feature Phase of Dinput')
if savefigs
    save_figure('/stelmo/demetris/analysis/figures/DLCtracking/', [figname animal 'jaweye'])
end



