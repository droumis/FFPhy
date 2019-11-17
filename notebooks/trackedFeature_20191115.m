%{


plot the 2d position of the tracked features over time



                            \\\// 
                           -(o o)- 
========================oOO==(_)==OOo=======================
%}



clockrate=30000;
likethresh = .999;
animal = 'Lotus';

%% load data
data_dir='/stelmo/demetris/lick_tracking/Lotus/AJ_processed/20191016/20191016_Lotus_03_lineartrack/';
filename = '20191016_Lotus_03_lineartrackDeepCut_resnet50_lineartrackOct17shuffle1_1030000.csv';
DLC_tracking = csvread([data_dir filename], 3);
%% Get timestamps from .rec file
timedir = '20191016_Lotus_03_lineartrack.time/';
basenmrec='20191016_Lotus_03';
ptp_ctime_filename=[data_dir timedir basenmrec '_lineartrack.continuoustime.dat'];

ptp_ctime=readTrodesExtractedDataFile(ptp_ctime_filename);
tt = double(ptp_ctime.fields(1).data);
st = double(ptp_ctime.fields(2).data);
%% Get camera timestamps from videoTimeStamps.cameraHWSync
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

%% Get lick and reward DIO.. How did i make the dio? trodes extracter?
diodir = [data_dir '20191016_Lotus_03_lineartrack.DIO/'];
dio = struct;
dct = {'in', 'out'};
c = 0;
for d = 1:length(dct)
    for i = 1:32
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
showfigs = 1;
savefigs = 0;

Pp=load_plotting_params({'defaults',figname});
init_plot(showfigs, 'position', Pp.position);
sf1 = subaxis(1,2,1,Pp.posparams{:});
% plot(cam_rt_fit(noseValid)', noseFilt(noseValid)-mean(noseFilt(noseValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
% plot(cam_rt_fit(noseValid)', nosePhase(noseValid)-mean(nosePhase(noseValid)), ...
%     'color', [.5 .5 .5 .5], 'linewidth', 1)
% hold on
%%
s = 697; % s start
e = 701; % s end
sIdx = knnsearch(cam_rt_fit, s);
eIdx = knnsearch(cam_rt_fit, e);
plot(nose(sIdx:eIdx,1), nose(sIdx:eIdx,2))
%% 
s = 700; % s start
e = 701; % s end
sIdx = knnsearch(cam_rt_fit, s);
eIdx = knnsearch(cam_rt_fit, e);

% scatter(smoothdata(tongue(sIdx:eIdx,1),'loess', 10), smoothdata(tongue(sIdx:eIdx,2),'loess', 10), 20,...
%     tonPhase(sIdx:eIdx,1), 'filled')

% scatter(smoothdata(eye(sIdx:eIdx,1),'loess', 10), smoothdata(eye(sIdx:eIdx,2),'loess', 10), 20,...
%     eyePhase(sIdx:eIdx,1), 'filled')

% scatter(smoothdata(nose(sIdx:eIdx,1),'loess', 10), smoothdata(nose(sIdx:eIdx,2),'loess', 10), 20,...
%     nosePhase(sIdx:eIdx,1), 'filled')
% hold on;
% scatter(smoothdata(jaw(sIdx:eIdx,1),'loess', 10), smoothdata(jaw(sIdx:eIdx,2),'loess', 10), 20,...
%     jawPhase(sIdx:eIdx,1), 'filled')

% scatter(smoothdata(ear(sIdx:eIdx,1),'loess', 10), smoothdata(ear(sIdx:eIdx,2),'loess', 10), 20,...
%     earPhase(sIdx:eIdx,1), 'filled')

ax = gca;
ax.YDir = 'reverse';
hold off

% scatter(cam_rt_fit(noseValid)', nose(noseValid,1)-mean(nose(noseValid,1)), Pp.Msz, ...
%     '.', 'MarkerEdgeColor', [0 0 1])
% axis tight
% ylabel('nose')
% line([lickTimes';lickTimes'], ylim, 'color', Pp.lickClr, 'linewidth', 1, 'linestyle', '-');
% line([rewTimes';rewTimes'], ylim, 'color', Pp.rewClr, 'linewidth', 1.5, 'linestyle', '--');
% xticks([])
% ylim([-20 20])
% 

%% Corr matrix of the tracked features

% cohmatrixc : chronux, multitaper 
% wcoherence : matlab, waveleti






















