

%{
XP tracking phase consistency

%}
clockrate=30000;
likethresh = .95;
animal = 'Lotus';

andef = animaldef(animal);
date = '20191016';
epoch = 3;

Label = 'headPhaseXP';
plotfigs = 1;
showfigs = 0;
pausefigs = 0;
savefigs = 1;
savefigas = {'png', 'pdf'};

plot_polarHist = 1;


% FIRPM Parks-McClellan optimal equiripple FIR filter design
Fstop1 = 4;                % First Stopband Frequency
Fpass1 = 7;                % First Passband Frequency
Fpass2 = 9;               % Second Passband Frequency
Fstop2 = 12;               % Second Stopband Frequency
Dstop1 = 0.0031622776602;  % First Stopband Attenuation
Dpass  = 0.057501127785;   % Passband Ripple
Dstop2 = 0.0031622776602;  % Second Stopband Attenuation
dens   = 20;               % Density Factor


%% make continuoustime.dat
pp_idx = find(cell2mat(cellfun(@(x) contains(x,'preprocessing'), ...
    andef, 'un', 0)));
preprocdaydir =  sprintf('%s%s/',andef{pp_idx}, date);

timefile=dir([preprocdaydir '*.time/*continuoustime.dat']);
% if isempty(timefile.name)
%     % ensures 5 minute gaps between epochs
%     ptp_adjust_timestamps(preprocdaydir)
%     timefile=dir([preprocdaydir '*.time/*continuoustime.dat']);
% end

%% Get trodestime and systime
ptp_ctime_filename = sprintf('%s/%s',timefile.folder, timefile.name);
ptp_ctime=readTrodesExtractedDataFile(ptp_ctime_filename);
tt_idx = find(cell2mat(cellfun(@(x) strcmp(x,'trodestime'), ...
    {ptp_ctime.fields.name}, 'un', 0)));
st_idx = find(cell2mat(cellfun(@(x) strcmp(x,'systime'), ...
    {ptp_ctime.fields.name}, 'un', 0)));
tt = double(ptp_ctime.fields(tt_idx).data);
st = double(ptp_ctime.fields(st_idx).data);

%% Get camera timestamps
raw_idx = find(cell2mat(cellfun(@(x) contains(x,'raw'), ...
    andef, 'un', 0)));
camHWlog_file = dir(sprintf('%s%s/*.cameraHWSync',andef{raw_idx}, date));
camHWlog_filename = sprintf('%s/%s',camHWlog_file.folder,camHWlog_file.name);
camHWlog = readTrodesExtractedDataFile(camHWlog_filename);

t_camHW_idx = find(cell2mat(cellfun(@(x) contains(x,'HWTimestamp'), ...
    {camHWlog.fields.name}, 'un', 0)));
n_camHW_idx = find(cell2mat(cellfun(@(x) contains(x,'HWframeCount'), ...
    {camHWlog.fields.name}, 'un', 0)));

t_camHW_ns = double(camHWlog.fields(t_camHW_idx).data);
n_camHW_framecount = camHWlog.fields(n_camHW_idx).data;
n_dropped_frames = sum(diff(n_camHW_framecount)-1);

t_camHW_s = t_camHW_ns/1e9;

fprintf('# of CamHW timestamps: %d \n', length(t_camHW_ns));
fprintf('# of dropped frames: %d \n', n_dropped_frames);
est_framerate = median(1./diff(t_camHW_s));
fprintf('Estimated frame rate from camHW: %0.3f \n', est_framerate);
fprintf('First/Last record: %0.3f %0.3f (%0.3f elapsed) \n', t_camHW_s(1), t_camHW_s(end), diff(t_camHW_s([1 end])));

%% regress systime to trodes time
XX = tt;
Y = st;
bls = polyfit(XX,Y,1);

% 'reconstructed' sysClock times of each Trodes packet
% st_fit = tt * bls(1) + bls(2);
cam_tt_fit = (t_camHW_ns -bls(2)) ./ bls(1);
cam_rt_fit=cam_tt_fit./clockrate;
% median(diff(cam_rt_fit(1:1000))) %= 8 ms


%% get DLC dir
% before running this, ensure that the desired cld .csv is in an/dlc/date/
DAE = sprintf('%s_%s_%02d%s', date, animal, epoch);
dlc_idx = find(cell2mat(cellfun(@(x) contains(x,'dlc'), ...
    andef, 'un', 0)));
file = dir([andef{dlc_idx} '/' date '/' DAE '*.csv']);
fprintf('\n importing::: %s \n',file.name)

%% load DLC result CSV
f = sprintf('%s/%s', file.folder, file.name);
DLC_tracking = csvread(f, 3);
opts = detectImportOptions(f,'NumHeaderLines',1); % number of header lines which are to be ignored
DLC_fields = opts.VariableNames;

%% Filter based on likelihood
% inclTimes = [1:length(cam_rt_fit)]';
eyex_idx = find(cell2mat(cellfun(@(x) strcmp(lower(x),'eye'), ...
    DLC_fields, 'un', 0)));
eyey_idx = find(cell2mat(cellfun(@(x) strcmp(lower(x),'eye_1'), ...
    DLC_fields, 'un', 0)));
eyeL_idx = find(cell2mat(cellfun(@(x) strcmp(lower(x),'eye_2'), ...
    DLC_fields, 'un', 0)));

eye = DLC_tracking(:,[eyex_idx eyey_idx eyeL_idx]);

L = 3; %likelihood
% eyeValid = all([[eye(:,L) > likethresh] inclTimes],2);
eyeValid = all(eye(:,L) > likethresh,2);
fprintf('eye valid for %.02f pct of time \n', sum(eyeValid)/length(eyeValid))

%% % get 4-12 Hz phase of eye tracking
Fs = est_framerate;  % Sampling Frequency
% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
    0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});

eyeFilt_x = filtfilt(b,1,double(eye(:,1)));
eyePhase_x = angle(hilbert(eyeFilt_x));

eyeFilt_y = filtfilt(b,1,double(eye(:,2)));
eyePhase_y = angle(hilbert(eyeFilt_y));

%% Get XP and reward DIO..

diodir = dir([preprocdaydir '*.DIO']);
dio_path = [preprocdaydir diodir.name '/'];
% "adjusted" timestamps are created by ptp_adjust_timestamps
% DONT USE ADJUSTED TIMESTAMPS!!
% adjdio = dir([preprocdaydir diodir.name '/*.adj.dat']);
diodat = dir([preprocdaydir diodir.name '/*.dat']);

dio = struct;
% dct = {'in', 'out'};
c = 0;

for i = 1:length(diodat)
    if strfind(diodat(i).name, 'adj')
        continue
    end
    di = readTrodesExtractedDataFile([dio_path diodat(i).name]);
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
% gather the dio
% numtimes = cell2mat({dio.numtimes}');
% validDIO = find(numtimes > 1);

% find lick input DIO's by the mean diff
dfrex = [];
for d = 1:length(dio)
    v = (double(dio(d).times(dio(d).state==1)) / clockrate);
    dfrex(d) = median(1./diff(v));
end
n = find(all([(dfrex >= 2)' (dfrex <= 20)'],2));

% find DIO with Video Tracking by the min diff from tracked vid time
eyeValid = all(eye(:,L) >= 1,2);
trackedTimes = cam_rt_fit(eyeValid);
fst = [];
for d = 1:length(n)
    % median diff between nearest tracked timepoint and XP
    v = double(dio(n(d)).times(dio(n(d)).state==1)) ./ clockrate;
    for t = 1:length(v)
        dnear(t) = min(abs(v(t) - trackedTimes));
    end
    fst(d) = median(dnear);
end
[~,i] = min(fst);
trackedXP = n(i);
fprintf('determined DIO # %d for XP with tracking \n', trackedXP)
XPtimes = double(dio(n(i)).times(dio(n(i)).state==1))./clockrate;
%% Feature phase at port crossing times
% the DLC frame times are the Camera Frame times, so use the regressed
% camera times along with the XPtimes to extract XP phase
XPCamIdx = knnsearch(cam_rt_fit, XPtimes);

XP_eyeXphase = eyePhase_x(XPCamIdx);
XP_eyeYphase = eyePhase_y(XPCamIdx);
%% PLOT polar histogram

% polarhistogram(jawDInPhase, 32, 'Normalization', 'pdf', 'edgealpha', .1, 'FaceColor', [.8, .4, .05]);
% hold on
% polarhistogram(tonDInPhase, 32, 'Normalization', 'pdf', 'edgealpha', .1);


% [Rp, z] = circ_rtest(XP_eyeXphase);
% fprintf('EYE X==== Rp:%.04f z:%.04f\n', Rp, z)
%
% [Rp, z] = circ_rtest(XP_eyeYphase);
% fprintf('EYE Y==== Rp:%.04f z:%.04f\n', Rp, z)

% title('Feature Phase of Dinput')
% if savefigs

%     save_figure('demetris', stit, ...
%         'savefigas', savefigas, 'subdir', figname)
% end
if plotfigs
    if plot_polarHist
        figname = sprintf('%s-pAn',Label);
        Pp=load_plotting_params({'defaults',figname});
        polarhistogram(XP_eyeXphase, 20, 'Normalization', 'pdf', ...
            'edgealpha', .1, 'FaceColor', [.1, .9, 1]);
        hold on
        polarhistogram(XP_eyeYphase, 20, 'Normalization', 'pdf', ...
            'edgealpha', .1, 'FaceColor', [.9, .1, 1]);
        legend({'eyeX' 'eyeY'})
        hold off
        title('eyeXY')
        % super
        stit = sprintf('%s %s', figname, animal);
        setSuperAxTitle(stit);
        hold off
        if pausefigs
            pause
        end
        if savefigs
            strsave = save_figure('demetris', stit, 'savefigas', savefigas, ...
                'subdir', figname);
        end
    end
    
end