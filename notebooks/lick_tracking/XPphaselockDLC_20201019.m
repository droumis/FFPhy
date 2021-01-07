

%{
XP tracking phase consistency
currently made for only 1 date, epoch per animal
%}
clockrate=30000;
likethresh = .95;
subjects = struct;

animals = {'GW1'}; %{'Lotus'}; %, 'GW1', 'GW27', 'Roqui'}; %'ym2', 'ym3'
Label = 'headPhaseXP';

get_XPphaseHead = 1;

plotfigs = 0;
showfigs = 1;
pausefigs = 0;
savefigs = 1;
savefigas = {'png', 'pdf'};

plot_pct_XPbobs = 0;
plot_polarHist_pAn = 0;
plot_polarHist = 0;

Shleft = -.2;
Shright = .2;

% FIRPM Parks-McClellan optimal equiripple FIR filter design
Fstop1 = 1;                % First Stopband Frequency
Fpass1 = 2;                % First Passband Frequency
Fpass2 = 8;               % Second Passband Frequency
Fstop2 = 9;               % Second Stopband Frequency
Dstop1 = 0.0031622776602;  % First Stopband Attenuation
Dpass  = 0.057501127785;   % Passband Ripple
Dstop2 = 0.0031622776602;  % Second Stopband Attenuation
dens   = 20;               % Density Factor

subjects(1).animal = 'ym3';
subjects(1).date = '20191107';
subjects(1).epoch = 1;

subjects(2).animal = 'ym2';
subjects(2).date = '20191107';
subjects(2).epoch = 1;

subjects(3).animal = 'Lotus';
subjects(3).date = '20191016';
subjects(3).epoch = 3;

subjects(4).animal = 'GW1';
subjects(4).date = '20201007';
subjects(4).epoch = 1;

subjects(5).animal = 'GW27';
subjects(5).date = '20201007';
subjects(5).epoch = 1;

subjects(6).animal = 'Roqui';
subjects(6).date = '20191107';
subjects(6).epoch = 1;
%%
if get_XPphaseHead
    
    XP_eyeXphase = {};
    XP_earXphase = {};
    XP_noseXphase = {};
    XP_eyeXphase_Sh = {};
    XP_earXphase_Sh = {};
    XP_noseXphase_Sh = {};
    
    meanMRVmag_eye = {};
    vecang_eye = {};
    meanMRVmag_ear = {};
    vecang_ear = {};
    meanMRVmag_nose = {};
    vecang_nose = {};
    
    phasemod_eyeX_Rp = {};
    phasemod_eyeX = {};
    phasemod_eyeX_sh = {};
    
    phasemod_earX_Rp = {};
    phasemod_earX = {};
    phasemod_earX_sh = {};
    
    phasemod_noseX_Rp = {};
    phasemod_noseX = {};
    phasemod_noseX_sh = {};
    
    numXP = {};
    pctBobsXP = {};
    for a = 1:length(animals)
        anidx = find(cell2mat(cellfun(@(x) strcmp(animals(a), x), ...
            {subjects.animal}, 'un', 0)));
        animal = animals{a};
        date = subjects(anidx).date;
        epoch = subjects(anidx).epoch;
        andef = animaldef(animal);
        %% get continuoustime.dat
        pp_idx = find(cell2mat(cellfun(@(x) contains(x,'preprocessing'), ...
            andef, 'un', 0)));
        preprocdaydir =  sprintf('%s%s/',andef{pp_idx}, date);
        
        timefile=dir([preprocdaydir '*.time/*continuoustime.dat']);
        if isempty(timefile)
            ptp_adjust_timestamps(preprocdaydir)
            timefile=dir([preprocdaydir '*.time/*continuoustime.dat']);
        end
        
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
        fprintf('importing %s \n',file.name)
        
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
        
        earx_idx = find(cell2mat(cellfun(@(x) strcmp(lower(x),'ear'), ...
            DLC_fields, 'un', 0)));
        eary_idx = find(cell2mat(cellfun(@(x) strcmp(lower(x),'ear_1'), ...
            DLC_fields, 'un', 0)));
        earL_idx = find(cell2mat(cellfun(@(x) strcmp(lower(x),'ear_2'), ...
            DLC_fields, 'un', 0)));
        ear = DLC_tracking(:,[earx_idx eary_idx earL_idx]);
        
        nosex_idx = find(cell2mat(cellfun(@(x) strcmp(lower(x),'nose'), ...
            DLC_fields, 'un', 0)));
        nosey_idx = find(cell2mat(cellfun(@(x) strcmp(lower(x),'nose_1'), ...
            DLC_fields, 'un', 0)));
        noseL_idx = find(cell2mat(cellfun(@(x) strcmp(lower(x),'nose_2'), ...
            DLC_fields, 'un', 0)));
        nose = DLC_tracking(:,[nosex_idx nosey_idx noseL_idx]);
        
        L = 3; %likelihood
        eyeValid = all(eye(:,L) > likethresh,2);
        fprintf('eye valid for %.02f pct of time \n', sum(eyeValid)/length(eyeValid))
        
        earValid = all(ear(:,L) > likethresh,2);
        fprintf('ear valid for %.02f pct of time \n', sum(earValid)/length(earValid))
        
        noseValid = all(nose(:,L) > likethresh,2);
        fprintf('nose valid for %.02f pct of time \n', sum(noseValid)/length(noseValid))
        
        %% % get 4-12 Hz phase of feature tracking
        Fs = est_framerate;  % Sampling Frequency
        % Calculate the order from the parameters using FIRPMORD.
        [N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
            0], [Dstop1 Dpass Dstop2]);
        
        % Calculate the coefficients using the FIRPM function.
        b  = firpm(N, Fo, Ao, W, {dens});
        
        eyeFilt_x = filtfilt(b,1,double(eye(:,1)));
        eyePhase_x = angle(hilbert(eyeFilt_x));
        
        earFilt_x = filtfilt(b,1,double(ear(:,1)));
        earPhase_x = angle(hilbert(earFilt_x));
        
        noseFilt_x = filtfilt(b,1,double(nose(:,1)));
        nosePhase_x = angle(hilbert(noseFilt_x));
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
        n = find(all([(dfrex >= 1)' (dfrex <= 20)'],2));
        
        % find DIO with Video Tracking by the min diff from tracked vid time
        eyeValid = all(eye(:,L) >= likethresh,2);
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
        
        XP_eyeXphase{a,1} = eyePhase_x(XPCamIdx);
        XP_earXphase{a,1} = earPhase_x(XPCamIdx);
        XP_noseXphase{a,1} = nosePhase_x(XPCamIdx);
        %% shuff
        N = length(XPtimes);
        
        rjitt = Shleft + (Shright-Shleft).*rand(N,1);
        XPtimes_Sh = XPtimes + rjitt;
        XPCamIdx_Sh = knnsearch(cam_rt_fit, XPtimes_Sh);
        XP_eyeXphase_Sh{a,1} = eyePhase_x(XPCamIdx_Sh);
        
        rjitt = Shleft + (Shright-Shleft).*rand(N,1);
        XPtimes_Sh = XPtimes + rjitt;
        XPCamIdx_Sh = knnsearch(cam_rt_fit, XPtimes_Sh);
        XP_earXphase_Sh{a,1} = earPhase_x(XPCamIdx_Sh);
        
        rjitt = Shleft + (Shright-Shleft).*rand(N,1);
        XPtimes_Sh = XPtimes + rjitt;
        XPCamIdx_Sh = knnsearch(cam_rt_fit, XPtimes_Sh);
        XP_noseXphase_Sh{a,1} = nosePhase_x(XPCamIdx_Sh);
        %% MRV
        
        [Rp, z] = circ_rtest(XP_eyeXphase{a,1});
        phasemod_eyeX_Rp{a,1} = Rp;
        phasemod_eyeX{a,1} = log(z);
        phasemod_eyeX_sh{a,1} = compute_phasemod_shuffle(length(XP_eyeXphase{a,1}),...
            phasemod_eyeX{a,1}, 1000);
        
        meanvec = mean(exp(1i*XP_eyeXphase{a,1})); % get mean resultant vector
        meanMRVmag_eye{a,1} = abs(meanvec); % vector magnitude
        vecang_eye{a,1} = angle(meanvec);
        
        
        [Rp, z] = circ_rtest(XP_earXphase{a,1});
        phasemod_earX_Rp{a,1} = Rp;
        phasemod_earX{a,1} = log(z);
        phasemod_earX_sh{a,1} = compute_phasemod_shuffle(length(XP_earXphase{a,1}),...
            phasemod_earX{a,1}, 1000);
        
        meanvec = mean(exp(1i*XP_earXphase{a,1})); % get mean resultant vector
        meanMRVmag_ear{a,1} = abs(meanvec); % vector magnitude
        vecang_ear{a,1} = angle(meanvec);
        
        
        [Rp, z] = circ_rtest(XP_noseXphase{a,1});
        phasemod_noseX_Rp{a,1} = Rp;
        phasemod_noseX{a,1} = log(z);
        phasemod_noseX_sh{a,1} = compute_phasemod_shuffle(length(XP_noseXphase{a,1}),...
            phasemod_noseX{a,1}, 1000);
        
        meanvec = mean(exp(1i*XP_noseXphase{a,1})); % get mean resultant vector
        meanMRVmag_nose{a,1} = abs(meanvec); % vector magnitude
        vecang_nose{a,1} = angle(meanvec);
        
        
        %% get pct of head bobs detected
        if 1
            maxIXPIthresh = 1; % seconds
            minBoutXP = 2; % # XP
            
            % get burst intervals
            ili = diff(XPtimes);
            g = [ili < maxIXPIthresh];
            bi = vec2list(g, 1:length(g));
            % include bouts with at least minBoutXP XP
            boutIntvIncl = (diff(bi,[],2)+1) >= minBoutXP;
            bii = bi(boutIntvIncl,:);
            
            % convert XP burst bound indices into time intervals
            boutIntvStart = XPtimes(bii(:,1));
            boutIntvEnd = XPtimes(bii(:,2)+1);
            boutTimes = [boutIntvStart boutIntvEnd];
            
            % get eye tracking times within XP bursts
            tILB = isIncluded(cam_rt_fit, boutTimes);
            eyeFilt_x_cp = eyeFilt_x;
            eyeFilt_x_cp(~tILB) = 0;
           
            %% plot 
            figname = sprintf('%s-pAn-pctXPbobs',Label);
            Pp=load_plotting_params({'defaults',figname});
            ifig = init_plot(showfigs, Pp.position);
            
            % raw eye
            cla
            eye_cp = eye(:,1);
            eye_cp(~tILB) = 0;
            eye_cp(eye(:,3)<.95) = 0;
            eye_cp(ear(:,3)<.95) = 0;
            eye_cp(nose(:,3)<.95) = 0;
            medeye = median(eye_cp(eye_cp>0));
            eye_cp(eye(:,1)<(medeye-25)) = 0;
            eye_cp(eye(:,1)>(medeye+25)) = 0;
            eye_cp = eye_cp - medeye;
            eyePhase_x_cp = eyePhase_x;
            excl = eye_cp==-medeye;
            eyePhase_x_cp(excl) = 0;
            eye_cp(excl) = 0;
            plot(cam_rt_fit, eye_cp(:,1))
            hold on;
            axis tight
            ylim([-100 100])
            xlim([boutTimes(4,1) boutTimes(9,2)])
            
            % filt eye
            plot(cam_rt_fit, eyePhase_x_cp)
            
            % eye peaks
            % count eye oscillation periods within XP burst
            [pks, pklocs] = findpeaks(eyePhase_x_cp);
            scatter(cam_rt_fit(pklocs), zeros(length(pks),1), '*k')
            
            % XP
            excl2 = vec2list(excl, cam_rt_fit);
            xpinc = isIncluded(XPtimes, boutTimes);
            xpinc2 = ~isIncluded(XPtimes, excl2);
            ilbXPtimes = XPtimes(all([xpinc xpinc2], 2));
            xplocs = lookup(ilbXPtimes, cam_rt_fit);
            scatter(cam_rt_fit(xplocs), zeros(length(xplocs),1), 'om')
            hold off
            
            % get pct of inta burst eye peaks detected by XP
            numXP{a,1} = length(ilbXPtimes);
            pctBobsXP{a,1} = numXP{a,1}/length(pks)*100;
            title(sprintf('%s %dpct XP captured',animal, round(pctBobsXP{a,1})))
            
%             figure
% %             plot(eye(:,1))
%             fs =  125;
% %             m = length(eye(:,1));       % original sample length
% %             n = pow2(nextpow2(m));  % transform length
% %             y = fft(eye(:,1),n);        % DFT of signal
% %             
% %             f = (0:n-1)*(fs/n)/10;
% %             power = abs(y).^2/n;
% %             eye_cp = eye(:,1);
% %             eye_cp(~eyeValid,1) = 0;
%             x = eye(eyeValid,1)';
%             n = length(x);          % number of samples
%             f = (0:n-1)*(fs/n);     % frequency range
%             power = abs(x).^2/n;    % power of the DFT
% 
%             plot(f,smoothdata(power, 'loess', 500))
%             xlabel('Frequency')
%             ylabel('Power')
% %             
% %             plot(f(1:floor(n/2)),power(1:floor(n/2)))
% %             xlabel('Frequency')
% %             ylabel('Power')
%             xlim([1 30])
% %             ylim([0 ])
        end
    end
end

               
%% PLOT polar histogram
if plotfigs
    if plot_pct_XPbobs
        
    end
    
    if plot_polarHist_pAn
        figname = sprintf('%s-pAn',Label);
        Pp=load_plotting_params({'defaults',figname});
        nrows = 2;
        ncols = 3;
        for a = 1:length(animals)
            ifig = init_plot(showfigs, Pp.position);
            animal = animals{a};
            
            % plot eye polar
            col = 1;
            row = 1;
            rc = (row-1)*(ncols)+1+(col-1);
            subaxis(nrows,ncols,rc,Pp.posparams{:},'SpacingHoriz',Pp.SpHz,...
                'SpacingVert', Pp.SpVt);
            
            polarhistogram(XP_eyeXphase_Sh{a}, Pp.bins, 'Normalization', ...
                'pdf', 'edgealpha', .1, 'facecolor', Pp.ShClr);
            hold on
            polarhistogram(XP_eyeXphase{a}, Pp.bins, 'Normalization', ...
                'pdf', 'edgealpha', .1, 'FaceColor', Pp.eyeClr, 'facealpha', .8);
            th = vecang_eye{a,1};
            r = meanMRVmag_eye{a,1};
            polarplot([0 th], [0 r], 'k')
            polarscatter(th, r, Pp.MRVsz,'ok', 'filled')
            rlim([0 Pp.rl])
            set(gca, 'FontSize', 10)
            thetaticks(0:45:315)
            rticks([.3 .6])
            title('eye')
            
            % plot eye phasemod, 95% intervals, real eye phasemod
            col = 1;
            row = 2;
            rc = (row-1)*(ncols)+1+(col-1);
            subaxis(nrows,ncols,rc,Pp.posparams{:},'SpacingHoriz',Pp.SpHz, ...
                'SpacingVert', Pp.SpVt);
            pmSh = phasemod_eyeX_sh{a}.phasemodSh;
            histogram(pmSh, 100, 'edgealpha', 0, 'facecolor', 'k');
            hold on;
            sigpct = .95;
            mrvsort = sort(pmSh);
            idxsig = round(sigpct*length(pmSh));
            line([mrvsort(idxsig) mrvsort(idxsig)], ylim, 'color', ...
                [0 0 0 .8], 'linestyle', '--', 'linewidth', 2);
            pm = phasemod_eyeX{a,1};
            line([pm pm], ylim, 'color', Pp.eyeClr, 'linewidth', 2);
            modp = phasemod_eyeX_sh{a}.modPctRank;
            pval = phasemod_eyeX_Rp{a,1};
            title(sprintf('p = %.02e', pval));
            ylabel('count')
            xlabel('log(Rayleigh Z)')
            axis tight
            xl = xlim;
            xlim([xl(1) xl(2)+diff(xlim)*.1])
            hold off
            set(gca, 'FontSize', 10)
            
            
            % plot ear polar
            col = 2;
            row = 1;
            rc = (row-1)*(ncols)+1+(col-1);
            subaxis(nrows,ncols,rc,Pp.posparams{:},'SpacingHoriz',Pp.SpHz, ...
                'SpacingVert', Pp.SpVt);
            polarhistogram(XP_earXphase_Sh{a}, Pp.bins, 'Normalization', 'pdf', ...
                'edgealpha', .1, 'facecolor', Pp.ShClr);
            hold on
            polarhistogram(XP_earXphase{a}, Pp.bins, 'Normalization', 'pdf', ...
                'edgealpha', .1, 'FaceColor', Pp.earClr, 'facealpha', .8);
            th = vecang_ear{a,1};
            r = meanMRVmag_ear{a,1};
            polarplot([0 th], [0 r], 'k')
            polarscatter(th, r, Pp.MRVsz,'ok', 'filled')
            rlim([0 Pp.rl])
            set(gca, 'FontSize', 10)
            thetaticks(0:45:315)
            rticks([.3 .6])
            title('ear')
            
            
            % plot ear phasemod, 95% intervals, real ear phasemod
            col = 2;
            row = 2;
            rc = (row-1)*(ncols)+1+(col-1);
            subaxis(nrows,ncols,rc,Pp.posparams{:},'SpacingHoriz',Pp.SpHz, ...
                'SpacingVert', Pp.SpVt);
            pmSh = phasemod_earX_sh{a}.phasemodSh;
            histogram(pmSh, 100, 'edgealpha', 0, 'facecolor', 'k');
            hold on;
            sigpct = .95;
            mrvsort = sort(pmSh);
            idxsig = round(sigpct*length(pmSh));
            line([mrvsort(idxsig) mrvsort(idxsig)], ylim, 'color', ...
                [0 0 0 .8], 'linestyle', '--', 'linewidth', 2);
            pm = phasemod_earX{a,1};
            line([pm pm], ylim, 'color', Pp.earClr, 'linewidth', 2);
            modp = phasemod_earX_sh{a}.modPctRank;
            pval = phasemod_earX_Rp{a,1};
            title(sprintf('p = %.02e', pval));
%             ylabel('count')
            xlabel('log(Rayleigh Z)')
            axis tight
            xl = xlim;
            xlim([xl(1) xl(2)+diff(xlim)*.1])
            hold off
            set(gca, 'FontSize', 10)            
            
            % plot nose polar
            col = 3;
            row = 1;
            rc = (row-1)*(ncols)+1+(col-1);
            subaxis(nrows,ncols,rc,Pp.posparams{:},'SpacingHoriz',Pp.SpHz, ...
                'SpacingVert', Pp.SpVt);
            
            polarhistogram(XP_noseXphase_Sh{a}, Pp.bins, 'Normalization', 'pdf', ...
                'edgealpha', .1, 'facecolor', Pp.ShClr);
            hold on
            polarhistogram(XP_noseXphase{a}, Pp.bins, 'Normalization', 'pdf', ...
                'edgealpha', .1, 'FaceColor', Pp.noseClr, 'facealpha', .8);
            
            th = vecang_nose{a,1};
            r = meanMRVmag_nose{a,1};
            polarplot([0 th], [0 r], 'k')
            polarscatter(th, r, Pp.MRVsz,'ok', 'filled')
            
            hold off
            rlim([0 Pp.rl])
            hold off
            set(gca, 'FontSize', 10)
            thetaticks(0:45:315)
            rticks([.3 .6])
            title('nose')
            
            % plot nose phasemod, 95% intervals, real nose phasemod
            col = 3;
            row = 2;
            rc = (row-1)*(ncols)+1+(col-1);
            subaxis(nrows,ncols,rc,Pp.posparams{:},'SpacingHoriz',Pp.SpHz, ...
                'SpacingVert', Pp.SpVt);
            pmSh = phasemod_noseX_sh{a}.phasemodSh;
            histogram(pmSh, 100, 'edgealpha', 0, 'facecolor', 'k');
            hold on;
            sigpct = .95;
            mrvsort = sort(pmSh);
            idxsig = round(sigpct*length(pmSh));
            line([mrvsort(idxsig) mrvsort(idxsig)], ylim, 'color', ...
                [0 0 0 .8], 'linestyle', '--', 'linewidth', 2);
            pm = phasemod_noseX{a,1};
            line([pm pm], ylim, 'color', Pp.noseClr, 'linewidth', 2);
            modp = phasemod_noseX_sh{a}.modPctRank;
            pval = phasemod_noseX_Rp{a,1};
            title(sprintf('p = %.0e', pval));
%             ylabel('count')
            xlabel('log(Rayleigh Z)')
            axis tight
            xl = xlim;
            xlim([xl(1) xl(2)+diff(xlim)*.1])
            hold off
            set(gca, 'FontSize', 10)
            
            % super
            stit = sprintf('%s %s', figname, animal);
            setSuperAxTitle(stit);
            text(0.01, 0.99, sprintf('pctBobsXP: %.02f; numXP: %d',pctBobsXP{a}, numXP{a}), 'FontSize', 8)
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
    
    if plot_polarHist
        figname = sprintf('%s',Label);
        Pp=load_plotting_params({'defaults',figname});
        ifig = init_plot(showfigs, Pp.position);
        
        % plot eye polar
        sf1 = subaxis(1, 3, 1, Pp.posparams{:});
        nk = 200;
        sub_XP_eyeXphase_Sh = cellfun(@(x) x(randsample(length(x),nk)), XP_eyeXphase_Sh,'un',0);
        polarhistogram(cell2mat(sub_XP_eyeXphase_Sh), Pp.bins, 'Normalization', 'pdf', ...
            'edgealpha', .1, 'facecolor', Pp.ShClr, 'DisplayName', 'shuffled');
        hold on
        sub_XP_eyeXphase = cellfun(@(x) x(randsample(length(x),nk)), XP_eyeXphase,'un',0);
        polarhistogram(cell2mat(sub_XP_eyeXphase), Pp.bins, 'Normalization', 'pdf', 'edgealpha', .1,...
            'facealpha', .8, 'facecolor', Pp.eyeClr, 'DisplayName', 'phasemod /Day');
        % overlay meanMRV of each an as point
        th = cell2mat(vecang_eye);
        r = cell2mat(meanMRVmag_eye);
        polarplot([zeros(length(th),1) th]', [zeros(length(r),1) r]', 'k')
        polarscatter(th, r, Pp.MRVsz, ...
            'ok', 'filled', 'markeredgecolor', 'w','markeredgealpha', .5);
        title(sprintf('eye range:%ddeg',round(rad2deg(range(th)))))
        rlim([0 Pp.rl])
        thetaticks(0:45:315)
        rticks([.3 .6])
        set(gca, 'FontSize', 10)
        
        % plot ear polar
        sf2 = subaxis(1, 3, 2, Pp.posparams{:});
        sub_XP_earXphase_Sh = cellfun(@(x) x(randsample(length(x),nk)), XP_earXphase_Sh,'un',0);
        polarhistogram(cell2mat(sub_XP_earXphase_Sh), Pp.bins, 'Normalization', 'pdf', ...
            'edgealpha', .1, 'facecolor',Pp.ShClr, 'DisplayName', 'shuffled');
        hold on
        sub_XP_earXphase = cellfun(@(x) x(randsample(length(x),nk)), XP_earXphase,'un',0);       
        polarhistogram(cell2mat(sub_XP_earXphase), Pp.bins, 'Normalization', 'pdf', 'edgealpha', .1,...
            'facealpha', .8, 'facecolor', Pp.earClr, 'DisplayName', 'phasemod /Day');
        % overlay meanMRV of each an as point
        th = cell2mat(vecang_ear);
        r = cell2mat(meanMRVmag_ear);
        polarplot([zeros(length(th),1) th]', [zeros(length(r),1) r]', 'k')
        polarscatter(th, r, Pp.MRVsz, ...
            'ok', 'filled', 'markeredgecolor', 'w','markeredgealpha', .5);
        title(sprintf('eye range:%ddeg',round(rad2deg(range(th)))))
        rlim([0 Pp.rl])
        thetaticks(0:45:315)
        rticks([.3 .6])
        set(gca, 'FontSize', 10)
        
        % plot nose polar
        sf3 = subaxis(1, 3, 3, Pp.posparams{:});
        sub_XP_noseXphase_Sh = cellfun(@(x) x(randsample(length(x),nk)), XP_noseXphase_Sh,'un',0);       
        polarhistogram(cell2mat(sub_XP_noseXphase_Sh), Pp.bins, 'Normalization', 'pdf', ...
            'edgealpha', .1, 'facecolor', Pp.ShClr, 'DisplayName', 'shuffled');
        hold on
        sub_XP_noseXphase = cellfun(@(x) x(randsample(length(x),nk)), XP_noseXphase,'un',0);       
        polarhistogram(cell2mat(sub_XP_noseXphase), Pp.bins, 'Normalization', 'pdf', 'edgealpha', .1,...
            'facealpha', .8, 'facecolor', Pp.noseClr, 'DisplayName', 'phasemod /Day');
        % overlay meanMRV of each an as point
        th = cell2mat(vecang_nose);
        r = cell2mat(meanMRVmag_nose);
        polarplot([zeros(length(th),1) th]', [zeros(length(r),1) r]', 'k')
        polarscatter(th, r, Pp.MRVsz, ...
            'ok', 'filled', 'markeredgecolor', 'w','markeredgealpha', .5);
        title(sprintf('eye range:%ddeg',round(rad2deg(range(th)))))
        rlim([0 Pp.rl])
        thetaticks(0:45:315)
        rticks([.3 .6 ])
        pax = gca;
        set(gca, 'FontSize', 10)
        
        % super
        stit = sprintf('%s %s', figname);
        setSuperAxTitle(stit);
        if pausefigs
            pause
        end
        if savefigs
            strsave = save_figure('demetris', stit, 'savefigas', savefigas, ...
                'subdir', figname);
        end
    end
end