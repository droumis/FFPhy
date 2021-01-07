
function out = dfa_XPtrigAvgRip(idx, timeFilter, varargin)
% get XP trig avg of rip power
% use with singledayanal
% DR_2020-09-30

fprintf('day %d \n',idx(1))
reqData = {'ca1rippleskons', 'lick', 'pos'};
check_required(reqData, varargin)
wienerFiltDiv = 1e1;
circperms = 25;
win = [-2 2];
eventType = 'ca1rippleskons';
srate = 1500;
swin = 1500;
stp = 375;
do_circperms = 0;
if ~isempty(varargin)
    assign(varargin{:});
end

out = init_out(); % init output
out.index = idx;
out.animal = animal;
day = idx(1);
eps = idx(2:end);

%% get rip pwr trace
evid = find(contains(varargin(1:2:end), eventType));
o = [1:2:length(varargin)]+1;
swr = varargin{o(evid)};

rippwr = [];
time = [];
XP = [];
position = [];
for e = 1:eps
    rippwr = [rippwr; swr{day}{eps(e)}{1}.powertrace'];
    time = [time; swr{day}{eps(e)}{1}.eegtimesvec_ref'];
    % Remove events with windows exceeding epoch bounds
    iXP = lick{day}{eps(e)}.starttime;
    epST = swr{day}{eps(e)}{1}.eegtimesvec_ref(1);
    epET = swr{day}{eps(e)}{1}.eegtimesvec_ref(end);
    iXP(any([iXP<(epST+abs(win(1))) iXP>(epET-abs(win(2)))],2)) = [];
    XP = [XP; iXP];    
%     position = [position; pos{day}{eps(e)}.data];
end
XPidxt = lookup(XP, time');
XPt = zeros(length(time),1);
XPt(XPidxt)=1;

%% get times when the animal is within X distance from a well
% alternatively for now just take the lick bout intervals

[~, boutTimes] = getLickBoutLicks(animal, ...
    [repmat(day, length(eps),1) eps'], varargin{:});
boutTimesall = vertcat(boutTimes{day}{eps});
% boutTimesall = [boutTimesall(:,1)-1500 boutTimesall(:,2)+1500];
%%
% figure; histogram(diff(boutTimesall, [], 2),100)
xcfft = 0;
acfft = 0;
% rippwr_mm = rippwr - nanmean(rippwr);
% XPt_mm = XPt - nanmean(XPt);
for b = 1:size(boutTimesall,1)
    tInB = isIncluded(time, boutTimesall(b,:));
    if sum(tInB) < 1500
            continue
    end
%     sum(tInB)/1500
    rip_bout = rippwr(tInB);
    rip_bout = rip_bout - nanmean(rip_bout);
    xp_bout = XPt(tInB);
    xp_bout = xp_bout - nanmean(xp_bout);
    lenbout = length(rip_bout);
    v = 1:stp:lenbout;
    v = v(v<lenbout-(swin-1)); % crop to swin units
    
    for iv = 1:length(v)
        iripBout = rip_bout(v(iv):v(iv)+(swin-1));
        ixpBout = xp_bout(v(iv):v(iv)+(swin-1));
        if any(isnan(iripBout)) || any(isnan(ixpBout))
            fprintf('nan found, skipping time chunk\n')
            continue
        end
        fft_iripBout = fft(iripBout);
        fft_ixpBout = fft(ixpBout);

        w = conj(fft_iripBout) .* fft_iripBout; % weiner filter 
        % xcorr and acorr
        xcfft = xcfft + (conj(fft_ixpBout) .* fft_iripBout);
        acfft = acfft + (conj(fft_ixpBout) .* fft_ixpBout); %+ nanmean(w)/wienerFiltDiv); 
    end
end

xcnorm = ifft(xcfft ./ acfft);
xcnorm = [xcnorm(length(xcnorm)/2:end); xcnorm(1:length(xcnorm)/2)];
% figure(1); plot(xcnorm, 'color', 'b', 'linewidth', 2); hold on;
out.xcnorm = xcnorm;

%% Convolve the xcnorm with the BB pulse signal
C = conv(XPt, xcnorm, 'same');
pconf = paramconfig;
flab = 'xcNormACpredCorr';
savefigs = 1;
showfigs= 0;
plotfigs = 0;
coef_all = [];
for b = 1:size(boutTimesall,1)
    tInB = isIncluded(time, boutTimesall(b,:));
    if sum(tInB) < 1500
       continue
    end
    if sum(XPt(tInB)) < (10 * diff(boutTimesall(b,:),[], 2))
        continue
    end
    
    rip_bout = rippwr(tInB);
    rip_bout = rip_bout - nanmean(rip_bout);
    C_bout = C(tInB);
    C_bout = C_bout - nanmean(C_bout);
    if plotfigs
        t = time(tInB);
        ifig = init_plot(showfigs, [.1 .1 .9 .3]);
        yyaxis left
        plot(t, C_bout);
        ylabel('predicted rippwr')
        yyaxis right
        plot(t, rip_bout);
        ylabel('actual rippwr')
        title(sprintf('%s %s', animal, flab))
        axis tight
        if savefigs
            strsave = save_figure([pconf.andef{4} '/' flab],...
                sprintf('%s_day%d_%s_t%.03f', animal, day, flab, t(1)), 'savefigas', 'png');
        end
    end
    [coef, pval] = corr(C_bout, rip_bout);
    coef_all = [coef_all; coef];
end
out.ripPowPredCorrCoef = coef_all;

%% Circ perm * 500
if do_circperms
    for cp = 1:circperms
        xcfft = 0;
        acfft = 0;
        for b = 1:size(boutTimesall,1)
            tInB = isIncluded(time, boutTimesall(b,:));
            if sum(tInB) < 1500
                continue
            end
            rip_bout = rippwr(tInB);
            rip_bout = rip_bout - nanmean(rip_bout);
            
            cutidx = round(length(rip_bout)*rand(1));
            while cutidx==0 % prevent 0 index
                cutidx = round(length(rip_bout)*rand(1));
            end
            rip_bout = [rip_bout(cutidx:end); rip_bout(1:cutidx)];
            xp_bout = XPt(tInB);
            xp_bout = xp_bout - nanmean(xp_bout);
            %         xp_bout = circshift(xp_bout, round(length(xp_bout)*rand(1)));
            lenbout = length(rip_bout);
            v = 1:stp:lenbout;
            v = v(v<lenbout-(swin-1)); % crop to swin units
            
            for iv = 1:length(v)
                iripBout = rip_bout(v(iv):v(iv)+(swin-1));
                ixpBout = xp_bout(v(iv):v(iv)+(swin-1));
                if any(isnan(iripBout)) || any(isnan(ixpBout))
                    fprintf('nan found, skipping time chunk\n')
                    continue
                end
                fft_iripBout = fft(iripBout);
                fft_ixpBout = fft(ixpBout);
                
                w = conj(fft_iripBout) .* fft_iripBout; % weiner filter
                % xcorr and acorr
                xcfft = xcfft + (conj(fft_ixpBout) .* fft_iripBout);
                acfft = acfft + (conj(fft_ixpBout) .* fft_ixpBout + nanmean(w)/wienerFiltDiv);
            end
        end
        xcnorm = ifft(xcfft ./ acfft);
        xcnorm = [xcnorm(length(xcnorm)/2:end); xcnorm(1:length(xcnorm)/2)];
        %     figure(1); plot(xcnorm, 'color', 'k'); hold on;
        %     axis tight
        %     line([750 750], ylim, 'color', 'k')
        out.circperm_xcnorm{cp,1} = xcnorm;
    end
end
%%
% for each bout:
% 1: subtract mean from rip and XP traces
% top = zeros(winlen,1);
% bot = zeros(winlen,1);
% for each time win:
    % 2: FFT each signal
    % 3: top = top + conj(fft_xp) .* fft_rip;
    % 4  bot = bot + conj(fft_xp) .* fft_xp + w;
% end
% res = top ./ bot
% tres = ifft(res)

%% Full version
% mean subtract
% XPt_msub = XPt - nanmean(XPt);
% rippwr_msub = rippwr - nanmean(rippwr);
% % fft
% fft_XPt_msub = fft(XPt_msub);
% fft_rippwr_msub = fft(rippwr_msub);
% % xcorr and acorr
% xcfft = conj(fft_XPt_msub).*fft_rippwr_msub;
% acfft = conj(fft_XPt_msub).*fft_XPt_msub+nanmean(fft_XPt_msub)/10;
% % div norm
% xcfft_acdiv = xcfft ./ acfft;
% % ifft
% ifft_xc_acdiv = ifft(xcfft_acdiv);
% ifft_xc_acdiv = ifft_xc_acdiv';
% % stack peri-event trace
% wint = win(1)*srate:win(2)*srate;
% xptimeIdx = lookup(XP, time);
% ifft_xc_acdiv_XPtrig = cell2mat(arrayfun(@(x) ifft_xc_acdiv(cell2mat(x)), ...
%     arrayfun(@(y) y+wint, xptimeIdx, 'un', 0), 'un', 0));
% % output real
% ifft_xc_acdiv_XPtrig = real(ifft_xc_acdiv_XPtrig);
% out.ifft_xc_acdiv_XPtrig = ifft_xc_acdiv_XPtrig;
% figure; plot(nanmean(ifft_xc_acdiv_XPtrig)); hold on;
%% time slide version
% twins = 1500; % 1.5k is 1 sec
% tstep = 150; % 150 is .1 sec 

% rippwr = rippwr';
% stack peri-event trace
% rippwrstack = cell2mat(arrayfun(@(x) rippwr(cell2mat(x))', arrayfun(@(y) ...
%     y+wint, xptimeIdx, 'un', 0), 'un', 0));
% XPtstack = cell2mat(arrayfun(@(x) XPt(cell2mat(x))', arrayfun(@(y) y+wint,...
%     xptimeIdx, 'un', 0), 'un', 0));

% rippwrstack = rippwrstack';
% rippwrstackcat  = rippwrstack(:);

% XPtstack = XPtstack';
% XPtstackcat = XPtstack(:);
% steps = 1:length(wint):length(XPtstackcat);

% while ~(steps(end) + twins < length(rippwr))
%     steps(end) = [];
% end
% mean subtract using mean per window
% rippwrstack = rippwrstack - nanmean(rippwrstack,2);
% XPtstack = XPtstack - nanmean(XPtstack,2);
%%
% swin = 1500;
% stp = 500;
% lwin = size(rippwrstack,2);
% v = 1:stp:lwin;
% v = v(v<lwin-(swin-1));
% % nanitialize
% xcfft = nan(size(XPtstack,1), swin, length(v));
% acfft = nan(size(XPtstack,1), swin, length(v));
% % vectorized fft 
% % fft_rippwr = fft(rippwrstack, [], 2);
% % fft_XPt = fft(XPtstack, [], 2);
% 
% %%
% for ist = 1:length(xptimeIdx)
%     fprintf('%d\n', ist)
%     iri = rippwrstack(ist,:);
%     ixp = XPtstack(ist,:);
% 
%     for iv = 1:length(v)
%         iv_iri = iri(v(iv):v(iv)+(swin-1));
%         iv_ixp = ixp(v(iv):v(iv)+(swin-1));
%     
%         iv_iri = iv_iri - nanmean(iv_iri);
%         iv_ixp = iv_ixp - nanmean(iv_ixp);
%         
%         istfft_rippwr = fft(iv_iri);
%         istfft_XPt = fft(iv_ixp);
% %     fprintf('%d\n', ist)
% %     istfft_rippwr = fft_rippwr(ist,:);
% %     istfft_XPt = fft_XPt(ist,:);
%         % xcorr and acorr
%         w = conj(istfft_rippwr) .* istfft_rippwr;
% %         oi = (iv-1)*750+1;
%         xcfft(ist,:,iv) = conj(istfft_XPt) .* istfft_rippwr;
%         acfft(ist,:,iv) = conj(istfft_XPt) .* istfft_XPt + nanmean(w)/10; %nanmean(istfft_XPt)/1e2; %1e-10;
%     end
% end
%%
% % div norm
% % d = smoothdata(acfft,2, 'rlowess',150);
% xcfft_acdiv = xcfft ./ acfft;
% % nanmean the sliding windows per segment before ifft
% % xcfft_acdiv_2d = nanmean(xcfft_acdiv,3);
% xc_ifft_acdiv = ifft(xcfft_acdiv, [], 2);
% 
% o = nan(size(XPtstack,1), lwin, length(v));
% for s = 1:size(xc_ifft_acdiv,3)
%     % transform to o matrix with correct time offset
%     oi = (s-1)*stp+1;
%     o(:,oi:oi+(swin-1),s) = xc_ifft_acdiv(:,:,s);
% end
% d = nanmean(o,3);
% % real_xcfft_acdiv = xcfft_acdiv);
% % xcfft_acdiv = xcfft_acdiv';
% % step_ifft_xc_acdiv_XPtrig = cell2mat(arrayfun(@(x) real_xcfft_acdiv(cell2mat(x)), arrayfun(@(y) y+wint, xptimeIdx, 'un', 0), 'un', 0));
% out.xc_acdiv_pXP_mean = nanmean(d);
% figure; plot(xc_acdiv_pXP_mean);

 %%
% % figure; plot(m_ifft_xc_acdiv_XPtrig); line([3000 3000], ylim)
% % figure; imagesc(real(ifft_xc_acdiv_XPtrig));
% 
% % % apply timefilter to rippwr
% % incpos = ~isIncluded(position(:,1), timeFilter);
% % figure; plot(position(incpos,2), position(incpos,3))
% 
% % zscore
% rippwr = rippwr';
% % zrippwr = ((rippwr - nanmean(rippwr))./nanstd(rippwr));
% 
% %% get mean of all low velocity periods
% 
% Fp.params = {'<4cm/s', 'wtrackdays', 'excludePriorFirstWell', 'excludeAfterLastWell'};
% Fp = load_filter_params(Fp);
% F = createfilter('animal', {animal}, 'epochs', ...
%     Fp.epochfilter, 'excludetime', Fp.timefilter);
% bstf = [];
% epsidx = find(ismember(F.epochs{1}(:,1), day));
% for de = 1:length(epsidx)
%     % events in current timefilter
%     bstf = [bstf; F.excludetime{1}{epsidx(de)}];
% end
% 
% baseline_include_idx = ~isExcluded(time, bstf);
% meanbaseline = nanmean(rippwr(baseline_include_idx));
% stdbaseline = nanstd(rippwr(baseline_include_idx));
% % sembaseline = stdbaseline / sqrt(length(baseline_include_idx));
% 
% 
% %% Get licks in lick-burst intervals
% % [intraBoutXP, boutTimes] = getLickBoutLicks(animal, [repmat(day,length(eps),1) eps'], ...
% %     varargin{:});
% % boutTimes = cell2mat({boutTimes{day}{eps}}');
% % intraBoutXP = cell2mat({intraBoutXP{day}{eps}}');
% % intraBoutXP = intraBoutXP(:,1);
% % fprintf('%d XP within %d bursts \n', numel(intraBoutXP), size(boutTimes,1))
% 
% % for each iLB XP, get a zrippwr trace within window
% % return this stack
% % also return the nanmean and nanstd of the stack
% 
% pre = abs(win(1)*srate)-1;
% post = abs(win(2)*srate);
% cols = pre+post+1;
% zrippwr_XPtrig = nan(size(XP,1), cols);
% rippwr_XPtrig = nan(size(XP,1), cols);
% 
% 
% % stack peri-event trace
% wint = win(1)*srate:win(2)*srate;
% xptimeIdx = lookup(XP, time);
% rippwr_XPtrig = cell2mat(arrayfun(@(x) rippwr(cell2mat(x)), arrayfun(@(y) y+wint, xptimeIdx, 'un', 0), 'un', 0));
% 
% 
% % for iXP = 1:length(XP)
% %     xptimeIdx = min(find(time>=XP(iXP)));
% %     if diff([XP(iXP) time(xptimeIdx)]) > .001
% %         sprintf('%.03f diff, skipping',diff([XP(iXP) time(xptimeIdx)]))
% %         continue
% %     end
% %     startIdx = xptimeIdx-pre;
% %     endIdx = xptimeIdx+post;
% %     try
% % %         zrippwr_XPtrig(iXP,:) = zrippwr(startIdx:endIdx);
% %         rippwr_XPtrig(iXP,:) = rippwr(startIdx:endIdx);
% %     catch
% %         disp('window exceeds bounds, skipping')
% %         continue
% %     end
% % end
% %% output
% 
% out.meanbaseline = meanbaseline;
% out.stdbaseline = stdbaseline;
% 
% 
% % out.rippwr = rippwr;
% % out.intraBoutXP = XP;
% 
% % out.zrippwr_XPtrig = zrippwr_XPtrig;
% % out.mean_zrippwr_XPtrig = nanmean(zrippwr_XPtrig,1);
% % out.sem_zrippwr_XPtrig = nanstd(zrippwr_XPtrig,1)/sqrt(size(zrippwr_XPtrig,1));
% 
% mAll = nanmean(nanmean(rippwr_XPtrig,1));
% rippwr_XPtrig_mcenter = rippwr_XPtrig - mAll;
% 
% out.rippwr_XPtrig = rippwr_XPtrig_mcenter;
% out.mean_rippwr_XPtrig = nanmean(rippwr_XPtrig_mcenter,1);
% out.sem_rippwr_XPtrig = nanstd(rippwr_XPtrig_mcenter,1)/sqrt(size(rippwr_XPtrig_mcenter,1));
% 
% 
% 
% % XPt = XPt - nanmean(XPt);
% [xc,~] = xcorr(XPt,win(2)*srate,'coeff');
% out.XPacorr = xc;
% 
% % rippwr = rippwr-nanmean(rippwr);
% [xc,~] = xcorr(rippwr-nanmean(rippwr),win(2)*srate,'coeff');
% out.ripAcorr = xc;
% 
% incRPidx = ~isIncluded(time, timeFilter); 
% rippwrNearWell = rippwr(incRPidx);
% rippwrNearWell = rippwrNearWell-nanmean(rippwrNearWell);
% rippwrAwayWell = rippwr(~incRPidx);
% rippwrAwayWell = rippwrAwayWell-nanmean(rippwrAwayWell);
% 
% % out.boutTimes = boutTimes;
% % boutTimesVec = isIncluded(time, boutTimes);
% [xc,~] = xcorr(rippwrNearWell,win(2)*srate,'coeff'); 
% out.nearW_ripAcorr = xc;
% % figure; plot(xc); hold on
% 
% [xc,~] = xcorr(rippwrAwayWell,win(2)*srate,'coeff'); 
% out.awayW_ripAcorr = xc;
% % plot(xc); 
% 
% 
% 
% % [xc,~] = xcorr(rippwr(~boutTimesVec),win(2)*srate,'coeff');
% % out.eLB_ripAcorr = xc;
% % plot(xc);
out.win = win;
out.srate = srate;
% out.time = linspace(win(1),win(2),length(out.xcnorm));

end

function out = init_out()
out.index = [];
out.animal = [];
% out.time = [];
% 
% out.meanbaseline = [];
% out.stdbaseline = [];
% 
% % out.rippwr = [];
% % out.intraBoutXP = [];
% 
% % out.zrippwr_XPtrig = [];
% % out.mean_zrippwr_XPtrig = [];
% % out.sem_zrippwr_XPtrig = [];
% out.rippwr_XPtrig = [];
% out.mean_rippwr_XPtrig = [];
% out.sem_rippwr_XPtrig = [];

out.win = [];
out.srate = [];
out.time = [];

% out.XPacorr = [];
% out.ripAcorr = [];
% 
% out.awayW_ripAcorr = [];
% out.nearW_ripAcorr = [];

out.xcnorm = [];
out.ripPowPredCorrCoef = [];
end