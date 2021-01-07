function xcNormAc = xcNormAc(sA, sB, varargin)
% calculate cross correlation normalized by autocorrelation
% using sliding window FFT approach
%
% input1: sA (series m x 1; this series is the one used also for autocorr)
% input2: sB (series m x 1)
%
% output: xcNormAc
%
% DKR 2021-01-01
%%
bin = 1500; % ms window for individual xc, ac
step = 250;
0; % ms to shift for each bin start
use_wienerFilt = 0;
wFiltDiv = 10; % weiner filter divisor

if ~isempty(varargin)
   assign(varargin{:}) 
end

% mean norm 
sA = sA - nanmean(sA);
sB = sB - nanmean(sB);

v = 1:step:length(sA);
v = v(v<(length(sA)-(bin-1))); % drop incomplete bin at end

xcfft = 0;
acfft = 0;
% slide along sA and sB, collecting xc and ac in freq space
for iv = 1:length(v)
    staIdx = v(iv);
    endIdx = v(iv)+(bin-1);
    isA = sA(staIdx:endIdx);
    isB = sB(staIdx:endIdx);
    if any([isnan(isA) isnan(isB)])
        fprintf('nan detected, omitting bin %d \n', staIdx)
        continue
    end
    isA_fft = fft(isA);
    isB_fft = fft(isB);
    % cross correlation
    xcfft = xcfft + (conj(isA_fft) .* isB_fft);
    % auto correlation
    if use_wienerFilt
        w = nanmean(conj(isB_fft) .* isB_fft)./wFiltDiv;
        acfft = acfft + (conj(isA_fft) .* isA_fft) + w;
    else
        acfft = acfft + (conj(isA_fft) .* isA_fft);
    end
end
% div norm xc with ac then inverse fft 
xcNormAc = ifft(xcfft ./ acfft);
% reorder second half in front of first
hlen = length(xcNormAc)/2;
xcNormAc = [xcNormAc(hlen+1:end); xcNormAc(1:hlen)];
% plot(xcNormAc)
% line([hlen hlen], ylim,'color','k')
end