function [S,t,f,R,Serr]=mtspecgramtrigpt(data,E,win,movingwin,params,fscorr)
% Multi-taper event triggered time-frequency spectrum - point process times
%
% Usage:
%
% [S,t,f,R,Serr]=mtspecgramtrigpt(data,E,win,movingwin,params,fscorr)
% Input: 
%       data        (univariate spike data -1d structure array of spike times; also accepts 1d array of spike times) -- required
%       E           (event times) - required
%       win         (in the form [winl,winr] i.e window around each event
%                                                 required
%       movingwin         (in the form [window winstep] i.e length of moving
%                                                 window and step size) -
%                                                 required
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       - optional
%           tapers (precalculated tapers from dpss, or in the form [NW K] e.g [3 5]) -- optional. If not 
%                                                 specified, use [NW K]=[3 5]
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over events when 1, don't average when 0) - optional. Default 0
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
%
% Output:
%       S       (Spectrogram with dimensions time x frequency x events if trialave=0; dimensions time x frequency if trialave=1)
%       t       (times)
%       f       (frequencies)
%       R       (spike rate)
%       Serr    (error bars)-only if err(1)>=1

if nargin < 4; error('Need data, events and parameters for the windows'); end;
if nargin < 5; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
clear tapers pad Fs fpass err trialave
if nargin < 6 || isempty(fscorr); fscorr=0; end;

data=createdatamatpt(data,E,win);
if nargout==5; 
    [S,t,f,R,Serr]=mtspecgrampt(data,movingwin,params,fscorr);
else 
    [S,t,f,R]=mtspecgrampt(data,movingwin,params,fscorr);
end;
