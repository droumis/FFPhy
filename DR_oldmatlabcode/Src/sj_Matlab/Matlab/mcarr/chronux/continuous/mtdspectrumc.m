function [dS,f]=mtdspectrumc(data,phi,params)
% Multi-taper frequency derivative of the spectrum - continuous process
%
% Usage:
%
% [dS,f]=mtdspectrumc(data,phi,params)
% Input: 
%   Note that all times can be in arbitrary units. But the units have to be
%   consistent. So, if E is in secs, win, t have to be in secs, and Fs has to
%   be Hz. If E is in samples, so are win and t, and Fs=1. In case of spike
%   times, the units have to be consistent with the units of data as well.
%       data        (in form samples x channels/trials or a single vector) -- required
%       phi         (angle for evaluation of derivative) -- required.
%                       e.g. phi=[0,pi/2] gives the time and frequency derivatives
%       params: structure with fields tapers, pad, Fs, fpass, trialave
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
%           trialave (average over trials/channels when 1, don't average when 0) - optional. Default 0
% Output:
%       dS       (spectral derivative in form phi x frequency x channels/trials if trialave=0 or in form phi x frequency if trialave=1)
%       f        (frequencies)

if nargin < 2; error('Need data and angle'); end;
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
clear err params
data=change_row_to_column(data);
N=size(data,1);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers
K=size(tapers,2);
J=mtfftc(data,tapers,nfft,Fs);
J=J(findx,:,:);
A=sqrt(1:K-1);
A=repmat(A,[size(J,1) 1]);
A=repmat(A,[1 1 size(J,3)]);
S=squeeze(mean(J(:,1:K-1,:).*A.*conj(J(:,2:K,:)),2));
if trialave; S=squeeze(mean(S,2));end;
nphi=length(phi);
for p=1:nphi;
    dS(p,:,:)=real(exp(i*phi(p))*S);
end;
dS=squeeze(dS);
dS=change_row_to_column(dS);
