function [C,phi,S12,S1,S2,t,f,zerosp,confC,phistd,Cerr]=cohgrampt(data1,data2,movingwin,params,fscorr)
% Multi-taper time-frequency coherence - two point processes given as times
% process 
%
% Usage:
%
% [C,phi,S12,S1,S2,t,f,zerosp,confC,phistd,Cerr]=cohgrampt(data1,data2,movingwin,params,fscorr)
% Input: 
% Note units have to be consistent. Thus, if movingwin is in seconds, Fs
% has to be in Hz. see chronux.m for more information.
%
%       data1  (structure array of spike times with dimension trials; also accepts 1d array of spike times) -- required
%       data2  (structure array of spike times with dimension trials; also accepts 1d array of spike times) -- required
%       movingwin (in the form [window winstep] -- required
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
%           trialave (average over trials when 1, don't average when 0) - optional. Default 0
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
% Output:
%       C (magnitude of coherency time x frequencies x trials for trialave=0; time x frequency for trialave=1)
%       phi (phase of coherency time x frequencies x trials for no trial averaging; time x frequency for trialave=1)
%       S12 (cross spectrum - time x frequencies x trials for no trial averaging; time x frequency for trialave=1)
%       S1 (spectrum 1 - time x frequencies x trials for no trial averaging; time x frequency for trialave=1)
%       S2 (spectrum 2 - time x frequencies x trials for no trial averaging; time x frequency for trialave=1)
%       t (time)
%       f (frequencies)
%       zerosp (1 for windows and trials where spikes were absent (in either channel),zero otherwise)
%       confC (confidence level for C at 1-p %) - only for err(1)>=1
%       phistd - theoretical/jackknife (depending on err(1)=1/err(1)=2) standard deviation for phi - Note that 
%                phi + 2 phistd and phi - 2 phistd will give 95% confidence
%                bands for phi - only for err(1)>=1 
%       Cerr  (Jackknife error bars for C - use only for Jackknife - err(1)=2)

if nargin < 3; error('Need data1 and data2 and window parameters'); end;
if nargin < 4; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargin < 5 || isempty(fscorr); fscorr=0; end;

if nargout > 10 && err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end;
if nargout > 8 && err(1)==0;
    error('Errors computed only if err(1) is not equal to zero');
end;

[N,Ch]=check_consistency(data1,data2);
[mintime1,maxtime1]=minmaxsptimes(data1);
[mintime2,maxtime2]=minmaxsptimes(data2);
mintime=min(mintime1,mintime2);
maxtime=max(maxtime1,maxtime2);

tn=mintime+movingwin(1)/2:movingwin(2):maxtime-movingwin(1)/2;
Nwin=round(Fs*movingwin(1)); % number of samples in window
% Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=max(2^(nextpow2(Nwin)+pad),Nwin);
f=getfgrid(Fs,nfft,fpass); Nf=length(f);
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers
nw=length(tn);
if trialave;
   C=zeros(nw,Nf);
   S12=zeros(nw,Nf);
   S1=zeros(nw,Nf);
   S2=zeros(nw,Nf);
   phi=zeros(nw,Nf);
   Cerr=zeros(2,nw,Nf);
%    phierr=zeros(2,nw,Nf);
   phistd=zeros(nw,Nf);
else
   C=zeros(nw,Nf,Ch);
   S12=zeros(nw,Nf,Ch);
   S1=zeros(nw,Nf,Ch);
   S2=zeros(nw,Nf,Ch);
   phi=zeros(nw,Nf,Ch);
   Cerr=zeros(2,nw,Nf,Ch);
%    phierr=zeros(2,nw,Nf,Ch);
   phistd=zeros(nw,Nf,Ch);
end;
zerosp=zeros(nw,Ch);

for n=1:nw;
   t=linspace(tn(n)-movingwin(1)/2,tn(n)+movingwin(1)/2,Nwin);
   datawin1=extractdatapt(data1,[t(1) t(end)]);datawin2=extractdatapt(data2,[t(1) t(end)]);
   if nargout==11;
     [c,ph,s12,s1,s2,f,zsp,confc,phie,cerr]=coherencypt(datawin1,datawin2,params,fscorr,t);
%      phierr(1,n,:,:)=squeeze(phie(1,:,:));
%      phierr(2,n,:,:)=squeeze(phie(2,:,:));
     phistd(n,:,:)=phie;
     Cerr(1,n,:,:)=squeeze(cerr(1,:,:));
     Cerr(2,n,:,:)=squeeze(cerr(2,:,:));
   elseif nargout==10;
     [c,ph,s12,s1,s2,f,zsp,confc,phie]=coherencypt(datawin1,datawin2,params,fscorr,t);
%      phierr(1,n,:,:)=squeeze(phie(1,:,:));
%      phierr(2,n,:,:)=squeeze(phie(2,:,:));
     phistd(n,:,:)=phie;
   else
     [c,ph,s12,s1,s2,f,zsp]=coherencypt(datawin1,datawin2,params,fscorr,t);
   end;
   C(n,:,:)=c;
   phi(n,:,:)=ph;
   S12(n,:,:)=s12;
   S1(n,:,:)=s1;
   S2(n,:,:)=s2;
   zerosp(n,:)=zsp;
end;
t=tn;
C=squeeze(C); phi=squeeze(phi);S12=squeeze(S12); S1=squeeze(S1); S2=squeeze(S2);zerosp=squeeze(zerosp);
if nargout > 9; confC=confc; end;
if nargout==11;Cerr=squeeze(Cerr);end;
% if nargout==10; phierr=squeeze(phierr);end
if nargout==10; phistd=squeeze(phistd);end
