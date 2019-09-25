[xc, S_unit, S_lfp, f] = spike_lfp_frequency_shift(unit,lfp,timeranges,params)
%SPIKE_LFP_FREQUENCY_SHIFT Estimate and compare spectra of simultaneous single-unit spiking and LFP.
%
%   [XC, S_UNIT, S_LFP, F] =
%   SPIKE_LFP_FREQUENCY_SHIFT(UNIT,LFP,TIMERANGES,PARAMS)
%   
%   where PARAMS = [W T P], W is bandwidth, T is timestep == window duration, P
%   is number of tapers to drop such that (2*W*T - P) tapers are used
%
%
%Depends on:
%   IS_UNIT (written by smk)
%   IS_CONTINUOUS (written by smk)
%
%Written by SMK, 2009 November 23.
%

TS_PER_SEC = 1e4;


for i = 1:numel(unit)

  Fs = lfp(i).Fs;

  N = numel(lfp(i).samples);

  mintime = double(unit(i).timerange)/TS_PER_SEC;
  maxtime = double(unit(i).timerange)/TS_PER_SEC;
  tn = (mintime+params(2)/2):params(2):(maxtime-params(2));
  Nwin = round(Fs*params(2));
  nfft = 2^ceil(log2(Nwin));
  tapers = dpsschk([params(1)*params(2), 2*params(1)*params(2)-params(3)],Nwin,Fs);

  winstart = 1:Nwin:N-Nwin+1;
  num_windows = length(winstart); 

  for j = 1:length(winstart)
    indx = winstart(j):winstart(j)+

    t = linspace(tn(j)-params(2)/2,tn(j)+params
  end

for n=1:nw;
   indx=winstart(n):winstart(n)+Nwin-1;
   datawin=data(indx,:);
   if nargout==4
     [s,f,serr]=mtspectrumc(datawin,params);
     Serr(1,n,:,:)=squeeze(serr(1,:,:));
     Serr(2,n,:,:)=squeeze(serr(2,:,:));
   else
     [s,f]=mtspectrumc(datawin,params);
   end
   S(n,:,:)=s;
end;
S=squeeze(S); 
if nargout==4;Serr=squeeze(Serr);end;
winmid=winstart+round(Nwin/2);
t=winmid/Fs;

  [S_unit{i},t{i},f] = mtspecgramc(lfp(i).samples,movingwin,params);
  [S_lfp{i},t{i},f] = mtspecgramc(lfp(i).samples,movingwin,params);


end

% Pick out windows that fall within timeranges and average over these

