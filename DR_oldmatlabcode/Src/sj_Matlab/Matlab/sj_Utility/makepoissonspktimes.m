function spiketimes=makepoissonspktimes(length,resolution,varargin)
%usage: spiketimes=makepoissonspktimes(length,resolution,varargin)
%
%input is length of desired spike train in ms and resolution of spiketimes in ms
%
%varargin can be equal to:
%'sine' which overlays on top of the poisson spiketrain a probability of
%   firing which changes as a sine function
%output is column vector of spike times

if exist('resolution','var')
    %do nothing
else
    resolution=.1;
end

spiketrain=random('poisson',.001,length/resolution,1);

sum(spiketrain==1)
sum(spiketrain==2)
sum(spiketrain==3)

spiketrain(spiketrain>1)=1;

if any(strcmp(varargin,'sine'))
   period=30000; %in bins, not ms
       sinetrain=sin((1:size(spiketrain,1))*2*pi/period)/2+.5;
   sinetrain=(sinetrain+(rand(1,size(spiketrain,1))-.5))';
   spiketrain=spiketrain.*(sinetrain>.5);
end

spiketimes=spktrain2spktime(spiketrain,resolution);

return

