function spikes=sss_splice_waves(spikes)

spikes.waveforms_ch1(:,1:2)=[]; spikes.waveforms_ch1(:,end-1:end)=[];
spikes.waveforms_ch2(:,1:2)=[]; spikes.waveforms_ch2(:,end-1:end)=[];
spikes.waveforms_ch3(:,1:2)=[]; spikes.waveforms_ch3(:,end-1:end)=[];
spikes.waveforms_ch4(:,1:2)=[]; spikes.waveforms_ch4(:,end-1:end)=[];

spikes.waveforms=[spikes.waveforms_ch1 spikes.waveforms_ch2 spikes.waveforms_ch3 spikes.waveforms_ch4];
