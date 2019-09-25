function [ampl] = sj_do_ampl(spikes, nch)


channel_lth = size(spikes.waveforms,2)/nch;
nSpikes = size(spikes.waveforms,1);

for ch=1:nch
    if isfield(spikes,'waveforms_ch1')
        cmd=sprintf('ampl(:,ch)=abs(max(spikes.waveforms_ch%d,[],2)) + abs(min(spikes.waveforms_ch%d,[],2));',ch,ch); eval(cmd);
    else
        w = spikes.waveforms(:,channel_lth*(ch-1)+1:channel_lth*ch);
        ampl(:,ch) = abs(max(w,[],2)) + abs(min(w,[],2));
    end
end


