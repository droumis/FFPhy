function sss_guiTC(spikes, show_clu, plot_ampl, arfp)

spktimes=spikes.ftimes;
assigns=spikes.hierarchy.assigns;
unqstim=unique(spikes.igorstim);
stimonset=spikes.stimonset*10;  %% Convert to 0.1ms resolution

% onset_res=zeros(length(show_clu),length(unqstim))
% for plot_clu=1:length(show_clu)
%     for stim=1:length(unqstim)
%         sweep_idxs=find((assigns==unqcs)&(spikes.igorstim==unqstim));
%         clu_spktimes=spktimes(sweep_idxs);
%         onsettimes=clu_spktimes(find ((clu_spktimes>=stimonset) & (clu_spktimes<=stimonset+1000)));
%         bcktimes=clu_spktimes(find ((clu_spktimes<=stimonset) & (clu_spktimes>=stimonset-1000)));
%         onset_res(plot_clu,stim)= (length(onsettimes)-length(bcktimes))/length(sweep_idxs);
%     end
%
%     onset_prob=onset_res(plot_clu,:)/max(onset_res(plot_clu,:));
% end


ctr=1;
for lpfs=1:ceil(length(show_clu)/5)
    figure ([49+lpfs])

    if length(show_clu) >= ctr+4
        tempcls=show_clu(ctr:ctr+4);
    else
        tempcls=show_clu(ctr:end);
    end
    tctr=1;
    for lptempcls=1:length(tempcls)

        sp_idx=find (assigns==tempcls(lptempcls));
        m = mean(spikes.waveforms(sp_idx,:));
        s = std(spikes.waveforms(sp_idx,:));

        % RAW WAVEFORM
        subplot (5,3,tctr); hold on;
        plot (spikes.waveforms(sp_idx,:)', 'k'); axis tight
        plot(m, 'b', 'LineWidth', 3);
        plot(m+2*s, 'r.', 'LineWidth', 2);
        plot(m-2*s, 'r.', 'LineWidth', 2);
        axis([0 size(spikes.waveforms,2) -plot_ampl plot_ampl])
        if lptempcls ==1; title ('Raw Waves');end
        ylabel (['CL- ' num2str(tempcls(lptempcls))])
        dtimes=diff(spikes.fstimes(sp_idx)/1000);
        yaxisval=get(gca,'Ylim');
        text (1,yaxisval(1,2)*0.6,['Sp<ISI=' num2str(length(find ((dtimes >= 0) & (dtimes <= arfp/1000))))]);
        text (1,yaxisval(1,2)*0.8,['Nsp=' num2str(length(sp_idx))]);

        % TUNING CURVE
        subplot (5,3,tctr+1); hold on;

        for stim=1:length(unqstim)
            sweep_idxs=find((assigns==tempcls(lptempcls))&(spikes.igorstim==stim));
            stim_swidxs=find(spikes.Alligorstim==stim);
            stim_spktimes=spktimes(sweep_idxs);
            onsettimes=stim_spktimes(find ((stim_spktimes>=stimonset) & (stim_spktimes<=stimonset+1000)));
            bcktimes=stim_spktimes(find ((stim_spktimes<=stimonset) & (stim_spktimes>=stimonset-1000)));
            onset_res(stim)= (length(onsettimes)-length(bcktimes))/length(stim_swidxs);
        end
        bar(onset_res);
        if lptempcls ==1; title ([ 'Onset Spks/stim: ' num2str(length(stim_swidxs)) ' Sweeps' ]); end
        set(gca, 'Xtick', [1:length(unqstim)], 'Xticklabel', [unqstim]);
        clear onset_res;

        % RASTER
        subplot (5,3,tctr+2);
        ttime = spikes.ftimes(sp_idx);
        rastersw = spikes.sweep(sp_idx);
        % rastersw = collapse_sweeps(rastersw);   DO NOT COLLAPSE SWEEPS
        plot (ttime, rastersw, 'b.')
        ylimit = length(rastersw);
        if ylimit <= 1, ylimit=5; end;
        if (exist('spikes.sweepall')), sweeps=spikes.sweepall; else, sweeps=rastersw; end;
        axis ([stimonset-2000 stimonset+2000 min(spikes.sweepall) max(spikes.sweepall)])
        if lptempcls == 1; title ('Overall'); elseif lptempcls == length(tempcls); xlabel ('Time (ms)'); end
        tctr=tctr+3;
    end
    ctr=ctr+5;

    redimscreen70s;
end

