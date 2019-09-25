
function gui_waveplot(spikes,clusters,unqcs,arfp, plot_ampl, plot_sd, plot_type)

%--------------------------------------------------------------------------

ctr=1;
for plotclu=1:ceil(length(unqcs)/10)
    figure ([999+plotclu]);
    redimscreen100(100);
    %figure (2000)
    if length(unqcs) >= ctr+9
        tempcls=unqcs(ctr:ctr+9);
    else
        tempcls=unqcs(ctr:end);
    end

    tctr=1;
    for n_tempcls=1:length(tempcls)
        sp_idx=find (clusters==tempcls(n_tempcls));
        
         m = mean(spikes.waveforms(sp_idx,:));
         s = std(spikes.waveforms(sp_idx,:));

        % raw waveform
        subplot (5,2,tctr); hold on;  %subplot (length(tempcls),5,tctr);
        if plot_type==0
            plot (spikes.waveforms(sp_idx,:)', 'k'); axis tight
        else
            [n,x,y] = hist2d(spikes.waveforms(sp_idx,:));
            imagesc(x,y,n); axis xy; colormap hot;
        end
        plot(m, 'r', 'LineWidth', 2);
        if plot_sd==1
            plot(m+2*s, 'r.', 'LineWidth', 3);
            plot(m-2*s, 'r.', 'LineWidth', 3);
        end
        axis([0 size(spikes.waveforms,2) -plot_ampl plot_ampl])
        ylabel (['CLU ' num2str(tempcls(n_tempcls))])
        dtimes=diff(spikes.fstimes(sp_idx)/1000);
        yaxisval=get(gca,'Ylim');
        text (1,yaxisval(1,2)*0.6,['Sp<ISI=' num2str(length(find ((dtimes >= 0) & (dtimes <= arfp/1000))))],'Color','g');
        text (1,yaxisval(1,2)*0.8,['Nsp=' num2str(length(sp_idx))],'Color','g');
        set(gca,'xtick',[], 'xticklabel',[]);
        %set(gca,'ytick',[], 'yticklabel',[]);
        
        tctr=tctr+1;
        clear sp_idx
    end    
    ctr=ctr+10; clear tctr
end

i=1;