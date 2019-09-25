function spikes = sss_cleancut_tet(spikes, n_ch, clusters)

%% Modified from sss_windowcut_tet
allwaveforms = spikes.waveforms;
stimes = spikes.fstimes';
spikes.noiseidx=[];
amplitude=400;

for nclu=1:length(clusters)
    currclu = clusters(nclu)

    assigns = spikes.hierarchy.assigns;  %%% UPDATE AT EACH STEP AFTER CLEANING EACH CLUSTER
   
    for i=3:3
        Tlines={'Amplitude?'};
        defent={'0.5';'1'}; % default entries
        infonames= {'amplitude'};
        info = inputdlg(Tlines, 'What Amplitude to use?' , 1, defent);

        if ~isempty(info)              %see if user hit cancel
            info = cell2struct(info,infonames);
            %n_ch = str2num(info.n_channels);   %convert string to number
            amplitude = str2num(info.amplitude);   %convert string to number
        end

        % Plot the current cluster and then ask user what should be the windowcut limits
        curridxs = find(assigns==currclu);   %%% Assigns updated for each channel too
        cmd=sprintf('swaveforms = spikes.waveforms_ch%d(curridxs,:);',i); eval(cmd);

        NF = figure; hold on;
        plot (swaveforms','b'); grid on;
        axis([0 size(swaveforms,2) -amplitude amplitude]);
        plot (mean (swaveforms), 'y', 'Linewidth', 4)
        redimscreen

        disp ('A graph of the raw waveforms for current channel are being plotted please consult the figure and return to the command window')
        [rng,lims]=ginput;

        if exist('lims')==0
            lims = input ('What are the limits of the WINDOW DISCRIMINATOR that you want to use?.....');
        else
            lims=sort(lims); lims=lims';
        end

        hold off

        lim1= min(lims); % lower boundary of the window
        lim2= max(lims); % upper boundary

        if exist('rng')==0
            rng = input ('What are the sampling point ranges you want to scan for the thresholds (i.e. [14 16])...');
        else
            rng=sort(round(rng));
        end
                
        jbfill(rng(1):rng(2),lim1*ones(size([rng(1):rng(2)])),-amplitude*ones(size([rng(1):rng(2)])),[0,1,0],[0,1,0],1,0.7);
        jbfill(rng(1):rng(2),amplitude*ones(size([rng(1):rng(2)])),lim2*ones(size([rng(1):rng(2)])),[0,1,0],[0,1,0],1,0.7);

        disp('_______________________________________________________________________________________________________')
        disp(['The voltage values you entered: ' num2str((lims)) ' uV at A/D point: ' num2str((rng'))]);
        disp('_______________________________________________________________________________________________________')

        keyboard;
        
        [spks, dpts] = find ( (swaveforms(:,(rng(1):rng(2)))>lim2) | (swaveforms(:,(rng(1):rng(2))) <lim1) );
        
        cutidxs=curridxs(spks);
        
        if (length(spks)~=0)
            swaveforms (unique(spks),:)=[];
            stimes (unique(spks))=[];
            assigns(unique(cutidxs))=[];

            spikes.waveforms(unique(cutidxs),:)=[];
            spikes.fstimes(unique(cutidxs))=[];
            spikes.hierarchy.assigns(unique(cutidxs))=[];

            for ch=1:n_ch
                cmd=sprintf('nwaveforms = spikes.waveforms_ch%d;',ch); eval(cmd);
                nwaveforms(unique(cutidxs),:)=[];
                cmd=sprintf('spikes.waveforms_ch%d = nwaveforms;',ch); eval(cmd);
            end

            spikes.ftimes(unique(cutidxs))=[];
            spikes.swtimes= spikes.fstimes/1000;
            spikes.spiketimes= spikes.ftimes/1000;
            spikes.noiseidx=[spikes.noiseidx; unique(cutidxs)];
            
            curridxs(spks)=[];
        end

        close(NF);
    end %%% END Nchannel

%     NF = figure; hold on;
%     plot (swaveforms','b');
    %axis tight; grid on
    plot (mean (swaveforms), 'y', 'Linewidth', 4);
    ylabel ('Voltage (uV)');
    xlabel ('Sample (32 pts/ms)')
    text (0.7*size(swaveforms,2), 0.8*max(max(swaveforms)), ['Clu: ' num2str(currclu)  '; #Sp: ' num2str(length(curridxs))]);

    clear spks dpts z t opt

    
end  %%% END Clusters



