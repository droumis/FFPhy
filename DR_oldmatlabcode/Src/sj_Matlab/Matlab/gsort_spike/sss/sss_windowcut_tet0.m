function spikes = sss_windowcut_tet(spikes)

allwaveforms = spikes.waveforms;
stimes = spikes.fstimes';
spikes.noiseidx=[];

Tlines={'How Many Channels?',['Max Amplitude on Plot?']};
defent={'4','400'}; % default entries
infonames= {'n_channels','amplitude'};
info = inputdlg(Tlines, 'How many channels? And Amplitude?', 1, defent);

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    n_ch = str2num(info.n_channels);   %convert string to number
    amplitude = str2num(info.amplitude);   %convert string to number
else
    return
end

for i=1:n_ch

    % Plot the waveforms and then ask user what should be the windowcut limits
    cmd=sprintf('swaveforms = spikes.waveforms_ch%d;',i); eval(cmd);
    NF = figure; hold on;
    plot (swaveforms'); grid on;
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

    disp('_______________________________________________________________________________________________________')
    disp(['The voltage values you entered: ' num2str(round(lims)) ' uV at A/D point: ' num2str(round(rng'))]);
    disp('_______________________________________________________________________________________________________')



    [spks, dpts] = find ( (swaveforms(:,(rng(1):rng(2)))>lim2) | (swaveforms(:,(rng(1):rng(2))) <lim1) );

    if (length(spks)~=0)
        swaveforms (unique(spks),:)=[];
        stimes (unique(spks))=[];
    end

    if isempty(swaveforms);
        disp ('WARNING! There is no spike within the window specified. Please rerun the program with new values')
        return
    end

    %spikes.waveforms=swaveforms;
    spikes.waveforms(unique(spks),:)=[];
    spikes.fstimes=stimes';
    %spikes.fstimes(unique(spks))=[];

    for ch=1:n_ch
        cmd=sprintf('nwaveforms = spikes.waveforms_ch%d;',ch); eval(cmd);
        nwaveforms(unique(spks),:)=[];
        cmd=sprintf('spikes.waveforms_ch%d = nwaveforms;',ch); eval(cmd);
        %spikes.waveforms_ch1(unique(spks),:)=[];
        %spikes.waveforms_ch2(unique(spks),:)=[];
        %spikes.waveforms_ch3(unique(spks),:)=[];
        %spikes.waveforms_ch4(unique(spks),:)=[];
    end

    spikes.ftimes(unique(spks))=[];

    spikes.sweep(unique(spks))=[];

    spikes.trial(unique(spks))=[];

    spikes.stimulus(unique(spks))=[];

    spikes.igorstim(unique(spks))=[];

    spikes.swtimes= spikes.fstimes/1000;
    spikes.spiketimes= spikes.ftimes/1000;
    spikes.noiseidx=[spikes.noiseidx; unique(spks)];

    % plot the new, noise-free, data
    close(NF);
end

NF = figure; hold on;
plot (spikes.waveforms');
axis tight; grid on
plot (mean (swaveforms), 'y', 'Linewidth', 4);
ylabel ('Voltage (uV)');
xlabel ('Sample (32 pts/ms)')
text (0.8*size(spikes.waveforms,2), 0.8*max(max(spikes.waveforms)), ['#Sp: ' num2str(length(stimes))]);


clear spks dpts z t opt


