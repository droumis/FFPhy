function spikes = sss_windowcut(spikes)

swaveforms = spikes.waveforms;
stimes = spikes.fstimes';



    % Plot the waveforms and then ask user what should be the windowcut limits
    
    NF = figure; hold on;
    plot (swaveforms'); grid on; 
    axis([0 size(swaveforms,2) -500 500]); 
    plot (mean (swaveforms), 'y', 'Linewidth', 4)
    redimscreen
    
    disp ('A graph of the raw waveforms are being plotted please consult the figure and return to the command window')
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

if isempty(swaveforms);
    disp ('WARNING! There is no spike within the window specified. Please rerun the program with new values')
    %return
end


if (length(spks)~=0)
    swaveforms (unique(spks),:)=[];
    stimes (unique(spks))=[];

    if(exist('spikes.waveforms_ch2'))
        spikes.waveforms_ch1(unique(spks),:)=[];
        spikes.waveforms_ch2(unique(spks),:)=[];
        spikes.waveforms_ch3(unique(spks),:)=[];
        spikes.waveforms_ch4(unique(spks),:)=[];
    end

    spikes.ftimes(unique(spks))=[];

    if(exist('spikes.sweep'))
        spikes.sweep(unique(spks))=[];
    end
    
    if(exist('spikes.trial'))
        spikes.trial(unique(spks))=[];
    end
    
    if(exist('spikes.stimulus'))
        spikes.stimulus(unique(spks))=[];
    end
    
    if(exist('spikes.igorstim'))
       spikes.igorstim(unique(spks))=[];
    end
    
end

%%% UPDATE WITH NEW %%%
spikes.waveforms=swaveforms;
spikes.fstimes=stimes';
%spikes.fstimes(unique(spks))=[];
spikes.swtimes= spikes.fstimes/1000;
spikes.spiketimes= spikes.ftimes/1000;


spikes.noiseidx=unique(spks);

%plot the new, noise-free, data
close(NF);
NF = figure; hold on;
plot (swaveforms');
axis tight; grid on
plot (mean (swaveforms), 'y', 'Linewidth', 4);
ylabel ('Voltage (uV)');
xlabel ('Sample (32 pts/ms)')
text (0.8*size(swaveforms,2), 0.8*max(max(swaveforms)), ['#Sp: ' num2str(length(stimes))]);


clear spks dpts z t opt 

