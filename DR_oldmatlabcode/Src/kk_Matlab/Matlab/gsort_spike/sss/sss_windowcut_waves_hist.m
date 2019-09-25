function [waves, timestamps, noiseidx] = sss_windowcut_waves_hist(waves, timestamps, n_ch)

noiseidx=[];
amplitude=400;
% Tlines={'How Many Channels?',['Max Amplitude on Plot?']};
% defent={'4','400'}; % default entries
% infonames= {'n_channels','amplitude'};
% info = inputdlg(Tlines, 'How many channels? And Amplitude?', 1, defent);


% main = figure; hold on;
% [n,x,y] = hist2d(spikes.waveforms);
% imagesc(x,y,n); axis xy; colormap hot;
% axis tight; grid on
% plot (mean (spikes.waveforms), 'y', 'Linewidth', 2);
% ylabel ('Voltage (uV)');
% xlabel ('Sample (32 pts/ms)')
% text (0.8*size(spikes.waveforms,2), 0.8*max(max(spikes.waveforms)), ['#Sp: ' num2str(length(stimes))]);

main=figure; colormap hot; hold on;
if (size (waves,3) > 1000)
    idxs = randperm (size (waves,3));
    unit_sm = waves(:,:,idxs(1:1000));
else
    unit_sm=waves;
end
Nsamp = size(unit_sm,1);
Nch = size(unit_sm,2);
Nspk = size(unit_sm,3);
unit_sm = reshape(unit_sm,Nsamp*Nch,Nspk);
plot(unit_sm,'b.-'); hold on;
plot (mean (unit_sm,2), 'y', 'Linewidth', 2);
text (0.8*size(unit_sm,1), double(0.8*max(max(unit_sm))), ['#Sp: ' num2str(size(waves,3))]);

clear unit_sm

% Tlines={'How Many Channels?'};
% defent={'4'}; % default entries
% infonames= {'n_channels'};
% info = inputdlg(Tlines, 'How many channels?' , 1, defent);
%
% if ~isempty(info)              %see if user hit cancel
%     info = cell2struct(info,infonames);
%     n_ch = str2num(info.n_channels);   %convert string to number
%     %amplitude = str2num(info.amplitude);   %convert string to number
% else
%     return
% end

for i=1:n_ch
    Tlines={'Amplitude?'};
    defent={'500'}; % default entries
    infonames= {'amplitude'};
    info = inputdlg(Tlines, 'What Amplitude to use?' , 1, defent);
    
    if ~isempty(info)              %see if user hit cancel
        info = cell2struct(info,infonames);
        %n_ch = str2num(info.n_channels);   %convert string to number
        amplitude = str2num(info.amplitude);   %convert string to number
    end
    
    % Plot the waveforms and then ask user what should be the windowcut limits
    cmd=sprintf('swaveforms = double(squeeze(waves(:,%d,:))\'');',i); eval(cmd);
    NF = figure; hold on;
    [n,x,y] = hist2d(swaveforms);
    imagesc(x,y,n); axis xy; colormap hot;
    grid on;
    axis([0 size(swaveforms,2) -amplitude amplitude]);
    plot (mean (swaveforms), 'y', 'Linewidth', 2)
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
        waves(:,:,unique(spks))=[];
        timestamps (unique(spks))=[];
    end
    
    if isempty(timestamps);
        disp ('WARNING! There is no spike within the window specified. Please rerun the program with new values')
        return
    end
    noiseidx=[noiseidx; unique(spks)];
    
    close(NF);
    
    
    %%Replot main
    close(main)
    main=figure; colormap hot; hold on;
    if (size (waves,3) > 1000)
        idxs = randperm (size (waves,3));
        unit_sm = waves(:,:,idxs(1:1000));
    else
        unit_sm=waves;
    end
    Nsamp = size(unit_sm,1);
    Nch = size(unit_sm,2);
    Nspk = size(unit_sm,3);
    unit_sm = reshape(unit_sm,Nsamp*Nch,Nspk);
    plot(unit_sm,'b.-'); hold on;
    plot (mean (unit_sm,2), 'y', 'Linewidth', 2);
    text (0.8*size(unit_sm,1), double(0.8*max(max(unit_sm))), ['#Sp: ' num2str(size(waves,3))]);
    
end

close(main)
% NF = figure; hold on;
% [n,x,y] = hist2d(spikes.waveforms);
% imagesc(x,y,n); axis xy; colormap hot;
% grid on;
% ylabel ('Voltage (uV)');
% xlabel ('Sample (32 pts/ms)')
% text (0.8*size(spikes.waveforms,2), 0.8*max(max(spikes.waveforms)), ['#Sp: ' num2str(length(stimes))]);


clear spks dpts z t opt


