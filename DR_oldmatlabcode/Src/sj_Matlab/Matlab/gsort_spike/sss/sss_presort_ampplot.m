function sss_presort_ampplot(spikes, nch, plotpca)
% Amp plots AND PC PLOTS for diff channels of ntrode: before Sorting
% If ntrode>=4: shows six subplots for first 4 channels
%  Shantanu - Apr 2011

if (nargin < 2), nch=4; end
if (nargin < 3), plotpca=0; end

%%numclusts = length(clusters);
%numclusts = length(show);

if plotpca==1
    if ~isfield(spikes,'pcadata')
        disp('PCAs do not exist.Skipping. Please hit "Do PCA" first');
        plotpca=0;
    end
end


%% Amplitude Plots

channel_lth = size(spikes.waveforms,2)/nch;
nSpikes = size(spikes.waveforms,1);

% Calculate Amplitude for mulit-channel data
if ~isfield(spikes,'ampl')
    
    for ch=1:nch
        if isfield(spikes,'waveforms_ch1')
            cmd=sprintf('ampl(:,ch)=abs(max(spikes.waveforms_ch%d,[],2)) + abs(min(spikes.waveforms_ch%d,[],2));',ch,ch); eval(cmd);
        else
            w = spikes.waveforms(:,channel_lth*(ch-1)+1:channel_lth*ch);
            ampl(:,ch) = abs(max(w,[],2)) + abs(min(w,[],2));
        end
    end    
else
    ampl=spikes.ampl;
end




if nch==1
    
    disp('Single Channel Only - Plotting Amplitude against time')  
    spktimes = spikes.fstimes./1000; % in sec
    
    afig=figure; hold on;
    %    redimscreen; orient(gcf,'landscape');
    set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
    ylabel('Amplitude (uV)' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Time (sec)','FontSize', 14, 'FontWeight', 'bold');
    
    hndl = plot(spktimes, ampl, '.');
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    title('Amplitude vs Time')
    hold off;
    
end



if nch==2
       
    afig=figure; hold on;
    %    redimscreen; orient(gcf,'landscape');
    set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
    ylabel('Electrode 2' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Electrode 1','FontSize', 14, 'FontWeight', 'bold');
    
    hndl = plot(ampl(:,1), ampl(:,2), '.');
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    title('Amp Plot for Stereotrode')
    hold off;
    
    set(afig, 'ButtonDownFcn', {@make_density, ampl});
    figure(afig);
    
    return
end

% Make Ampl_plot

if nch>=3
    
    for ch=1:nch     %n channels for ntrode
        if isfield(spikes,'waveforms_ch1')
            cmd=sprintf('ampl(:,ch)=abs(max(spikes.waveforms_ch%d,[],2)) + abs(min(spikes.waveforms_ch%d,[],2));',ch,ch); eval(cmd);
        else
            w = spikes.waveforms(:,channel_lth*(ch-1)+1:channel_lth*ch);
            ampl(:,ch) = abs(max(w,[],2)) + abs(min(w,[],2));
        end
    end
    
    % Make the plot.
    figure;  hfig = gcf;  hax = gca;    hold on; redimscreen100(100);
    %    redimscreen; orient(gcf,'landscape');
    set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
    
    %whitebg(hfig);
    
    subplot(3,2,1); hold on;
    hndl = plot(ampl(:,1), ampl(:,2), '.','Markersize',2);
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    
    title('AMP PLOTS','FontSize', 14, 'FontWeight', 'bold');
    ylabel('Electrode 2' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Electrode 1','FontSize', 14, 'FontWeight', 'bold');
    
    subplot(3,2,2); hold on;
    ylabel('Electrode 3' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Electrode 2','FontSize', 14, 'FontWeight', 'bold');
    hndl = plot(ampl(:,2), ampl(:,3), '.','Markersize',2);
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    
    
    subplot(3,2,3); hold on;
    ylabel('Electrode 2' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Electrode 1','FontSize', 14, 'FontWeight', 'bold');
    hndl = plot(ampl(:,1), ampl(:,2), '.','Markersize',2);
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    
    if nch>=4
        subplot(3,2,4); hold on;
        ylabel('Electrode 4' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Electrode 2','FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(ampl(:,2), ampl(:,4), '.','Markersize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        
        
        subplot(3,2,5); hold on;
        ylabel('Electrode 4' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Electrode 1','FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(ampl(:,1), ampl(:,4), '.','Markersize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        
        
        subplot(3,2,6); hold on;
        ylabel('Electrode 4' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Electrode 3','FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(ampl(:,3), ampl(:,4), '.','Markersize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    end
    
end

spikes.ampl=ampl;

if plotpca==1,
    
    
    if nch==1
        
        disp('Single Channel Only - Plotting PCA against time')
        spktimes = spikes.fstimes./1000; % in sec
        
        bfig=figure; hold on;
        %    redimscreen; orient(gcf,'landscape');
        set(gcf, 'PaperPositionMode', 'auto');
        set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
        ylabel('PCA1' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Time (sec)','FontSize', 14, 'FontWeight', 'bold');
        
        hndl = plot(spktimes, spikes.pcadata.pc1(:,1), '.');
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        title('PCA1 vs Time')
        hold off;
        
        return
        
    end
    
    if nch==2
        
        afig=figure; hold on;
        %    redimscreen; orient(gcf,'landscape');
        set(gcf, 'PaperPositionMode', 'auto');
        set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
        ylabel('PCA1: Ch2' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('PCA1: Ch1','FontSize', 14, 'FontWeight', 'bold');
        
        hndl = plot(spikes.pcadata.pc1(:,1), spikes.pcadata.pc1(:,2), '.');
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        title('PCA Plot for Stereotrode')
        hold off;
        
        
        return
    end
    
    
    if nch>=3
        
        % Make the PCA plot.
        figure;  hfig = gcf;     hax = gca;    hold on; redimscreen100(100);
        %whitebg(hfig);
        %    redimscreen; orient(gcf,'landscape');
        set(gcf, 'PaperPositionMode', 'auto');
        set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
        
        
        subplot(3,2,1); hold on;
        ylabel('PCA1: Ch2' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('PCA1: Ch1' ,'FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(spikes.pcadata.pc1(:,1), spikes.pcadata.pc1(:,2), '.','Markersize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        title('PCA PLOTS')
        
        subplot(3,2,2); hold on;
        ylabel('PCA1: Ch3' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('PCA1: Ch2','FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(spikes.pcadata.pc1(:,2), spikes.pcadata.pc1(:,3), '.','Markersize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        
        
        subplot(3,2,3); hold on;
        ylabel('PCA1: Ch3' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('PCA1: Ch1','FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(spikes.pcadata.pc1(:,1), spikes.pcadata.pc1(:,3), '.','Markersize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        
        if nch>=4
            subplot(3,2,4); hold on;
            ylabel('PCA1: Ch4' ,'FontSize', 14, 'FontWeight', 'bold');
            xlabel('PCA1: Ch2','FontSize', 14, 'FontWeight', 'bold');
            hndl = plot(spikes.pcadata.pc1(:,2), spikes.pcadata.pc1(:,4), '.','Markersize',2);
            %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
            
            
            subplot(3,2,5); hold on;
            ylabel('PCA1: Ch4' ,'FontSize', 14, 'FontWeight', 'bold');
            xlabel('PCA1: Ch1','FontSize', 14, 'FontWeight', 'bold');
            hndl = plot(spikes.pcadata.pc1(:,1), spikes.pcadata.pc1(:,4), '.','Markersize',2);
            %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
            
            
            subplot(3,2,6); hold on;
            ylabel('PCA1: Ch4' ,'FontSize', 14, 'FontWeight', 'bold');
            xlabel('PCA1: Ch3','FontSize', 14, 'FontWeight', 'bold');
            hndl = plot(spikes.pcadata.pc1(:,3), spikes.pcadata.pc1(:,4), '.','Markersize',2);
            %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        end
        
    end
end

%%%%%%%%%%% Density Plot for Stereotrode Amplitude  %%%%%%%%%%%%%%%%

function make_density(afig,event,ampl)
figure; hdens = gca;
histxy(ampl(:,1),ampl(:,2),250,1);
%set(gca,properties2d{:},'Color',[0 0 0.508]);
HGlogmenu(findobj(gca,'Type','Image'));
