
function sss_finalplot(spikes, useassigns, show, plot_ampl)
% figure; hold on; sss_finalplot(spikes, spikes.hierarchy.assigns, [1,2,15], 0.25);
%    OR CALL FROM Gsortm

if (nargin < 2)
    useassigns = spikes.hierarchy.assigns;
end
if (nargin < 3)
    show = unique(useassigns);
end
if (nargin < 4)
    plot_ampl = 300;
end

if isfield(spikes,'nsecs_reduced'),
    nsecs=spikes.nsecs_reduced;
elseif isfield(spikes,'nsecs'),
    nsecs=spikes.nsecs;
end

show = reshape(show, 1, []);
show(find(show==0))=[];

clf;
ylims = [min(spikes.waveforms(:)) max(spikes.waveforms(:))];
figure; hold on; 
redimscreen
%redimscreen_halfvert(0);

%set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold'); set(0,'defaultaxeslinewidth',2);
for clust = 1:length(show)
    
    
    members = find(useassigns == show(clust));
    memberwaves = spikes.waveforms(members,:);

    membertimes = sort(spikes.fstimes(members));

    tmin=1;
    tref = max(2, tmin*1.5);
    %isis = sort(diff(membertimes));
    isis = diff(membertimes);
    nrcsp = length(find ((isis <= tmin) & (isis>0))); % number of coincident spikes

    %figure(30); hold on; 
    subplot(length(show),4, 4 * (clust-1) + 1);
    [n,x,y] = hist2d(memberwaves);
    imagesc(x,y,n); axis xy; colormap hot;
    if (clust < length(show))
        xlabel([]);
    end
    %set(gca, 'YLim', ylims);
    set(gca, 'YLim', [-plot_ampl plot_ampl]);
    
    if clust==1
        msgstr = {};
        msgstr{end+1} = 'WAVEFORM';
       	msgstr{end+1} = ['NSpikes = ' num2str(length(membertimes)) ];   
        title(msgstr)
    else
        title(['NSpikes = ' num2str(length(membertimes)) ] );
    end
    
    if (show(clust) ~= 0)
        ylabel(['Cluster# ' num2str(show(clust))]);
    else
        ylabel('Outliers');
    end
    
     set(gca,'xtick',[1:2:7]*16, 'xticklabel',{['E1';'E2';'E3';'E4']});

%----------------------------------------------------
     
    a=50; b=200; c=1000;

    subplot(length(show),4,4 * (clust-1) + 2);
    [ta, scores] = isiQuality(membertimes/1000, membertimes/1000, tmin/1000, 0.010, tref/1000, spikes.Fs);
    aisis = isis(isis < a);

    if length(aisis)==0, aisis=a-1; end
    n = histc (aisis,0:0.2:a);
    bar(0:0.2:a,n,'histc');
    set(gca,'Xlim',[0 a]);
    if (clust == length(show))
        xlabel('ISI (ms)');
    end
    
    set(gca,'Xtick',[0:10:a],'Xticklabel',{num2str([0:10:a]')});
    if clust==1
        msgstr = {};
        msgstr{end+1} = 'HISTOGRAM';
        if exist('nsecs'),
            fr=length(membertimes)/nsecs;
            msgstr{end+1} = [' Firing Rate = ' num2str(roundn (fr,-1)) ' Hz' ];
        else
            msgstr{end+1} = [];
        end
        title(msgstr)
    else
         if exist('nsecs'),
            fr=length(membertimes)/nsecs;
            title([' Firing Rate = ' num2str(roundn (fr,-1)) ' Hz' ]);
        end
    end
    
%     if isfield(spikes,'nsecs'), fr=length(membertimes)/nsecs;
%         title([' Firing Rate = ' num2str(roundn (fr,-1)) ' Hz' ]);
%     else
%         title(['ISI score = ' num2str( roundn(scores(1),-2) )]);
%     end
    
    
    clear isis n

    subplot(length(show),4,4 * (clust-1) + 3); hold on;
    %%%%% GENERATE TS OBJECT %%%%%
    tMclust = membertimes*10;  % timestamps in 0.1ms resolution
    tsMcl = ts(tMclust,ts);
    acorr_bin_msec = 1; scale=30; %ms
    [histvals,x] = AutoCorr(Data(tsMcl), acorr_bin_msec, scale); % This is peter's C AutoCorr
    
    px=[flipud(-x);x]; phistvals=[flipud(histvals); histvals];
    bar(px,phistvals,'FaceColor','k');			% show acorr
    m=1.2*max(histvals);  arr=0:0.1:m;
    plot(1*ones(1,length(arr)),arr,'r-','Markersize',2,'Linewidth',2);
    plot(-1*ones(1,length(arr)),arr,'r-','Markersize',2,'Linewidth',2);

    axis([-scale scale 0 1.3*m]);   
    if (clust == length(show))
        xlabel(['ms (' num2str(acorr_bin_msec) 'ms binsize)']);
        ylabel('Rate');
    end
    if clust==1
        msgstr = {};
        msgstr{end+1} = 'AUTOCORR ';
       	msgstr{end+1} = ['Coinspks= ' num2str(nrcsp) '; %= ' num2str( roundn( 100*nrcsp/length(membertimes),-2))];   
        title(msgstr)
    else
        title(['Coinspks= ' num2str(nrcsp) '; %= ' num2str( roundn( 100*nrcsp/length(membertimes),-2))] );
    end
    drawnow;


    subplot(length(show),4,4 * (clust-1) + 4);
    %%%%% GENERATE TS OBJECT %%%%%
    tMclust = membertimes*10;  % timestamps in 0.1ms resolution
    tsMcl = ts(tMclust,ts);

    HistISI(tsMcl);
    %set(gca,'FontSize',text_font_size);
    
    if clust==1
        msgstr = {};
        msgstr{end+1} = 'HIST OF LOG(ISI)';
       	msgstr{end+1} = ['ISI score = ' num2str( roundn(scores(1),-2) )];   
        title(msgstr)
    else
        title(['ISI score = ' num2str( roundn(scores(1),-2) )]);
    end
    
    
    if (clust == length(show))
        ylabel('nSpikes','FontSize',18,'Fontweight','bold');
        xlabel('ISI, ms','FontSize',18,'Fontweight','bold');
    end
        drawnow

end
