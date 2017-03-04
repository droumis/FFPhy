warning('off', 'MATLAB:Figure:RecursionOnClose') %suppress this specific warning
%% gather the ripple window data from all the regions
for currrip=1:length(ripsStartTimes{srcRegionind}) %for each ripples from source region..
    skipripscount = 0;
    clear YripLFPdata
    for i = 1:length(regions); %for each region
        skiprip = 0;
        tetcnt = 0;
        for tet = selected_tets{i}
            tetcnt = tetcnt + 1;
            if strcmp(LFPtype,'eeg')
                YripLFPdata{i}(:,tetcnt) = eegstruct{i}{day}{epoch}{tet}.data(ripsStartTimesind{srcRegionind}(currrip)-(windowsize*samprate):ripsStartTimesind{srcRegionind}(currrip)+(windowsize*samprate));
            elseif strcmp(LFPtype,'ripple')
                YripLFPdata{i}(:,tetcnt) = eegstruct{i}{day}{epoch}{tet}.data(ripsStartTimesind{srcRegionind}(currrip)-(windowsize*samprate):ripsStartTimesind{srcRegionind}(currrip)+(windowsize*samprate),2); %use second col for rip data
            else
                disp('incompatible LFPtype')
                return
            end
            lfptraceregion{i} = ones(length(YripLFPdata{i}(1,:)),1)*i;
        end
    end
    YripLFPdataMAT = cell2mat(YripLFPdata); %stack all the traces next to each other
    Xwindowtimes = WindowStartEndTimes(currrip,1):1/samprate:WindowStartEndTimes(currrip,2);
    % lfptraceLUTregion = cell2mat(cellfun(@(x) repmat(size(x,2),1,size(x,2)),YripLFPdata,'UniformOutput',false)); %get the number of LFP traces for each region
    lfptraceLUTregion = cell2mat(lfptraceregion');
    %% plot params
    fcclr = ([1 1 1 ; 0.2 .6 .9; 0 0 1; 0 1 0]);
    edgclr = ([170,11,10; 80,215,167; 255,170,55;0,40,185])./255;
    patchclr = (1.-edgclr)/1.2 + edgclr;   % ([.9 .7 .7 ; 0.2 .6 .9; 0 0 1; 0 1 0]);
    % edgclr = ([1 0 0; 0.2 .6 .9; 0 0 1; 0 1 0]);
    lintyp = (['o' '^' '^' '^']);
    linwid = ([2 2 2 2]);
    marksz = ([8 7 7 7]);
    clear ripsinwin
    % ha = tight_subplot(length(regions),1,[.01 .03],[.1 .1],[.1 .05]);
    %% plot LFP
    close(fig)
    fig = figure;
    for p = 1:length(YripLFPdataMAT(1,:)); %for each lfp trace
        if p == 1;
            plot(Xwindowtimes,YripLFPdataMAT(:,p), 'Color',edgclr(lfptraceLUTregion(p),:), 'LineWidth',2); hold on;
            lfpoffset = 0;
        else
            lfpoffset(p) = lfpoffset(p-1) + abs(min(YripLFPdataMAT(:,p-1))) + abs(max(YripLFPdataMAT(:,p))); %use abs(max(YripLFPdataMAT(:,p)))/3 to overlay traces more
            plot(Xwindowtimes,YripLFPdataMAT(:,p) - lfpoffset(p), 'Color',edgclr(lfptraceLUTregion(p),:), 'LineWidth',2); hold on;
        end
    end
    %% Plot Ripple times and accessory stuff
    for i = 1:length(regions); %for each region
        ripsinwinInds = (ripsStartTimes{i}>WindowStartEndTimes(currrip,1) & ripsStartTimes{i}<WindowStartEndTimes(currrip,2));
        ripsinwin{currrip}{i} = [ripsStartTimes{i}(ripsinwinInds) ripsEndTimes{i}(ripsinwinInds)];
        if ~isempty(ripsinwin{currrip}{i}) %if there are any ripples from this region in this window
            for m = 1:length(ripsinwin{currrip}{i}(:,1)) %for each ripple within the current window
                Ylfpranges4region{i} = [-lfpoffset(find(lfptraceLUTregion == i,1,'last'))-(abs(min(YripLFPdataMAT(:,find(lfptraceLUTregion == i,1,'last'))))) -lfpoffset(find(lfptraceLUTregion == i,1,'first'))+max(YripLFPdataMAT(:,find(lfptraceLUTregion == i,1,'first')))];
                %                 Ylfpranges4region{i} = [-lfpoffset(find(lfptraceLUTregion == i,1,'last'))-(abs(min(YripLFPdataMAT(find(lfptraceLUTregion == i,1,'last'))))) abs(max(YripLFPdataMAT(:,find(lfptraceLUTregion == i,1,'first'))))];
                line([ripsinwin{currrip}{i}(m,1) ripsinwin{currrip}{i}(m,1)], Ylfpranges4region{i},'Color',edgclr(i,:),'LineWidth',1.5);
                Xpatch = [ripsinwin{currrip}{i}(m,1) ripsinwin{currrip}{i}(m,2) ripsinwin{currrip}{i}(m,2) ripsinwin{currrip}{i}(m,1)];
                %                 Ypatch = [Ylfpranges4region{i}(1) Ylfpranges4region{i}(1)+diff(Ylfpranges4region{i})/4 Ylfpranges4region{i}(2)-diff(Ylfpranges4region{i})/4 Ylfpranges4region{i}(2)]; %trapezoidal patch that decays toward rip end
                Ypatch = [Ylfpranges4region{i}(1) Ylfpranges4region{i}(1) Ylfpranges4region{i}(2) Ylfpranges4region{i}(2)];
                patch(Xpatch, Ypatch, patchclr(i,:), 'edgecolor','none'); %triggering-ripple patch
            end
        end
    end
    centerripStartTime=ripsStartTimes{srcRegionind}(currrip);
    Yextent = [-lfpoffset(end)-(abs(min(YripLFPdataMAT(:,end)))) max(max(YripLFPdataMAT(:,1)))];
    line([centerripStartTime centerripStartTime], Yextent ,'Color',edgclr(srcRegionind,:), 'LineStyle', '--','LineWidth',1.5) %line for the center trigger-ripple
    set(gca,'children',flipud(get(gca,'children'))) %send the patch behind the LFP traces
    ylim([Yextent(1) Yextent(2)])
    xlim([WindowStartEndTimes(currrip,1) WindowStartEndTimes(currrip,2)])
    xl = xlim;
    for k=1:length(regions)
        text(xl(1)-diff(WindowStartEndTimes(currrip,:))/10, -lfpoffset(find(lfptraceLUTregion == k,1,'first')), regions{k}, 'Color', edgclr(k,:),'FontSize',13,'FontWeight','bold')
    end
    set(gca, 'YTick', []);
    set(gca,'XTick',[WindowStartEndTimes(currrip,1):diff(WindowStartEndTimes(currrip,:))/10:WindowStartEndTimes(currrip,2)], 'FontSize',10,'FontWeight','bold')
    set(gca, 'XTickLabel', [-windowsize:windowsize/5:windowsize])
    xlabel('seconds from rip start','FontSize',12,'FontWeight','bold','Color','k')
    title({[sprintf('%s d%de%d %sRip(%d) Triggered LFP',animal, day, epoch, rippleregion, currrip)];[sprintf('%16.f', centerripStartTime)]},'FontSize',12,'FontWeight','bold')
    pause
end