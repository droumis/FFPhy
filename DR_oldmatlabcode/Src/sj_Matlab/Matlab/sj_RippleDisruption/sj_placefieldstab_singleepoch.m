function [allcorr,all2dcorr] = sj_placefieldstab_singleepoch (animdirect,prefix,days,allepoch,type,tetlist,cellarray,includeStates, minvel, binsizelin, binsize2d, figopt1,saveg1,savedata)

%% Generate and Plot Place fields for cells on tetrode
% sj_placefieldstab_singleepoch('/data25/sjadhav/RippleInterruption/sjc_direct','sjc',4,[2],'Run',{[5]}, {{[2:3 9:10]}}, 3, 6, 2, 1, 1,0,0);

%%

if nargin<8,
    includeStates = 3; % for getbehavestate
end
if nargin<9,
    minvel = 6; % for getbehavestate
end

if nargin<10,
    binsizelin = 2; % for calclinfields
end
if nargin<11,
    binsize2d = 1; % for openfieldrate
end

if nargin<12,
    figopt1 = 0; % Plot figures for each cell
end
if nargin<13,
    saveg1 = 0;
end

if nargin<14,
    savedata = 0;
end

%% Fixed parameters
directoryname = animdirect;
clr = {'b','r','g','m','c','y','k','r'};


%% Loop over days
alldays_stim_spkshist_mat =[];
alldays_rip_spkshist_mat =[];

cntneu=0;  % Count neurons across days
cntneu_ep2=0;
for d = 1:length(days)
    
    day = days(d)
    
    % Load pos file
    posfile = sprintf('%s/%spos%02d.mat', directoryname, prefix, day);
    load(posfile);
    
    % Load linpos file
    linposfile = sprintf('%s/%slinpos%02d.mat', directoryname, prefix, day);
    load(linposfile);
    
    % Load spike file
    spkfile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
    load(spkfile);
    
    allepoch
    
    statefile = sprintf('%s/%sstate%02d.mat', directoryname, prefix, day);
    if ~exist(statefile,'file');
        [state, lindist] = getbehavestate(linpos, day, allepoch, includeStates, 'minlinvelocity', minvel);
        behavestate{day}{allepoch}.state=state;
        behavestate{day}{allepoch}.lindist=lindist;
        save(statefile,'behavestate');
    else
        load(statefile);
        state = behavestate{day}{allepoch}.state;
        lindist = behavestate{day}{allepoch}.lindist;
    end
    
    
    % Divide time in current epoch
    timepts = pos{day}{allepoch}.data(:,1);
    sep_epoch_idx=floor(length(timepts)/2);
    idx{1}=1:sep_epoch_idx;
    idx{2}=sep_epoch_idx+1:min([length(state),length(timepts)]);
    
    for ep=1:2,
        
        % Split epoch
        curridxs = idx{ep}; 
        currtime = timepts(curridxs);
        % Behave state
            %currstate =  state(curridxs);
            %currlindist = lindist(curridxs);
        currstate=state;
        currlindist=lindist;
        % Pos
        currpos=pos{day}{allepoch};
        currpos.data=currpos.data(curridxs,:);
        
        tet = tetlist{d}; % tetrodes for current day
        celllist = cellarray{d}; % cell list for terodes of current day
        
        for elec = 1:length(tet)      
            currtet = tet(elec);
            cells = celllist{elec};
            
            for neuron = 1:length(cells)
                currcell = cells(neuron);
                
                 % Spikes
                tmpspikes = spikes{day}{allepoch}{currtet}{currcell};
                      %Can use time
                %tmpspkidxs = find( (tmpspikes.data(:,1)>=currtime(1)) & (tmpspikes.data(:,1)<=currtime(end)) );
                %tmpspikes.data = currspikes.data(currspkidxs,:);
                      % or the position index field
                tmpspkidxs = find( (tmpspikes.data(:,7)>=curridxs(1)) & (tmpspikes.data(:,7)<=curridxs(end)) );
                tmpspikes.data = tmpspikes.data(tmpspkidxs,:);
                currspikes{day}{allepoch}{currtet}{currcell}=tmpspikes;
                 
                % Linearized traj
                trajdata = calclinfields(currspikes,currstate,currlindist,linpos, [day allepoch currtet currcell],binsizelin);
                %cmd = sprintf('trajdata_day%02d_ep%02d_tet%02d_cell%02d=trajdata;',day,allepoch,currtet,currcell); eval(cmd);
                
                % 2d map
                [placemap] = openfieldrate(tmpspikes.data,currpos.data,binsize2d,2);
                %cmd = sprintf('placemap_day%02d_ep%02d_tet%02d_cell%02d=placemap;',day,allepoch,currtet,currcell); eval(cmd);
                
                % Get traj with maximum peak fir rate
                for i=1:length(trajdata),
                    maxrate(i) = max(trajdata{i}(5:end-5,5));
                end
                
                if figopt1==1,
                    
                    figure(str2num([num2str(day) num2str(currtet) num2str(currcell)]));
                    redimscreen;
                    hold on;
                    
                    subplot(2,2,ep); hold on;
                    for i=1:length(trajdata),
                        plot(trajdata{i}(:,5),[clr{i} '.-'],'Linewidth',2);
                        %maxrate(i) = max(trajdata{i}(5:end-5,5));
                    end
                    if ep==1
                        title(['Day ' num2str(day) '  Epoch ' num2str(allepoch) '- Part 1'  'Tet ' num2str(currtet) '  Cell ' num2str(currcell)],'FontSize',14,'Fontweight','bold');
                        legend('InRight','OutRight','InLeft','OutLeft');
                    else
                        title(['Epoch ' num2str(allepoch) '- Part 2'],'FontSize',14,'Fontweight','bold');
                    end
                    
                    subplot(2,2,2+ep); hold on;
                    imagesc(placemap.smoothedspikerate); colorbar
                    title(['Epoch ' num2str(allepoch) '- Part ' num2str(ep)], 'FontSize',14,'Fontweight','bold');
                    
                end
                
                if ep==1
                    cntneu=cntneu+1;
                    
                    allplacemaps_ep1{cntneu} = placemap;
                    alllintrajs_ep1{cntneu} = trajdata;
                    allpeakrates_ep1{cntneu} = maxrate;
                end
                
                if ep==2
                    cntneu_ep2=cntneu_ep2+1;
                    allplacemaps_ep2{cntneu_ep2} = placemap;
                    alllintrajs_ep2{cntneu_ep2} = trajdata;
                    allpeakrates_ep2{cntneu_ep2} = maxrate;
                end
                
                trajplace{day}{ep}{currtet}{currcell}=trajdata;
                mapplace{day}{ep}{currtet}{currcell}=placemap;
                
                clear trajdata placemap maxrate
                clear tmpspikes currspikes
                
            end % end neuron
        end % end elec
        
    end % end ep
    
    
    %dayplacefieldfile = sprintf('%s/%placefield%02d.mat', directoryname, prefix, day);
    %save(dayplacefieldfile,'trajplace','mapplace');
    clear trajplace mapplace
    
end % end day




% Find which trajectory had the max firing rate in all epochs
for n=1:cntneu
    
    peakrate=[]; peaktraj=[];
    
    [peakrate(1),peaktraj1] = max(allpeakrates_ep1{n});
    [peakrate(2),peaktraj2] = max(allpeakrates_ep2{n});
    
    
    if peaktraj1(1)==peaktraj2(1),
        allmaxtraj{n}=peaktraj1(1);
    else
        allmaxtraj{n}=[peaktraj1(1) peaktraj2(1)];
    end
    
end

% Get corrleations between place fields

single_traj=0; double_traj=0;
for n=1:cntneu
    n
    
    % 2d corr
    sizex = min([length(allplacemaps_ep1{n}.yticks),length(allplacemaps_ep2{n}.yticks)]);
    sizey = min([length(allplacemaps_ep1{n}.xticks),length(allplacemaps_ep2{n}.xticks)]);
    corr2d = corr2(allplacemaps_ep1{n}.smoothedspikerate(1:sizex,1:sizey), allplacemaps_ep2{n}.smoothedspikerate(1:sizex,1:sizey));
    
    
    % Lin Traj Corr
    maxtraj = allmaxtraj{n};
    if length(maxtraj)==1
        single_traj = single_traj+1;
        place_ep1 = alllintrajs_ep1{n}{maxtraj}(:,5);
        place_ep2 = alllintrajs_ep2{n}{maxtraj}(:,5);
        
        invalid1 = find(isnan(place_ep1)==1);
        invalid2 = find(isnan(place_ep2)==1);
        invalid=unique([invalid1; invalid2]);
        place_ep1(invalid)=[];
        place_ep2(invalid)=[];
        
        currcorr = corrcoef(place_ep1,place_ep2);
        currcorr = currcorr(1,2);
    else
        double_traj = double_traj+1;
        place1_ep1 = alllintrajs_ep1{n}{maxtraj(1)}(:,5);
        place1_ep2 = alllintrajs_ep2{n}{maxtraj(1)}(:,5);
        
        invalid1 = find(isnan(place1_ep1)==1);
        invalid2 = find(isnan(place1_ep2)==1);
        invalid=unique([invalid1; invalid2]);
        place1_ep1(invalid)=[];
        place1_ep2(invalid)=[];
        currcorr1 = corrcoef(place1_ep1,place1_ep2);
        currcorr1 = currcorr1(1,2);
        
        place2_ep1 = alllintrajs_ep1{n}{maxtraj(2)}(:,5);
        place2_ep2 = alllintrajs_ep2{n}{maxtraj(2)}(:,5);
        invalid1 = find(isnan(place2_ep1)==1);
        invalid2 = find(isnan(place2_ep2)==1);
        invalid=unique([invalid1; invalid2]);
        place2_ep1(invalid)=[];
        place2_ep2(invalid)=[];
        currcorr2 = corrcoef(place2_ep1,place2_ep2);
        currcorr2 = currcorr2(1,2);
        
        currcorr = mean([currcorr1,currcorr2]);
        
        %         place_ep1=sum(place1_ep1,place2_ep1);
        %         place_ep2=sum(place1_ep2,place2_ep2);
        %         currcorr = corrcoef(place_ep1,place_ep2);
        %         currcorr = currcorr(1,2)
        
    end
    
    allcorr(n)=currcorr;
    all2dcorr(n)=corr2d;
    
end

disp('Mean LinTrajCorr');
nanmean(allcorr),
nanstd(allcorr),

disp('Mean 2dCorr');
nanmean(all2dcorr),
nanstd(all2dcorr),



% % Plot for day
% if figopt1==1
%
%     figure(day); hold on;
%     redimscreen_2versubplots;
%     orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%
%     subplot(2,1,1); hold on;
%     yplot = (1000/binsize_plot)*mean(day_stim_spkshist,1); %Multiunit fr (Hz) in binsize ms bins
%     taxis = [0:binsize:600];
%     taxis = taxis - 200;
%     plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
%     bar(taxis, yplot,'r');
%     %plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
%
%     ypts = 0:1.1*max(yplot); xpts = 0*ones(size(ypts));
%     plot(xpts , ypts, 'k--','Linewidth',2);
%     title(['MUStim-' type '-Day' num2str(day)],'FontSize',20,'Fontweight','bold');
%     ylabel('Instantaeous MU Firing Rate'); %xlabel('Time(ms)');
%
%     subplot(2,1,2); hold on;
%     yplot = (1000/binsize_plot)*mean(day_rip_spkshist,1); %Multiunit fr (Hz) in binsize ms bins
%     taxis = [0:binsize:600];
%     taxis = taxis - 200;
%     plot(taxis, yplot,['r.-'],'Linewidth',2,'Markersize',12);
%     bar(taxis, yplot,'r');
%     %plot(taxis, 2*ripnostim_pre(i,:),['r-'],'Linewidth',2,'Markersize',6);
%
%     ypts = 0:1.1*max(yplot); xpts = 0*ones(size(ypts));
%     plot(xpts , ypts, 'k--','Linewidth',2);
%     title(['MURip'],'FontSize',20,'Fontweight','bold');
%     ylabel('Instantaeous MU Firing Rate');
%     xlabel('Time(ms)');
%
%     if saveg1==1,
%         orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%         saveas(gcf,[prefix 'MultiFRstimrip_' type '_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'fig');
%         saveas(gcf,[prefix 'MultiFrstimrip_' type '_' num2str(day) '_tet' num2str(tet) '_SEPstd' num2str(std)],'jpg');
%     end
%
% end  % end figure
%
% end % end day
%
%
% %%%%%% Plot Across Days  %%%%%%%%
%
%
%
%
% if savedata==1
%     savefile = sprintf('%s/ProcessedData/%sMultiFRstimrip_Days%01dto%01d_Tet%02d.mat', animdirect, prefix, min(days), max(days), tet);
%     save(savefile);
% end







