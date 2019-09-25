function [allcorr,all2dcorr] = sj_placefieldstab_epochs_1_onlymaps (animdirect,prefix,days,allepochs,type,includeStates, minvel, binsizelin, binsize2d, figopt1,saveg1,savedata)
%% Generate and Plot Place fields for cells on tetrode

% Shantanu: Change to get tet list and cell list automatically from day file
% 06/02/2011

% sj_placefieldstab_epochs_1_onlymaps('/data25/sjadhav/RippleInterruption/REf_direct','REf',3,[2 4],'Run', 3, 6, 2, 1, 1,0,0);

% With tetlist and cellarray fed in
% sj_placefieldstab_epochs_1_onlymaps('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',[3 7 8],[2 4],'Run',{[1 5 7],1,7}, {{[2:3],[2 3 5],[2:4]},{[2:4]},{1:2}}, 3, 6, 2, 1, 1,0,0);
% sj_placefieldstab_epochs_1_onlymaps('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',3,[2 4],'Run',{[1 5 7]}, {{[2:3],[2 3 5],[2:5]}}, 3, 6, 2, 1, 1,0,0);
% sj_placefieldstab_epochs_1_onlymaps('/data25/sjadhav/RippleInterruption/REd_direct','REd',3,[2 4],'Run',{[10]}, {{1}}, 3, 6, 2, 1, 1,0,0);
% sj_placefieldstab_epochs_1_onlymaps('/data25/sjadhav/RippleInterruption/REf_direct','REf',9,[2 4],'Run',{[9]}, {{1}}, 3, 6, 2, 1, 1,0,0);
% sj_placefieldstab_epochs_1('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',3,[1 3 5],'Sleep',5, {2:5}, 3,6,2,1, 1,0,0);
%

if nargin<6,
    includeStates = 3; % for getbehavestate
end
if nargin<7,
    minvel = 6; % for getbehavestate
end

if nargin<8,
    binsizelin = 2; % for calclinfields
end
if nargin<9,
    binsize2d = 1; % for openfieldrate
end

if nargin<10,
    figopt1 = 0; % Plot figures for each cell
end
if nargin<11,
    saveg1 = 0;
end

if nargin<12,
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
    %     linposfile = sprintf('%s/%slinpos%02d.mat', directoryname, prefix, day);
    %     load(linposfile);
    
    % Load spike file
    spkfile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
    load(spkfile);
    
    
    %% Get tetlist and celllist for current day
    %  Get tetrodes and cells automatically from epoch 1 of day   
    epoch=allepochs(1);
    tet = []; % New day, so restart tetlist and celllist
    
    ntottet = length(spikes{day}{epoch});
    for nt=1:ntottet
        if length(spikes{day}{epoch}{nt})~=0
            tet = [tet; nt];
        end
    end
    
    celllist=[]; % New day, so restart celllist
    celllist=cell(length(tet),1);
    
    
    
    %     for ep = 1:length(allepochs)
    %         epoch = allepochs(ep);
    
    %         statefile = sprintf('%s/%sstate%02d.mat', directoryname, prefix, day);
    %         if ~exist(statefile,'file');
    %             [state, lindist] = getbehavestate(linpos, day, epoch, includeStates, 'minlinvelocity', minvel);
    %             behavestate{day}{epoch}.state=state;
    %             behavestate{day}{epoch}.lindist=lindist;
    %             if ep==length(allepochs)
    %                 save(statefile,'behavestate');
    %             end
    %         else
    %             load(statefile);
    %             state = behavestate{day}{epoch}.state;
    %             lindist = behavestate{day}{epoch}.lindist;
    %         end
    %
    %         cmd = sprintf('state_day%02d_ep%02d=state;',day,epoch); eval(cmd);
    %         cmd = sprintf('lindist_day%02d_ep%02d=lindist;',day,epoch); eval(cmd);
    
    
    
    %tet = tetlist{d}; % tetrodes for current day
    %celllist = cellarray{d}; % cell list for tetrodes of current day
    
    for elec = 1:length(tet)
        
        % Get cells for tet
        currtet = tet(elec);
        epoch=allepochs(1);
        ntotcells = length(spikes{day}{epoch}{currtet});
        cells=[];
        for nc=1:ntotcells
            if length(spikes{day}{epoch}{currtet}{nc})~=0
                if isfield(spikes{day}{epoch}{currtet}{nc},'tag')
                    if (strcmp(spikes{day}{epoch}{currtet}{nc}.tag,'PyrSR') || strcmp(spikes{day}{epoch}{currtet}{nc}.tag,'PyrS') )
                        celllist{elec} = [celllist{elec}; nc];
                        cells = [cells;nc];
                    end
                else
                    celllist{elec} = [celllist{elec}; nc]; % If tag does not exist, add the cell
                    cells = [cells;nc];
                end
            end
        end
        
        %cells = celllist{elec};
        
        for neuron = 1:length(cells)
            currcell = cells(neuron);
            
            for ep = 1:length(allepochs)
                epoch=allepochs(ep);
                [placemap] = openfieldrate(spikes{day}{epoch}{currtet}{currcell}.data,pos{day}{epoch}.data,binsize2d,2);
                %cmd = sprintf('placemap_day%02d_ep%02d_tet%02d_cell%02d=placemap;',day,epoch,currtet,currcell); eval(cmd);
         
                if figopt1==1,
                    
                    figure(str2num([num2str(day+1) num2str(currtet) num2str(currcell)]));
                    redimscreen_2horsubplots;
                    hold on;                 
                    
                    %                     subplot(2,length(allepochs),ep); hold on;
                    %                     for i=1:length(trajdata),
                    %                         plot(trajdata{i}(:,5),[clr{i} '.-'],'Linewidth',2);
                    %                         %maxrate(i) = max(trajdata{i}(5:end-5,5));
                    %                     end
                    %                     if ep==1
                    %                         title(['Day ' num2str(day) '  Epoch ' num2str(epoch) '  Tet ' num2str(currtet) '  Cell ' num2str(currcell)],'FontSize',14,'Fontweight','bold');
                    %                         legend('InRight','OutRight','OutLeft','InLeft');
                    %                     else
                    %                         title(['Epoch ' num2str(epoch)],'FontSize',14,'Fontweight','bold');
                    %                     end
                               
                    subplot(1,length(allepochs),ep); hold on;
                    %figure; hold on;
                    imagesc(flipud(placemap.smoothedspikerate)); colorbar
                    set(gca, 'XTickLabel',[]);
                    set(gca, 'YTickLabel',[]);
                    %title(['Epoch ' num2str(epoch)], 'FontSize',14,'Fontweight','bold');
                    
                end
                
                if ep==1
                    cntneu=cntneu+1;
                    
                    allplacemaps_ep1{cntneu} = placemap;
                    % alllintrajs_ep1{cntneu} = trajdata;
                    % allpeakrates_ep1{cntneu} = maxrate;
                end
                
                if ep==2
                    cntneu_ep2=cntneu_ep2+1;
                    allplacemaps_ep2{cntneu_ep2} = placemap;
                    % alllintrajs_ep2{cntneu_ep2} = trajdata;
                    % allpeakrates_ep2{cntneu_ep2} = maxrate;
                end
                
                %trajplace{day}{ep}{currtet}{currcell}=trajdata;
                mapplace{day}{ep}{currtet}{currcell}=placemap;
                
                %clear trajdata placemap maxrate
                clear placemap
                
            end % end epoch
            
        end % end neuron
    end % end elec
    %   end % end epoch
    
    %dayplacefieldfile = sprintf('%s/%placefield%02d.mat', directoryname, prefix, day);
    %save(dayplacefieldfile,'trajplace','mapplace');
    
    %clear trajplace mapplace
    clear placemap
    
end % end day




% Find which trajectory had the max firing rate in all epochs
% for n=1:cntneu
%
%     peakrate=[]; peaktraj=[];
%
%     [peakrate(1),peaktraj1] = max(allpeakrates_ep1{n});
%     [peakrate(2),peaktraj2] = max(allpeakrates_ep2{n});
%
%
%     if peaktraj1(1)==peaktraj2(1),
%         allmaxtraj{n}=peaktraj1(1);
%     else
%         allmaxtraj{n}=[peaktraj1(1) peaktraj2(1)];
%     end
%
% end
%
% % Get corrleations between place fields
%
% single_traj=0; double_traj=0;
% for n=1:cntneu
%     n
%
%     % 2d corr
%     sizex = min([length(allplacemaps_ep1{n}.yticks),length(allplacemaps_ep2{n}.yticks)]);
%     sizey = min([length(allplacemaps_ep1{n}.xticks),length(allplacemaps_ep2{n}.xticks)]);
%     corr2d = corr2(allplacemaps_ep1{n}.smoothedspikerate(1:sizex,1:sizey), allplacemaps_ep2{n}.smoothedspikerate(1:sizex,1:sizey));
%
%
%     % Lin Traj Corr
%     maxtraj = allmaxtraj{n};
%     if length(maxtraj)==1
%         single_traj = single_traj+1;
%         place_ep1 = alllintrajs_ep1{n}{maxtraj}(:,5);
%         place_ep2 = alllintrajs_ep2{n}{maxtraj}(:,5);
%
%         invalid1 = find(isnan(place_ep1)==1);
%         invalid2 = find(isnan(place_ep2)==1);
%         invalid=unique([invalid1; invalid2]);
%         place_ep1(invalid)=[];
%         place_ep2(invalid)=[];
%
%         currcorr = corrcoef(place_ep1,place_ep2);
%         currcorr = currcorr(1,2);
%     else
%         double_traj = double_traj+1;
%         place1_ep1 = alllintrajs_ep1{n}{maxtraj(1)}(:,5);
%         place1_ep2 = alllintrajs_ep2{n}{maxtraj(1)}(:,5);
%
%         invalid1 = find(isnan(place1_ep1)==1);
%         invalid2 = find(isnan(place1_ep2)==1);
%         invalid=unique([invalid1; invalid2]);
%         place1_ep1(invalid)=[];
%         place1_ep2(invalid)=[];
%         currcorr1 = corrcoef(place1_ep1,place1_ep2);
%         currcorr1 = currcorr1(1,2);
%
%         place2_ep1 = alllintrajs_ep1{n}{maxtraj(2)}(:,5);
%         place2_ep2 = alllintrajs_ep2{n}{maxtraj(2)}(:,5);
%         invalid1 = find(isnan(place2_ep1)==1);
%         invalid2 = find(isnan(place2_ep2)==1);
%         invalid=unique([invalid1; invalid2]);
%         place2_ep1(invalid)=[];
%         place2_ep2(invalid)=[];
%         currcorr2 = corrcoef(place2_ep1,place2_ep2);
%         currcorr2 = currcorr2(1,2);
%
%         currcorr = mean([currcorr1,currcorr2]);
%
% %         place_ep1=sum(place1_ep1,place2_ep1);
% %         place_ep2=sum(place1_ep2,place2_ep2);
% %         currcorr = corrcoef(place_ep1,place_ep2);
% %         currcorr = currcorr(1,2)
%
%     end
%
%     allcorr(n)=currcorr;
%     all2dcorr(n)=corr2d;
%
% end
%
% disp('Mean LinTrajCorr');
% nanmean(allcorr),
% nanstd(allcorr),
%
% disp('Mean 2dCorr');
% nanmean(all2dcorr),
% nanstd(all2dcorr),












