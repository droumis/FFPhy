function [allcorr,all2dcorr] = sj_placefieldstab_epochs_1 (animdirect,prefix,days,allepochs,includeStates, minvelstate, minabsvel, binsizelin, binsize2d, figopt1,saveg1,savedata)
%% Generate and Plot Place fields for cells on tetrode

% Shantanu: Change to get tet list and cell list automatically from day file
% 06/02/2011

% sj_placefieldstab_epochs_1('/data25/sjadhav/RippleInterruption/REf_direct','REf',[1:8],[2 4], 6, 5, 5, 2, 1, 0,0,0);
% sj_placefieldstab_epochs_1('/data25/sjadhav/RippleInterruption/REf_direct','REf',[8],[2 4], 6, 3, 3, 2, 1, 1,0,0);

% With tetlist and cellarray fed in
% sj_placefieldstab_epochs_1('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',[3 7 8],[2 4],{[1 5 7],1,7}, {{[2:3],[2 3 5],[2:4]},{[2:4]},{1:2}}, 3, 6, 2, 1, 1,0,0);
% sj_placefieldstab_epochs_1('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',3,[2 4],{[1 5 7]}, {{[2:3],[2 3 5],[2:5]}}, 3, 6, 2, 1, 1,0,0);
% sj_placefieldstab_epochs_1('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',3,[1 3 5],'Sleep',5, {2:5}, 3,6,2,1, 1,0,0);
% sj_placefieldstab_epochs_1('/data25/sjadhav/RippleInterruption/REf_direct','REf',9,[2 4],'Run',{[9]}, {{1}}, 3, 6, 2, 1, 1,0,0);

if nargin<5,
    includeStates = 6; % for sj_getbehavestate
end
if nargin<6,
    minvelstate = 3; % for sj_getbehavestate
end
if nargin<7,
    minabsvel = 3; % for sj_getbehavestate
end
if nargin<8,
    binsizelin = 2; % for sj_calclinfields
end
if nargin<9,
    binsize2d = 1; % for sj_openfieldrate
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


novel=1:3;
fam=4:8;

cntneu=0;  % Count neurons across days
cntneu_ep2=0;
cnt_allneu=0;
allcorr=[]; all2dcorr=[];
allcorr_day=[]; all2dcorr_day=[];
for d = 1:length(days)
    
    day = days(d)
    cnt_dayneu = 0;
    
    % Load pos file
    posfile = sprintf('%s/%spos%02d.mat', directoryname, prefix, day);
    load(posfile);
    
    % Load linpos file
    linposfile = sprintf('%s/%slinpos%02d.mat', directoryname, prefix, day);
    load(linposfile);
    
    % Load spike file
    spkfile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
    load(spkfile);
    
    % You will load or calculate state file later
    
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
            cnt_dayneu = cnt_dayneu + 1;
            
            for ep = 1:length(allepochs)
                
                epoch = allepochs(ep);
                state = linpos{day}{epoch}.statematrix.traj;
                lindist = linpos{day}{epoch}.statematrix.lindist;
                
                cmd = sprintf('state_day%02d_ep%02d=state;',day,epoch); eval(cmd);
                cmd = sprintf('lindist_day%02d_ep%02d=lindist;',day,epoch); eval(cmd);
                
                % Abs velocity filter shoul be applied in getbehavestate, and Check here again
                trajdata = sj_calclinfields(spikes,state,lindist,linpos,pos,[day epoch currtet currcell],binsizelin,minabsvel);
                %cmd = sprintf('trajdata_day%02d_ep%02d_tet%02d_cell%02d=trajdata;',day,epoch,currtet,currcell); eval(cmd);
                
                std=2;  % std for Gaussian smoothing of spike rate
                [placemap] = sj_openfieldrate(spikes,state,pos,[day epoch currtet currcell],binsize2d,std,minabsvel);
                %[placemap] = openfieldrate(spikes{day}{epoch}{currtet}{currcell}.data,pos{day}{epoch}.data,binsize2d,2);
                %cmd = sprintf('placemap_day%02d_ep%02d_tet%02d_cell%02d=placemap;',day,epoch,currtet,currcell); eval(cmd);
                
                % Get traj with maximum peak fir rate
                for i=1:length(trajdata),
                    maxrate(i) = max(trajdata{i}(5:end-5,5));
                end
                
                if figopt1==1,
                    
                    figure(str2num([num2str(day+1) num2str(currtet) num2str(currcell)]));
                    redimscreen_2x2subplots;
                    hold on;
                    
                    subplot(2,length(allepochs),ep); hold on;
                    for i=1:length(trajdata),
                        plot(trajdata{i}(:,5),[clr{i} '.-'],'Linewidth',2);
                        %maxrate(i) = max(trajdata{i}(5:end-5,5));
                    end
                    if ep==1
                        title(['Day ' num2str(day) '  Epoch ' num2str(epoch) '  Tet ' num2str(currtet) '  Cell ' num2str(currcell)],'FontSize',14,'Fontweight','bold');
                        legend('InRealRight','OutRealLeft','OutRealRight','InRealLeft');
                        xlabel ('Position along linear trajectory (cm)','FontSize',14,'Fontweight','bold');
                        ylabel ('Firing Rate (Hz)','FontSize',14,'Fontweight','bold');
                        
                    else
                        title(['Epoch ' num2str(epoch)],'FontSize',14,'Fontweight','bold');
                        xlabel ('Position along linear trajectory (cm)','FontSize',14,'Fontweight','bold');
                    end
                    
                    subplot(2,length(allepochs),length(allepochs)+ep); hold on;
                    % Flipping image ud OR lr AFTER axis is set by hold on gives correct direction
                    imagesc(flipud(placemap.smoothedspikerate)); colorbar
                    title(['Epoch ' num2str(epoch)], 'FontSize',14,'Fontweight','bold');
                    if ep==1
                        xlabel ('X-position (cm)','FontSize',14,'Fontweight','bold');
                        ylabel ('Y-position (cm)','FontSize',14,'Fontweight','bold');
                    else
                        xlabel ('X-position (cm)','FontSize',14,'Fontweight','bold');
                    end
                    
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
                
            end % end epoch
            
        end % end neuron
    end % end elec
    
    %dayplacefieldfile = sprintf('%s/%placefield%02d.mat', directoryname, prefix, day);
    %save(dayplacefieldfile,'trajplace','mapplace');
    %clear trajplace mapplace
    
    % Get correlation data for current day neurons
    for n=1:cnt_dayneu
        
        cnt_allneu = cnt_allneu + 1;  % Neuron Count across days - Not reset
        
        % Find which trajectory had the max firing rate in all epochs
        peakrate=[]; peaktraj=[];      
        [peakrate,peaktraj1] = max(allpeakrates_ep1{n});
        [peakrate,peaktraj2] = max(allpeakrates_ep2{n});
        
        if peaktraj1(1)==peaktraj2(1),
            allmaxtraj{n}=peaktraj1(1);
        else
            allmaxtraj{n}=[peaktraj1(1) peaktraj2(1)];
        end     
        
        % Get corrleations between place fields
        single_traj=0; double_traj=0;
        
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
            le = min(length(place_ep1),length(place_ep2));
            place_ep1 = place_ep1(1:le);
            place_ep2 = place_ep2(1:le);
            
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
            le = min(length(place1_ep1),length(place1_ep2));
            place1_ep1 = place1_ep1(1:le);
            place1_ep2 = place1_ep2(1:le);
            
            invalid1 = find(isnan(place1_ep1)==1);
            invalid2 = find(isnan(place1_ep2)==1);
            invalid=unique([invalid1; invalid2]);
            place1_ep1(invalid)=[];
            place1_ep2(invalid)=[];
          
            currcorr1 = corrcoef(place1_ep1,place1_ep2);
            currcorr1 = currcorr1(1,2);
            
                        
            place2_ep1 = alllintrajs_ep1{n}{maxtraj(2)}(:,5);
            place2_ep2 = alllintrajs_ep2{n}{maxtraj(2)}(:,5);
            le = min(length(place2_ep1),length(place2_ep2));
            place2_ep1 = place2_ep1(1:le);
            place2_ep2 = place2_ep2(1:le);
            
            
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
    
        allcorr(cnt_allneu)=currcorr;
        all2dcorr(cnt_allneu)=corr2d;
        
        allcorr_day{day}(n)=currcorr;
        all2dcorr_day{day}(n)=corr2d;
        
    end % end cnt_dayneu
   
    
end % end day

disp('Mean LinTrajCorr');
meancorr = nanmean(allcorr),
errcorr = nansem(allcorr),

% disp('Mean 2dCorr');
% nanmean(all2dcorr),
% nansem(all2dcorr),


%% Summarize correlation - 1) Bar plot for all, novel, familiar days, & 2) Corrln across days

set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','normal');
set(0,'defaultaxeslinewidth',2);

% 1) Bar plot for all, novel, familiar days
figure; hold on; redimscreen_figforppt1;
bar(1,meancorr,'k'); errorbar(1,meancorr,errcorr,'k');
set(gca,'XTick',[1],'XTickLabel',{'All Days'},'FontSize',16,'Fontweight','normal');
axis([0 2 0 1.0])

if ~isempty(intersect(days,novel))
    noveldays = intersect(days,novel);
    novelcorr = [];
    for nd=noveldays
        novelcorr = [novelcorr, allcorr_day{nd}];
    end
    
    bar(2,nanmean(novelcorr),'r'); errorbar(2,nanmean(novelcorr),nansem(novelcorr),'r');
    set(gca,'XTick',[1 2],'XTickLabel',{'All Days','Novel Days'},'FontSize',16,'Fontweight','normal');
    axis([0 3 0 1.0])   
end

if ~isempty(intersect(days,fam))
    famdays = intersect(days,fam);
    famcorr = [];
    for fd=famdays
        famcorr = [famcorr, allcorr_day{fd}];
    end
    
    if isempty(intersect(days,novel))
        bar(2,nanmean(famcorr),'g'); errorbar(2,nanmean(famcorr),nansem(famcorr),'g');
        set(gca,'XTick',[1 2],'XTickLabel',{'All Days','Fam Days'},'FontSize',16,'Fontweight','normal');
        axis([0 3 0 1.0])       
    else
        
        bar(3,nanmean(famcorr),'g'); errorbar(3,nanmean(famcorr),nansem(famcorr),'g');
        set(gca,'XTick',[1 2 3],'XTickLabel',{'All Days','Novel Days','Fam Days'},'FontSize',16,'Fontweight','normal');
        axis([0 4 0 1.0])
    end
end

title(['Place Field Stability - Days ' num2str(days) ],'FontSize',20,'Fontweight','normal');
ylabel('Place Field Corrln across epochs','FontSize',20,'Fontweight','normal');
    
% 2) Corrln across days

if length(days)>=4
    
    cnt=0;
    for nd = 1:length(days)-2
        
        d = days(nd);
        cnt=cnt+1;
        dayscorr = [];
        currdays = [d, d+1, d+2];
        dayscorr = [dayscorr, allcorr_day{d}, allcorr_day{d+1}, allcorr_day{d+2}]; 
        cumdaycorr{nd} = dayscorr;
        cumdaycorrmean(nd) = nanmean(dayscorr);
        cumdaycorrerr(nd) = nansem(dayscorr);
        startday(nd) = d;
    end
    
    figure; hold on; redimscreen_figforppt1;
    errorbar(startday, cumdaycorrmean, cumdaycorrerr,'ko-','Linewidth',2,'MarkerSize',8);
    title('Place Field Stability Across Days','FontSize',24,'Fontweight','normal');
    ylabel('Place Field Corrln across epochs','FontSize',20,'Fontweight','normal');
    xlabel('Start Day','FontSize',16,'Fontweight','normal');
    axis([startday(1)-0.5 startday(end)+0.5 0 1.0])
    
end
i=1;


%keyboard;




