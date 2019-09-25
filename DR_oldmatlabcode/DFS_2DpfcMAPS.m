
warning('off','all');
clear; close all;

runscript = 0;
savedata = runscript; % save data option - only works if runscript is also on
savefigs=0;
cyclemaps =1;
fonttype = 'Arial';
titlesize = 16;
axissize = 16;
plot2dPFC = 1;
savedir = '/mnt/data25/sjadhav/HPExpt/ProcessedDataDR/';
savefilename = '2DPFCmaps_040214_3occ_3spk_thresh002'; %has to match saved data name. specify GLM or Corr pairs and filters used: velocity filter <=.....nrip >= (#tetrodes ripples detected); std > of ripple detection power
savefigfilename = '2DPFCmaps_050214_3occ_3spike_thresh002'; %fig name
savefile = [savedir savefilename];
figdir = '/mnt/data25/sjadhav/HPExpt/Figures_DR/';
mkdir([figdir savefigfilename]);

% If runscript, run Datafilter and save data
if runscript == 1
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPa' 'HPb' 'HPc' 'nadal'};
    
    %Filter creation
    %-----------------------------------------------------
    % Epoch filter
    % -------------
    dayfilter = ''; %leave blank to take all days from HP animals and Ndl
    
    runepochfilter = 'isequal($type, ''run'') && ~isequal($environment, ''lin'')';
    
    % %Cell filter
    % %-----------
    placecellfilter = '( strcmp($area, ''PFC'')) && ($numspikes > 100)';
    
    % Time filter -
    %%-----------
    riptetfilter = '(isequal($descrip, ''riptet''))';
    timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 3))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    spatf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cells',placecellfilter,'excludetime', timefilter,'iterator',iterator);
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    twoDfields = setfilterfunction(spatf, 'DFAsj_twodoccupancy_DR', {'spikes', 'linpos', 'pos'});
    
    % Run analysis-----------------------
    twoDflds = runfilter(twoDfields);  % 2D Place Field Map by Trajs
    
    %     end
    disp('Finished running filter script');
    if savedata == 1
        clear savefigfilename runscript  savedata cyclemaps savefigs figdir savefilename savedir fonttype titlesize axissize plot2dPFC
        save(savefile);
    end
else
    load(savefile);
end % end runscript

if ~exist('savedata')
    return
end
%--------------------- Finished Filter Function Run -------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% make a subplot with the background grey in the upper 4 panels and the trajectory in the lower 4 panels... in ai just overlay
cmap = colormap;
% cmap((ceil(0.8*length(cmap(:,1))):end), :) = repmat(cmap(end,[1 2 3]),length(cmap(:,1))-(floor(0.8*length(cmap(:,1)))),1);  %saturate the highest 20% rate
cmap(1,:) = [.8 .8 .8]; %make backtrack gray
cmap(end,:) = [1 1 1]; %make background white
cellcount = 0;
for i = 1:length(twoDflds); %per anim
    for j = 1:length(twoDflds(i).output{1}); %per cell
        
        cellcount = cellcount +1;
        %create grey background from all trajs
        backgrey = [];
        skipp = 0;
        %create new xticks and yticks to accomidate the 2 cm bins..
                
        
        
        for ii = 1:4;
%             currYticks = [twoDflds(i).output{1}(j).yticks{ii}(1):1:twoDflds(i).output{1}(j).yticks{ii}(1)+length(twoDflds(i).output{1}(j).yticks{ii})-1];
%             currXticks = [twoDflds(i).output{1}(j).xticks{ii}(1):1:twoDflds(i).output{1}(j).xticks{ii}+length(twoDflds(i).output{1}(j).xticks{ii})-1];
            
          
            try% ~isempty(twoDflds(i).output{1}(j).smoothedoccupancy{ii});
%                 backgrey(twoDflds(i).output{1}(j).yticks{ii}, twoDflds(i).output{1}(j).xticks{ii}) =  twoDflds(i).output{1}(j).smoothedoccupancy{ii};
                backgrey(twoDflds(i).output{1}(j).xticks{ii},twoDflds(i).output{1}(j).yticks{ii}) =  twoDflds(i).output{1}(j).smoothedoccupancy{ii};
%                 backgrey(twoDflds(i).output{1}(j).yticks{ii}, twoDflds(i).output{1}(j).xticks{ii}) =  twoDflds(i).output{1}(j).occupancy{ii};
%                 backgrey(currYticks, currXticks) =  twoDflds(i).output{1}(j).smoothedoccupancy{ii};
            catch
                skipp = 1;
            end
        end
        
        if skipp == 1;
            twoDflds(i).output{1}(j).index{1}(:)
            break
        end
        
        backgrey(:,:) = 0;  %the above is a hack to get the dimensions of what backgrey will be so i don't overwrite previous traj vals in the below
        for ii = 1:4;
%                         currYticks = [twoDflds(i).output{1}(j).yticks{ii}(1):1:twoDflds(i).output{1}(j).yticks{ii}(1)+length(twoDflds(i).output{1}(j).yticks{ii})-1];
%             currXticks = [twoDflds(i).output{1}(j).xticks{ii}(1):1:twoDflds(i).output{1}(j).xticks{ii}+length(twoDflds(i).output{1}(j).xticks{ii})-1];
            
            twoDflds(i).output{1}(j).smoothedspikerate{ii}(twoDflds(i).output{1}(j).smoothedspikerate{ii}==-1) = 0;
            backgrey(twoDflds(i).output{1}(j).xticks{ii}, twoDflds(i).output{1}(j).yticks{ii}) = backgrey(twoDflds(i).output{1}(j).xticks{ii}, twoDflds(i).output{1}(j).yticks{ii}) + twoDflds(i).output{1}(j).smoothedoccupancy{ii};
%             backgrey(twoDflds(i).output{1}(j).yticks{ii}, twoDflds(i).output{1}(j).xticks{ii}) = backgrey(twoDflds(i).output{1}(j).yticks{ii}, twoDflds(i).output{1}(j).xticks{ii}) + twoDflds(i).output{1}(j).occupancy{ii};
            backgrey(backgrey==-1) = 0; %need to do this so the negatives don't summate
            
            %              get a maxrate to normalize in the next loop
            currtraj = twoDflds(i).output{1}(j).smoothedspikerate{ii};
            maxeachtraj(ii) = max(max(currtraj));
        end
        
        
        
        maxrate = max(max(maxeachtraj));
        backgrey(backgrey == 0) = -1; %some rates are more negative than -1
        backgrey(backgrey > 0) = 0;
        backgrey(backgrey ==-1) = 1.0001; % will be max
        
        
        figure; hold on;
        for ii = 1:4; %trajs
%             currYticks = [twoDflds(i).output{1}(j).yticks{ii}(1):1:twoDflds(i).output{1}(j).yticks{ii}(1)+length(twoDflds(i).output{1}(j).yticks{ii})-1];
%             currXticks = [twoDflds(i).output{1}(j).xticks{ii}(1):1:twoDflds(i).output{1}(j).xticks{ii}+length(twoDflds(i).output{1}(j).xticks{ii})-1];
            currplot = subplot(1,5,ii);
            fronttraj = zeros(length(backgrey(:,1)),length(backgrey(1,:)));
            currtraj = twoDflds(i).output{1}(j).smoothedspikerate{ii};
%             maxeachtraj(ii) = max(max(currtraj));
            currtraj(currtraj == -1) = 0;
            currtraj = (currtraj)./maxrate; %% this is wrong... i need to be dividing everything by the max rate of all trajs.. not just the current one..
            fronttraj(twoDflds(i).output{1}(j).xticks{ii}, twoDflds(i).output{1}(j).yticks{ii}) = currtraj;
            combtraj = fronttraj+backgrey;
% %             
%             fixedtraj = combtraj(2:2:end, 1:2:end);
%             imagesc(fixedtraj);

if i == 1 || i ==4; %HPA or nadal
            imagesc(combtraj');
                        widthvec = min(backgrey,[],2);
            heightvec = min(backgrey,[],1)';
            minX = find(widthvec==0,1, 'first');
            maxX = find(widthvec==0,1,'last');
            minY = find(heightvec==0,1, 'first');
            maxY = find(heightvec==0,1,'last');
            axis([minX maxX minY maxY]);
else
    imagesc(combtraj)
                widthvec = min(backgrey,[],1);
            heightvec = min(backgrey,[],2)';
            minX = find(widthvec==0,1, 'first');
            maxX = find(widthvec==0,1,'last');
            minY = find(heightvec==0,1, 'first');
            maxY = find(heightvec==0,1,'last');
            axis([minX maxX minY maxY]);
end
            colormap(cmap);
            
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            
            lastpos = get(currplot,'position');
            
        end
        
%         maxalltraj = max(maxeachtraj);
        subplot(1,5,5);
        blankp = backgrey(:,1);
        blankp(:,:) = 1.0001;
        plot(blankp);
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        
        set(gcf,'PaperPositionMode','auto');
        set(gcf, 'Position', [100 100 950 300])
        tit = sprintf('%s Peak(%0.2f) %d %d %d %d', twoDflds(i).animal{1}, maxrate,  twoDflds(i).output{1}(j).index{1}(1), twoDflds(i).output{1}(j).index{1}(2), twoDflds(i).output{1}(j).index{1}(3), twoDflds(i).output{1}(j).index{1}(4));
        %         tit = sprintf('%s __%d Peak(%0.5f)', twoDflds(i).animal{1}, twoDflds(i).output{1}(j).index{1}(:), maxrate);
        %         supertitle(tit);
        %         supertitle(tit,  'Fontsize', titlesize,'FontName',fonttype);
        %         set(currplot, 'position', lastpos)
        supertitle(tit);
        colorbar('location', 'EastOutside')
        if cyclemaps == 1;
            keyboard
        end
        figfile = [figdir,savefigfilename,'/',sprintf('%s_%d', twoDflds(i).animal{1}, cellcount)];
        if savefigs==1
            print('-dpng', figfile, '-r300');
            print('-dpdf', figfile, '-r300');
        end
        close all
    end
end
















