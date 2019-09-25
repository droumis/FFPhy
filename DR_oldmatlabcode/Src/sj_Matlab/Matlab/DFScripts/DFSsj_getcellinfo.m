
% Just using the cellinfo structure to get and plot cell properties


clear;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 1; % Figure Options

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
%savefile = [savedir 'HP_cellinfo_CA1']; %
savefile = [savedir 'HP_cellinfo_PFC']; %


% If runscript, run Datafilter and save data
if runscript == 1
    
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPa','HPb','HPc'};
    
    %Filter creation
    %--------------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    %epochfilter = 'isequal($type, ''run'')'; % All run environments
    epochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'')'; % Start with just W-tracks run
    
    % Cell filter
    % -----------
    
    % Any cell defined in at least one run in environment
    
    cells = '(strcmp($area, ''PFC'')) &&  ($numspikes > 0)';
    %cells = '( strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) &&  ($numspikes > 0)';

    
    % Time filter
    % -----------    
    % None
    
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    cellsf = createfilter('animal',animals,'days',dayfilter,'epochs',epochfilter,'cells',cells,'iterator', iterator);
    
    % Set analysis function
    % ----------------------
    cellsf = setfilterfunction(cellsf, 'DFAsj_getcellinfo', {'cellinfo'}); % You can also use "cellprop" which uses spikes as the variable and calls getcellprop

    
    % Run analysis
    % ------------
    cellsf = runfilter(cellsf);  % 
   
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript savedata 
        save(savefile);
    end
    
else
    
    load(savefile);
    
end  % end runscript

if ~exist('savedata')
    return
end


% --------------------------
% Post-filter run analysis
% --------------------------

% ---------------------------------------------------------
% Gather data and append animal index to it as well
% Unlike codes where you want to keep data for epochs separate for plotting (eg. place field code), here, combine everything for population analysis.
% ---------------------------------------------------------

alldata = []; alltags = []; allanimindex = []; % includes animal index 
cnt = 0; % For cnting across animals for tags
for an = 1:length(animals)
    for i=1:length(cellsf(an).output{1}),
        anim_index{an}(i,:)=cellsf(an).output{1}(i).index;
        % Only indexes
        animindex=[an cellsf(an).output{1}(i).index]; % Put animal index in front
        allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
        
        % Indexes + Data [anim day epoch tet cell rate csi propbursts numspikes spikewidth
        append = [an cellsf(an).output{1}(i).index cellsf(an).output{1}(i).meanrate cellsf(an).output{1}(i).csi ...
            cellsf(an).output{1}(i).propbursts cellsf(an).output{1}(i).numspikes cellsf(an).output{1}(i).spikewidth];
        alldata = [alldata; append];
         
        % Make a tag field with the same length
        cnt=cnt+1;
        allrawtags(cnt).tag = cellsf(an).output{1}(i).tag;
        
    end
end

% -----------------------------------------------------------------
% Consolidate single cells' across epochs in a day .celloutput field
% Can do for each animal and then store separately by using the index from output field - Do in loop for each animal, 
% or can use the above appended animal+index field to parse data 
% -----------------------------------------------------------------

% Save the consolidated data back in a structure - TO put back in filter
% structure, has to be animal-wise. Sell alternate method below
celloutput = struct;
dummyindex=allanimindex;  % all anim-day-epoch indices


cntcells=0;
for i=1:size(alldata,1)
    animdaytetcell=alldata(i,[1 2 4 5]);
    ind=[];
    while rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))~=0          % collect all rows (epochs)
        ind = [ind rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))];        % finds the first matching row
        dummyindex(rowfind(animdaytetcell,dummyindex(:,[1 2 4 5])),:)=[0 0 0 0 0]; % after adding index, remove the corresponding row
                                                                           % so you could find the next one 
    end
    
    % Gather everything for the current cell across epochs
    currrate=[]; currcsi=[]; currpropbursts=[]; currnumspikes=[]; currspikewidth=[]; currtag = [];
    for r=ind
        currrate=[currrate ; alldata(r,6)];
        currcsi=[currcsi; alldata(r,7)];
        currpropbursts=[currpropbursts; alldata(r,8)];
        currnumspikes=[currnumspikes;alldata(r,9)];
        currspikewidth=[currspikewidth;alldata(r,10)];
    end    
    
    if ~isempty(currrate)
        currtag = allrawtags(ind(1)).tag; % Tag is same for all epochs
        cntcells = cntcells + 1;
        celloutput(cntcells).index=animdaytetcell;
        celloutput(cntcells).meanrate=currrate; % Remember, this is seperated across epochs. Can take mean
        celloutput(cntcells).csi=currcsi;
        celloutput(cntcells).propbursts=currpropbursts;
        celloutput(cntcells).numspikes=currnumspikes;
        celloutput(cntcells).spikewidth=currspikewidth;
        celloutput(cntcells).tag=currtag;
    end   
end

% Alternate method - Keep Animal data separate. For each animal, look within dayeptetcell
% ----------------------------------------------------------------------------------------
% for an = 1:length(animals)
%     dummyindex = anim_index{an};     % collect all day-epoch-tet-cell indices    
%     for i=1:length(cellsf(an).output{1})
%         daytetcell=cellsf.output{1}(i).index([1 3 4]);
%         ind=[];
%         while rowfind(daytetcell,dummyindex(:,[1 3 4]))~=0          % collect all rows (epochs)
%             ind = [ind rowfind(daytetcell,dummyindex(:,[1 3 4]))];
%             dummyindex(rowfind(daytetcell,dummyindex(:,[1 3 4])),:)=[0 0 0 0];
%         end
%         % Gather everything for the current cell across epochs
%         currrate=[]; currcsi=[]; currpropbursts=[]; currnumspikes=[]; currspikewidth=[]; currtag = [];
%         for r=ind
%             currrate=[currrate ; cellsf.output{1}(r).meanrate];
%             currcsi=[currcsi; cellsf.output{1}(r).csi];
%             currpropbursts=[currpropbursts; cellsf.output{1}(r).propbursts];
%             currnumspikes=[currnumspikes;cellsf.output{1}(r).numspikes];
%             currspikewidth=[currspikewidth;cellsf.output{1}(r).spikewidth];
%         end
%         if ~isempty(currrate)
%             currtag = cellsf.output{1}(ind(1)).tag; % Tag is same for all epochs
%             cellsf(an).celloutput.index=daytetcell;
%             cellsf(an).celloutput.meanrate=currrate;
%             cellsf(an).celloutput.csi=csi;
%             cellsf(an).celloutput.propbursts=propbursts;
%             cellsf(an).celloutput.numspikes=numspikes;
%             cellsf(an).celloutput.spikewidth=spikewidth;
%             cellsf(an).celloutput.tag=ncurrtag;
%         end 
%     end
% end


% Values for cells - take means across epochs
allmeanrate = []; allcsi = []; allpropbursts = []; allnumspikes = []; allspikewidth = []; Intidx = [];
for i=1:length(celloutput)
    currtag = celloutput(i).tag;
    allcellidx(i,:) = celloutput(i).index;
    allmeanrate = [allmeanrate; mean(celloutput(i).meanrate)];
    allcsi = [allcsi; mean(celloutput(i).csi)];
    allpropbursts = [allpropbursts; mean(celloutput(i).propbursts)];
    allspikewidth = [allspikewidth, mean(celloutput(i).spikewidth)];
    if strcmp(currtag,'CA1Int') || strcmp(currtag,'iCA1Int')
        Intidx = [Intidx;i];
    end
end

% Gather tag indices - For Hippocmapus
% -------------------------------------
% Pyridx = []; Intidx = []; Partidx = []; CA1idx = []; iCA1idx = [];
% for i=1:length(celloutput)
%         
%     currtag = celloutput(i).tag;
%     if strcmp(currtag,'CA1Pyr') || strcmp(currtag,'iCA1Pyr') || strcmp(currtag,'CA1Run') || strcmp(currtag,'iCA1Run')
%         Pyridx = [Pyridx;i];
%     end
%     if strcmp(currtag,'CA1Int') || strcmp(currtag,'iCA1Int')
%         Intidx = [Intidx;i];
%     end
%     if strcmp(currtag,'CA1Pyrpr') || strcmp(currtag,'CA1Pyrps') || strcmp(currtag,'CA1Pyrp') || strcmp(currtag,'iCA1Pyrpr') || strcmp(currtag,'iCA1Pyrps') || strcmp(currtag,'iCA1Pyrp')
%         Partidx = [Partidx;i];
%     end
%     if strcmp(currtag,'CA1Pyr') || strcmp(currtag,'CA1Run') || strcmp(currtag,'CA1Pyrpr') || strcmp(currtag,'CA1Pyrps') || strcmp(currtag,'CA1Pyrp')
%         CA1idx = [CA1idx;i];
%     end
%     if strcmp(currtag,'iCA1Pyr') || strcmp(currtag,'iCA1Run') || strcmp(currtag,'iCA1Pyrpr') || strcmp(currtag,'iCA1Pyrps') || strcmp(currtag,'iCA1Pyrp')
%         iCA1idx = [iCA1idx;i];
%     end
%               
% end



% ---------
% Plotting
% ---------

ncells = length(allmeanrate)
ncells_1Hz = length(find(allmeanrate<1))
ncells_point1Hz = length(find(allmeanrate<0.1))
ncells_point01Hz = length(find(allmeanrate<0.01))

if ~exist('figopt1'), figopt1=0; end
if figopt1==1,
    
    %1) All Cells - Fir rates vs spikewidths
    figure; hold on;
    redimscreen_figforppt1;
    
    %scatter(allspikewidth, allmeanrate, [],'r');
    plot(allspikewidth, allmeanrate,'ro','MarkerSize',8,'LineWidth',1.5);
    %plot(allspikewidth(Intidx), allmeanrate(Intidx),'kx','MarkerSize',8,'LineWidth',1.5);

    xaxis = min(allspikewidth):0.1:max(allspikewidth);
    plot(xaxis,0.01*ones(size(xaxis)),'c-','Linewidth',1);
    plot(xaxis,0.1*ones(size(xaxis)),'k-','Linewidth',1);
    plot(xaxis,7*ones(size(xaxis)),'k-','Linewidth',1);   
    
    set(gca,'YLim',[0 1.05*max(allmeanrate)]);
    set(gca,'XLim',[0.95*min(xaxis) 1.05*max(xaxis)]);
    
    xlabel('Spike Widths','FontSize',24,'Fontweight','normal');
    ylabel('Firing Rates (Hz)','FontSize',24,'Fontweight','normal');
    title('Cell Properties','FontSize',24,'Fontweight','normal')

    
   %2) Numspikes 
   figure; hold on;
    redimscreen_figforppt1;
    for i=1:length(celloutput)
        currnspk = celloutput(i).numspikes;
        currtag = celloutput(i).tag;
        xaxis = i*ones(size(currnspk));
        %plot(xaxis,currnspk,'ro','MarkerSize',8,'LineWidth',1.5);
        %line([i i], [0 max(currnspk)]); % Line for current cell through all its values
        plot(i, mean(currnspk), 'gsq', 'MarkerSize',10,'LineWidth',2); % Mean value marked

        if strcmp(currtag,'CA1Int') || strcmp(currtag,'iCA1Int');
            plot(xaxis,currnspk,'ko','MarkerSize',8,'LineWidth',1.5);
            %line([i i], [0 max(currnspk)]); % Line for current cell through all its values
            plot(i, mean(currnspk), 'ksq', 'MarkerSize',8,'LineWidth',2); % Mean value marked
        end
  
    end
    xline = 1:length(celloutput);
    plot(xline,100*ones(size(xline)),'m-','Linewidth',1)
    
    
    %3) Numspikes in run epochs for cells, or meanrate separately for epochs
    figure; hold on;
    redimscreen_figforppt1;
    for i=1:length(celloutput)
        currnspk = celloutput(i).numspikes;
        currtag = celloutput(i).tag;
        xaxis = i*ones(size(currnspk));
        plot(xaxis,currnspk,'ro','MarkerSize',8,'LineWidth',1.5);
        %line([i i], [0 max(currnspk)]); % Line for current cell through all its values
        plot(i, mean(currnspk), 'gsq', 'MarkerSize',10,'LineWidth',2); % Mean value marked
        
        if strcmp(currtag,'CA1Int') || strcmp(currtag,'iCA1Int');
            plot(xaxis,currnspk,'ko','MarkerSize',8,'LineWidth',1.5);
            %line([i i], [0 max(currnspk)]); % Line for current cell through all its values
            plot(i, mean(currnspk), 'ksq', 'MarkerSize',8,'LineWidth',2); % Mean value marked
        end
  
    end
    % Make line at 100 spikes (0.1Hz rate for 20 min run session is ~120 spikes )
    xaxis = 1:length(celloutput);
    plot(xaxis,100*ones(size(xaxis)),'k-','Linewidth',1);
    
    
    %allcellidx(find(allmeanrate<0.1),:)
    %allmeanrate(find(allmeanrate<0.1),:)
    %celloutput(find(allmeanrate<0.1)).numspikes
    %allcellidx(Intidx,:)
    %allmeanrate(Intidx)
    
    
    
%       %3) All Cells - 3d. Fir rate vs Spike WIdth vs Burst Proportion
%     figure; hold on;
%     redimscreen_figforppt1;
%     
%     plot3(allspikewidth, allmeanrate, allpropbursts,'ro');
%     %plot3(allspikewidth(Intidx), allmeanrate(Intidx), allpropbursts(Intidx),'kx');
% 
%     %set(gca,'YLim',[0 1.05*max(allpropbursts)]);
%     %set(gca,'XLim',[0.95*min(xaxis) 1.05*max(xaxis)]);
%     
%     xlabel('Spike Widths','FontSize',24,'Fontweight','bold');
%     ylabel('Firing Rates (Hz)','FontSize',24,'Fontweight','bold');
%     zlabel('Burst Proportion','FontSize',24,'Fontweight','bold');
%     title('Cell Properties','FontSize',24,'Fontweight','bold');
    
    
%       %2) All Cells - CSI vs Burst Proportion
%     figure; hold on;
%     redimscreen_figforppt1;
%     
%     plot(allcsi, allpropbursts,'ro','MarkerSize',8,'LineWidth',1.5);
%     %plot(allcsi(Intidx), allpropbursts(Intidx),'kx','MarkerSize',8,'LineWidth',1.5);
% 
%     xaxis = min(allcsi):0.1:max(allcsi);    
%     set(gca,'YLim',[0 1.05*max(allpropbursts)]);
%     set(gca,'XLim',[0.95*min(xaxis) 1.05*max(xaxis)]);
%     
%     xlabel('CSI','FontSize',24,'Fontweight','bold');
%     ylabel('Burst Proportion','FontSize',24,'Fontweight','bold');
%     title('Cell Properties','FontSize',24,'Fontweight','bold')
    
    
%     
%     %4) CA1 vs iCA1 
%     figure; hold on;
%     redimscreen_figforppt1;
%     plot(allspikewidth(CA1idx), allmeanrate(CA1idx),'ro','MarkerSize',8,'LineWidth',1.5);
%     plot(allspikewidth(iCA1idx), allmeanrate(iCA1idx),'bx','MarkerSize',8,'LineWidth',1.5);
% 
%     xaxis = min(allspikewidth):0.1:max(allspikewidth);
%     plot(xaxis,0.1*ones(size(xaxis)),'k-','Linewidth',1);
%     plot(xaxis,7*ones(size(xaxis)),'k-','Linewidth',1);   
%     
%     set(gca,'YLim',[0 1.05*max(allmeanrate)]);
%     set(gca,'XLim',[0.95*min(xaxis) 1.05*max(xaxis)]);
%     
%     xlabel('Spike Widths','FontSize',24,'Fontweight','bold');
%     ylabel('Firing Rates (Hz)','FontSize',24,'Fontweight','bold');
%     title('Pyr Cells - dorsal and int CA1','FontSize',24,'Fontweight','bold')
%     
    
    
    
    

    %n) Partial vs defined in all
%     
%     figure; hold on;
%     redimscreen_figforppt1;
%     plot(allspikewidth(Pyridx), allmeanrate(Pyridx),'ro','MarkerSize',8,'LineWidth',1.5);
%     plot(allspikewidth(Partidx), allmeanrate(Partidx),'bx','MarkerSize',8,'LineWidth',1.5);
% 
%     xaxis = min(allspikewidth):0.1:max(allspikewidth);
%     plot(xaxis,0.1*ones(size(xaxis)),'k-','Linewidth',1);
%     plot(xaxis,7*ones(size(xaxis)),'k-','Linewidth',1);   
%     
%     set(gca,'YLim',[0 1.05*max(allmeanrate)]);
%     set(gca,'XLim',[0.95*min(xaxis) 1.05*max(xaxis)]);
%     
%     xlabel('Spike Widths','FontSize',24,'Fontweight','bold');
%     ylabel('Firing Rates (Hz)','FontSize',24,'Fontweight','bold');
%     title('Pyr Cells - All vs partial epochs','FontSize',24,'Fontweight','bold')
    
    
end % end figopt1
    
    

    
    

          

