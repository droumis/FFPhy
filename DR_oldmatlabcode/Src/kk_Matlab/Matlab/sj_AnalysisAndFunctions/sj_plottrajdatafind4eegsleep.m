function out = sj_plottrajdatafind4eegsleep(animdirect,prefix,index, varargin)

% sj_plottrajdatafind4eegsleep('/data25/sjadhav/RippleInterruption/REf_direct','REf',[9 5]);
% EEG  + Ripple band during sleep with spped. Kinda like sj_plottrajdatafind4eeg, but without parsing by trajectory

% 5: Dont go by trajectory - Just rolling times with 5 secs
% 4: Only look for eeg examples - forget spikes

% Do not use subplots. All on one figure and put linearized position on top
% of rasters. 1 eeg and ripple band above. And Speed above that


%Shantanu - CHange version 1.
% Look for only 1 kind of traj. Pick cells that you want. and change
% subplot dimensions

% Shantanu
% Load linpos data for day and epoch. Get traj start sequences and plot a
% giant figure for each trajectory / pair of trajectories in a loop so that
% you can pick example figures. You want to plot linearized position,
% linearvelocity/abs veloctiy (need pos field for that), place cell spike
% rasters (spike field), and LFP - either EEG or ripple band (eeg/ripple file),
% and DIO times marked (DIO structure - jbfill it).

% Can also load ripples structure to superimpose. Ideally, I should have a
% file which has looked for ripples combined across tetrodes

% After each figure, stop in a loop

day=index(1,1);
epoch = index(1,2);
out = 1;

Fspos = 30; %Hz
respos = 1/30; % sec
Fseeg = 1500; %Hz
reseeg = 1/1500; % sec
Fsspikes = 10000; %Hz
resspikes = 1/10000; %sec



% Use subset of eegtets for eeg and ripple: 2 each
tets = [1,5,9,10,11,12];
eegtets = [5, 9, 10, 12];
maineegtets = [9, 5];
maineegidxs = [2];


% Cell indices within spikes structure - can also get by doing
% setfiltercellswith the cellinfo file

% Trajectories:
% 1 = Out Left    well 1to2
% 2 = In Left     well 2to1
% 3 = Out Right   well 1to3
% 4 = In Right    well 3to1


posfile = sprintf('%s/%spos%02d.mat', animdirect, prefix, day);
load(posfile);
absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
posn = pos{day}{epoch}.data(:,2:3); % Can also use fields 6:7
postime = pos{day}{epoch}.data(:,1); % same as linpostime

DIOfile = sprintf('%s/%sDIO%02d.mat', animdirect, prefix, day);
load(DIOfile);
% Get DIO times
stim = DIO{day}{epoch}{16};
if isempty(stim)
    stim = DIO{day}{epoch}{15};
end
stim_starttime = stim.pulsetimes(:,1)./10000; %sec
stim_endtime = stim.pulsetimes(:,2)./10000; %sec


% EEG and if desired ripple band LFP data + theta band LFP data
cnt=0;
for tet=eegtets
    cnt=cnt+1;
    eegfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
    load(eegfile);
    eval(['eegdata',num2str(cnt),'= eeg{day}{epoch}{tet}.data;'])
    ripplefile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
    load(ripplefile);
    eval(['rippledata',num2str(cnt),'= ripple{day}{epoch}{tet}.data;'])
    thetafile = sprintf('%s/EEG/%stheta%02d-%01d-%02d.mat', animdirect, prefix, day,epoch,tet);
    load(thetafile);
    eval(['thetadata',num2str(cnt),'= theta{day}{epoch}{tet}.data;'])
    if cnt==1
        teeg = geteegtimes(eeg{day}{epoch}{tet});
        eegstart = eeg{day}{epoch}{tet}.starttime;
        eegsamprate = round(eeg{day}{epoch}{tet}.samprate);
        % Theta is downsampled 10 times
        ttheta = geteegtimes(theta{day}{epoch}{tet});
        thetastart = theta{day}{epoch}{tet}.starttime;
        thetasamprate = round(theta{day}{epoch}{tet}.samprate);
    end
end

%--------------------------------------------------------------------------

speedthrs1 = 2; %cm/sec;
% Set fig prop
set(0,'defaultaxesfontsize',12);set(0,'defaultaxesfontweight','normal');
set(0,'defaultaxeslinewidth',1);
figno = 1;

win = 4; % x sec windows with win/2 sec overlap
currstart = 1;
currend = lookup(postime(currstart)+win,postime);

% Plot in a loop
while currend < length(postime)
    
    starttime = postime(currstart);
    endtime = postime(currend);
    plotsttime = starttime;
    plotendtime = endtime;
    
    figure(figno); hold on;
    redimscreen_widehor;
    
    % Speed
    %--------
    % Speed ylimits in plot: 0-3
    
    figure(figno);
    hold on;
    currspeed = absvel(currstart:currend);
    currpostime = postime(currstart:currend);
    currpostime = currpostime - postime(currstart);
    
    % Change scale
    speedscale = max(currspeed)-min(currspeed);
    plotscale = 3-0;
    currspeed = currspeed.*(plotscale/speedscale);
    
    plot(currpostime, currspeed,'Linewidth',2);
    set(gca,'XLim',[0 max(currpostime)]);
    
    % Line for speedthrs
    yaxis = speedthrs1*(plotscale/speedscale)*ones(size([0:length(currspeed)]'));
    xaxis = [0:length(currspeed)]';
    plot(xaxis,yaxis,'r');
    speedthrscale = speedthrs1*(plotscale/speedscale);
    
    y1=min(currspeed); y2=max(currspeed);
    % Plot 1 sec scale bar
    xaxis = 0:0.1:1;
    yaxis=(y2)*0.9*ones(size(xaxis));
    plot(xaxis,yaxis,'k-','Linewidth',4);
    
    
    % 2) 2 EEGs
    baseline = 4; %
    
    % Get EEG Ind in current traj
    eind1 = lookup(starttime, teeg);
    eind2 = lookup(endtime, teeg);
    
    eegtimeaxis = teeg(eind1:eind2);
    eegtimeaxis = eegtimeaxis - teeg(eind1);
    
    tind1 = lookup(starttime, ttheta);
    tind2 = lookup(endtime, ttheta);
    thetatimeaxis = ttheta(tind1:tind2);
    thetatimeaxis = thetatimeaxis - ttheta(tind1);
    
    figure(figno); hold on;
    for i=1:length(maineegidxs)
        %                 % Get theta and ripple
        if i==1
            % Theta
%             eval(['thetadata = thetadata',num2str(maineegidxs(i)),';']);
%             curreegdata = double(thetadata(tind1:tind2,1));
%             eegscale = max(curreegdata)-min(curreegdata);
%             downeeg = baseline+1;
%             upeeg = downeeg+5;
%             plotscale = 5;
%             curreegdata = downeeg + (plotscale/2) + curreegdata.*(plotscale/eegscale);
%             plot(thetatimeaxis,curreegdata,'Linewidth',1);
%             
%             % Update baseline
%             baseline = upeeg;
            
            % Ripple
            eval(['rippledata = rippledata',num2str(maineegidxs(i)),';']);
            curreegdata = double(rippledata(eind1:eind2,1));
            eegscale = max(curreegdata)-min(curreegdata);
            downeeg = baseline;
            upeeg = downeeg+8;
            plotscale = 8;
            curreegdata = downeeg + (plotscale/2) + curreegdata.*(plotscale/eegscale);
            plot(eegtimeaxis,curreegdata,'k-','Linewidth',1);
            
            % Update baseline
            baseline = upeeg;
        end
        eval(['eegdata = eegdata',num2str(maineegidxs(i)),';']);
        curreegdata = eegdata(eind1:eind2);
        eegscale = max(curreegdata)-min(curreegdata);
        downeeg = baseline+1;
        upeeg = downeeg+12;
        plotscale = 12;
        curreegdata = downeeg + (plotscale/2) + curreegdata.*(plotscale/eegscale);
        plot(eegtimeaxis,curreegdata,'k-','Linewidth',1);
        
        % Update baseline
        baseline = upeeg;
        
        
        
    end
    
    % Return to user
    keyboard;
    close(figno);
    
    
    % Update Window Indices for next iteration
    
    currstart = lookup(starttime+win/2,postime);
    currend = lookup(postime(currstart)+win,postime);
    
    
end







