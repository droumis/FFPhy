% use this script to visually check theta epochs extracted by
% gethighthetatimes -- also use to empirically determine good values for
% powerratio

% 1st block: checks for BIMODALITY
% 2nd block: plots eeg, theta, powerratios, velocity for a chosen tetrode
% 3rd block: plots high theta times for all tetrodes in animal + ripple
% times


animaldir = '/data25/sjadhav/HPExpt/HPa_direct/';
animalprefix = 'HPa';
epochs = [4 4]; % this is day and epochs
tetlist = 1;

% set these to taste
powerratio1 = 4;             %% set manually so you can see
powerratio2 = .3;
velocitythresh = 3;
windowsize = 10;            % in seconds
smoothing_width = .25;      % for smoothing of powerratio
mindur = 0.5;
ntheta_thresh = 3;
consensus = 0;            % set to 1 (see below) if want to plot same acceptable period for all
% set to 0 if plot each individual tetrode's
% acceptable period

doplot=1; noplot=0;


eeg = loadeegstruct(animaldir, animalprefix, 'eeg', epochs(1), epochs(2),tetlist);
theta = loadeegstruct(animaldir, animalprefix, 'theta', epochs(1), epochs(2),tetlist);
delta = loadeegstruct(animaldir, animalprefix, 'delta', epochs(1), epochs(2),tetlist);
supratheta = loadeegstruct(animaldir, animalprefix, 'supratheta', epochs(1), epochs(2),tetlist);
pos = loaddatastruct(animaldir, animalprefix, 'pos', epochs(1));
eeg = eeg{epochs(1)}{epochs(2)}{tetlist};
theta = theta{epochs(1)}{epochs(2)}{tetlist};
supratheta = supratheta{epochs(1)}{epochs(2)}{tetlist};
delta = delta{epochs(1)}{epochs(2)}{tetlist};
tdratio = double(theta.data(:,3))./double(delta.data(:,3));
tsratio = double(theta.data(:,3))./double(supratheta.data(:,3));
pos = pos{epochs(1)}{epochs(2)};
windowsize_possamp = windowsize*29.97;
windowsize_thetasamp = floor(windowsize*theta.samprate);
windowsize_eegsamp = windowsize*1500;
nowindows = floor((pos.data(end,1)-pos.data(1,1))/windowsize);
Fs = round(eeg.samprate);
timevec = 0:(1/Fs):windowsize;

%Check for bimodality in theta, theta:delta, and velocity.

if noplot %0
    velocities = [];
    powers = [];
    ratios = [];
    
    for w = 1:nowindows
        starttime = pos.data(1,1) + (w-1)*windowsize;
        endtime = starttime + windowsize;
        %pos
        posindices = (lookup(starttime,pos.data(:,1))):(lookup(endtime,pos.data(:,1)));
        velocities = [ velocities ; mean(pos.data(posindices,5)) ];
        %theta
        thetatimes = geteegtimes(theta);
        thetaindices = (lookup(starttime,thetatimes)):(lookup(endtime,thetatimes));
        powers = [ powers ; double(mean(theta.data(thetaindices,3)))] ;
        %tdratio
        ratios = [ ratios ; double(mean(tdratio(thetaindices)))] ;
    end
    
    tempratios = ratios;
    tempratios(find(ratios==Inf))= 2;
    
    figure
    hist(velocities,100)
    title('velocities')
    figure
    hist(powers,100)
    title('powers')
    figure
    hist(tempratios,100)
    title('theta:delta ratios')
end


% Plot raw eeg, theta, delta, tdr, and velocity in scrolling window.
if noplot %0
    
    % find high theta times
    out = gethighthetatimes2(animaldir,animalprefix, epochs, tetlist, 'powerratio1',powerratio1,'powerratio2',powerratio2, ...
        'velocitythresh',velocitythresh,'mindur',mindur);
    out = out{epochs(1)}{epochs(2)};
    eegtimes=geteegtimes(eeg);
    thetaperiods=interp1(out.time,out.ntheta,eegtimes,'nearest');
    highthetatrace = eeg.data.*thetaperiods';
    highthetatrace(highthetatrace==0)=NaN;
    
    %smooth power ratio
    samprate = 150;
    kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
    tsratio_smooth = smoothvect(tsratio, kernel);
    tdratio_smooth = smoothvect(tdratio, kernel);
    
    for w = 1:nowindows
        starttime = pos.data(1,1) + (w-1)*windowsize;
        endtime = starttime + windowsize;
        %h = figure;
        figure('units','normalized','outerposition',[0 0 1 1])
        
        %pos
        subplot(5,1,5)
        posindices = (lookup(starttime,pos.data(:,1))):(lookup(endtime,pos.data(:,1)));
        velocity = pos.data(posindices,5);
        velocity(isnan(velocity))=0;
        threshindices = velocity > velocitythresh;
        hold on
        plot(velocity,'Color',[.8 .8 .8],'LineWidth',5)
        plot(velocity.*threshindices,'k.','MarkerSize',15)
        plot(velocitythresh*ones(1,length(velocity)),'k--')
        title('velocity')
        ylim([0 20])
        axis tight;
        
        %theta, delta
        subplot(10,1,3:4)
        thetatimes = geteegtimes(theta);
        thetaindices = (lookup(starttime,thetatimes)):(lookup(endtime,thetatimes));
        thetatrace = theta.data(thetaindices,1);
        deltatrace = delta.data(thetaindices,1);
        suprathetatrace = supratheta.data(thetaindices,1);
        hold on
        plot(thetatrace,'r','LineWidth',2)
        plot(deltatrace,'g','LineWidth',2)
        plot(suprathetatrace,'b','LineWidth',2)
        
        axis tight
        ylim([-700 700])
        
        %theta:supratheta ratio
        subplot(10,1,1)
        hold on
        plot(tsratio(thetaindices),'Color',[.9 .9 .9],'LineWidth',2)
        plot(tsratio_smooth(thetaindices),'k','LineWidth',2)
        plot(powerratio1*ones(1,length(thetaindices)),'b--','LineWidth',2)
        title('supratheta')
        axis tight
        ylim([0 30])
        
        %theta:delta ratio
        subplot(10,1,2)
        hold on
        plot(tdratio(thetaindices),'Color',[.9 .9 .9],'LineWidth',2)
        plot(tdratio_smooth(thetaindices),'k','LineWidth',2)
        plot(powerratio2*ones(1,length(thetaindices)),'g--','LineWidth',2)
        title('delta')
        axis tight
        ylim([0 4])
        
        %raweeg
        subplot(5,1,3:4)
        hold on
        eegtimes = geteegtimes(eeg);
        eegstartindex = lookup(starttime,eegtimes);
        eegtrace = eeg.data(eegstartindex:(eegstartindex+windowsize_eegsamp));
        plot(timevec,eegtrace,'k','LineWidth',2)
        ylim([-3000 3000])
        
        %highthetatrace
        plot(timevec,highthetatrace(eegstartindex:(eegstartindex+windowsize_eegsamp)),'b','LineWidth',2)
        
        pause
        
        close all
    end
end



% Compare high theta between tetrodes / regions + plot ripple times.

if doplot %1
    clear out; clear highthetatrace; clear eeg; clear rippletrace;
    tetlist = [1:30];
    tetinfo = loaddatastruct(animaldir, animalprefix,'tetinfo');
    ripples = loaddatastruct(animaldir, animalprefix,'ripples',epochs(1));
    % extract recording areas of each tetrode and determine whether valid to plot
    validtet = zeros(size(tetlist));
    tetlist_region = cell(size(tetlist));
    for tet=1:length(tetlist)
        try
            tetlist_region{tet} = tetinfo{epochs(1)}{epochs(2)}{tetlist(tet)}.area;            
            %if strcmp(tetlist_region{tet},'Reference') || strcmp(tetlist_region{tet},'cc') || strcmp(tetlist_region{tet},'ctx')
            
            if isfield(tetinfo{epochs(1)}{epochs(2)}{tetlist(tet)},'descrip')
                tetlist_descrip{tet} = tetinfo{epochs(1)}{epochs(2)}{tetlist(tet)}.descrip;
                if strcmp(tetlist_descrip{tet},'CA1Ref') || strcmp(tetlist_descrip{tet},'iCA1Ref') || strcmp(tetlist_descrip{tet},'ctx')
                    validtet(tet) = 1;
                    disp(sprintf('tetrode %d : reference',tetlist(tet)))
                    tetlist_region{tet} = 'CA1Ref'; % Mark region as CA1Ref for reference tetrode
                    continue
                end  
            end
            
            try
                nocellsontet = tetinfo{epochs(1)}{epochs(2)}{tetlist(tet)}.numcells;
                if nocellsontet >= 2         % at least 2 cells
                    validtet(tet) = 1;
                    disp(sprintf('tetrode %d : %d cells',tetlist(tet),nocellsontet))
                    continue
                else
                    disp(sprintf('tetrode %d : ignored, %d cells',tetlist(tet),nocellsontet))
                    continue
                end
            catch
                disp(sprintf('tetrode %d : no .numcells field',tetlist(tet)))
                continue
            end
        catch
            disp(sprintf('tetrode %d : no area field / other error',tetlist(tet)))
            continue
        end
    end
    
    
    % load full eeg + calculate high theta times + load ripple times
    for tet=1:length(tetlist)
        
        if validtet(tet)
            % eeg
            eeg{tet} = loadeegstruct(animaldir, animalprefix, 'eeg', epochs(1), epochs(2),tetlist(tet));
            eeg{tet} = eeg{tet}{epochs(1)}{epochs(2)}{tetlist(tet)};
            % high theta times
            
            if consensus == 1
                out{tet} = gethighthetatimes2(animaldir,animalprefix, epochs, [],'tetfilter','(isequal($area, ''CA1'') && ($numcells >= 2))', 'powerratio1',powerratio1,'powerratio2',powerratio2, ...
                    'mindur',mindur);
            else
                out{tet} = gethighthetatimes2(animaldir,animalprefix, epochs, tetlist(tet), 'powerratio1',powerratio1,'powerratio2',powerratio2, ...
                    'mindur',mindur);
            end
            
            out{tet}=out{tet}{epochs(1)}{epochs(2)};
            eegtimes=geteegtimes(eeg{tet});
            thetaperiods{tet}=interp1(out{tet}.time,out{tet}.ntheta,eegtimes,'nearest')';
            
            %filter for periods that exceed ntheta_thresh
            if range(thetaperiods{tet}) > 1       %  if ntheta is more than 0s and 1s
                highthetatrace{tet} = eeg{tet}.data.*(thetaperiods{tet} >= ntheta_thresh);
                highthetatrace{tet}(highthetatrace{tet} == 0)=NaN;
            else
                highthetatrace{tet} = eeg{tet}.data.*(thetaperiods{tet} == 1);
                highthetatrace{tet}(highthetatrace{tet} == 0)=NaN;
            end
            
            % ripple times
            rippletrace{tet} = nan(size(eegtimes))';
            noripples = length(ripples{epochs(1)}{epochs(2)}{tetlist(tet)}.startind);
            disp(sprintf('tetrode %d : %d cells, %d ripples',tetlist(tet),nocellsontet,noripples))
            
            for j=1:noripples
                startind(j) = ripples{epochs(1)}{epochs(2)}{tetlist(tet)}.startind(j);
                endind(j) = ripples{epochs(1)}{epochs(2)}{tetlist(tet)}.endind(j);
                rippletrace{tet}(startind(j):endind(j)) = 1;
            end
        end
    end
    
    for w = 1:nowindows
        starttime = pos.data(1,1) + (w-1)*windowsize;
        endtime = starttime + windowsize;
        %h = figure;
        figure('units','normalized','outerposition',[0 0 1 1])
        
        %pos
        subplot(6,1,6)
        posindices = (lookup(starttime,pos.data(:,1))):(lookup(endtime,pos.data(:,1)));
        velocity = pos.data(posindices,5);
        threshindices = velocity > velocitythresh;
        hold on
        plot(velocity,'Color',[.8 .8 .8],'LineWidth',5)
        plot(velocity.*threshindices,'k.','MarkerSize',15)
        plot(velocitythresh*ones(1,length(velocity)),'k--')
        title('velocity')
        ylim([0 20])
        axis tight;
        
        %raweeg
        subplot(6,1,1:5)
        hold on
        
        plotcount=0;
        %regions = {'cc' 'ctx' 'Reference' 'CA1' 'CA2' 'CA3'};
        regions = {'CA1' 'iCA1','CA1Ref','iCA1Ref'};
        tetrodes_plotted = [];
        
        for r=1:4     % iterate through regions
            region=regions{r};
            for tet=[find(validtet==1)]
                if strcmp(tetlist_region{tet},region)
                    tetrodes_plotted = [tetrodes_plotted tet];
                    eegtimes = geteegtimes(eeg{tet});
                    eegstartindex = lookup(starttime,eegtimes);
                    eegtrace = eeg{tet}.data(eegstartindex:(eegstartindex+windowsize_eegsamp));
                    %plot
                    %select colors
                    switch region
                        case 'CA1'
                            clr1 = [.9 .9 .9];
                            clr2 = [0 0 0];
                        case 'CA1Ref'
                            clr1 = [1 .9 .9];
                            clr2 = [1 0 0];
                        case 'iCA1'
                            clr1 = [0 .7 0];
                            clr2 = [.7 1 .7];
                        case 'iCA1Ref'
                            clr1 = [1 .9 1];
                            clr2 = [1 0 1];
                        case 'cc'
                            clr1 = [1 .9 1];
                            clr2 = [1 0 1];
                        case 'ctx'
                            clr1 = [1 .9 1];
                            clr2 = [1 0 1];
                    end
                    %eeg
                    plot(timevec,eegtrace-plotcount*1000,'Color',clr1,'LineWidth',1.2)
                    %highthetatrace
                    plot(timevec,highthetatrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000,'Color',clr2,'LineWidth',1.2)
                    %rippletrace
                    plot(timevec,rippletrace{tet}(eegstartindex:(eegstartindex+windowsize_eegsamp))-plotcount*1000+500,'b','LineWidth',5)
                    
                    plotcount=plotcount+1;
                end
            end
        end
        ylim([-plotcount*1000 400])
        line1 = sprintf('%s highthetatimes day %d, epoch %d',animalprefix,epochs(1),epochs(2));
        line2 = sprintf('powerratio-supra: %d, powerratio-theta: %d',powerratio1,powerratio2);
        line3 = sprintf('tetrodes: %s',mat2str(tetrodes_plotted));
        title({line1 line2 line3},'FontSize',14,'FontWeight','bold')
        pause
        close all
    end
    
end

