function plotspikeraster(animaldir,day,epoch,tet)
% Plots spike rasters for the chosen epoch.
% animaldir - name of directory containing processed data
% day - day of experiment
% epoch - epoch to be analysed
% tet - tetrode to be analysed
% remember to check velocity cutoff (velocutoff) and time boundries (bnd)

% get the current directory to go back to after the function is done
currentdir = pwd;

% ---- File loading ----
% See if day number needs 0
dsz = '';
   if (day < 10)
      dsz = '0';
   end
   
% Specify data directory and load the file 
animalname = animaldir(1:end-1);
datadir = '/data14/jai/';

% Converts the day and epoch into text
dayt = num2str(day);
epocht = num2str(epoch);
tett = num2str(tet);
sfilename = strcat(datadir,animaldir,'/',animalname,'spikes',dsz,dayt,'.mat');

switch epoch
    
    case {1,3,5,7}
        disp('Sleep session, no data analysed');
        
    case 2
        lname = strcat(datadir,animaldir,'/',animalname,'linpos',dsz,dayt,'_lineartrack.mat');
        eval(['load ', lname]);   
    
    case 6
        lname = strcat(datadir,animaldir,'/',animalname,'linpos',dsz,dayt,'_lineartrack.mat');
        eval(['load ', lname]); 
        
    case 4
        lname = strcat(datadir,animaldir,'/',animalname,'linpos',dsz,dayt,'_wtrack.mat');
        eval(['load ', lname]); 
end    

% Load the mat file
eval(['load ', sfilename]);

  
% ----Data extraction ----

% Extract spike data from spike file

spkmat = [spikes{day}{epoch}{tet}{1}.data(:,1) spikes{day}{epoch}{tet}{1}.data(:,6) spikes{day}{epoch}{tet}{1}.data(:,3)];

% expand time boundries on either side of time cutoff
bnd = 0

% Extract position file from linpos file

switch epoch
    
    case {2,6}
    posmat = [linpos{day}{epoch}.statematrix.time floor(linpos{day}{epoch}.statematrix.linearVelocity)];
    pos2 = find(posmat(:,2));
    posmat = [posmat(pos2,1) posmat(pos2,2)]; 
    well = linpos{day}{epoch}.wellSegmentInfo.wellCoord(:,2);

    m = size(posmat,1);
    runbound = [min(posmat(:,1))];
    result = zeros(2,1);

    for i= 1:m-1
        for j = i+1

            if posmat(i,2)*posmat(j,2)<0
              result(1,1) = posmat(i,1);
              result(2,1) = posmat(j,1);  
              runbound = [runbound;result];
            end
        i=i+1;

        end     
    end

    m=size(runbound,1);
    % Stores [runstart runend runvelocity]
    runpeak=[];
    tmptb=[];
    result=zeros(1,3);

    for i= 1:m-1
        for j=i+1
            [lower]=find(posmat(:,1) == runbound(i,1));
            [higher]=find(posmat(:,1) == runbound(j,1));
            results(1,1)=runbound(i,1);
            results(1,2)=runbound(j,1);
            tmptb = [posmat(lower:higher,2)];
            tmptb = abs(tmptb);
            results(1,3) = max(tmptb)*(sign(posmat(lower,2)));
            runpeak = [runpeak;results];
            i=i+1;
        end
    end

    % removes runs where velocity is less than cutoff
    velocutoff = 10

    [lower] = find(abs(runpeak(:,3))>=velocutoff);
    runpeak = runpeak(lower,:,:);

    % get the spikes for each run

    spkplot=[];
    tmptb=[];

    figure;
    % Generates a title from information provided
    title(strcat('Spike raster plot for ',animalname,': Day ',dayt,' epoch ',epocht, ' Tetrode ', tett));
    xlabel('cm');
    ylabel('trial');
    hold on;
    
    
    for i=1:size(runpeak,1)-1;
        
        [list]=find(spkmat(:,1) >= runpeak(i,1)-bnd & spkmat(:,1) <= runpeak(i,2)+bnd);
        tmpindex = zeros(size(list));
        tmpindex = tmpindex + i;
        tmptb=[spkmat(list,1)-runpeak(i,1) tmpindex spkmat(list,3) spkmat(list,2)];
        spkplot = [spkplot;tmptb];
        for k = 1:length(tmptb(:,1))
            line([tmptb(k,3) tmptb(k,3)], [i-1 i], [tmptb(k,4) tmptb(k,4)], 'LineWidth', 0.5, 'Color', 'k');
        end

        ylim([0 i]);
        line([well(1,1) well(1,1)], [0 i], 'LineWidth', 2, 'Color', 'r');
        line([well(2,1) well(2,1)], [0 i], 'LineWidth', 2, 'Color', 'r');
        
        i=i+1;

    end 
        % ----Saving----
        % Saves figure as pdf
        % First checks if a folder called Plot exists in the processed data folder,
        % if not, creates it.

        cd(strcat(datadir,animaldir,'/'));

        plotdir = dir('Plot');

        if (isempty(plotdir))
            %an a plot folder needs to be created
            !mkdir Plot
        end

        % change to that directory and saves the figure with file name
        % animal_day_epoch
        cd(strcat(datadir,animaldir,'/','Plot/'));
        figurename = strcat(animalname,'_ras_',dayt,'e',epocht,'t',tett);
        saveas(gcf, figurename, 'pdf');
        
        % Closes the figure
        close;
    
    case 4
    
    for u=1:3 
        
        posmat = [linpos{day}{epoch}.statematrix.time floor(linpos{day}{epoch}.statematrix.linearVelocity)];
        pos2 = find(posmat(:,1+u));
        posmat = [posmat(pos2,1) posmat(pos2,u+1)]; 
        % not sure about well locations on W track
        %well = linpos{day}{epoch}.wellSegmentInfo.wellCoord(i,1);
        
        m = size(posmat,1);
        runbound = [min(posmat(:,1))];
        result = zeros(2,1);

        for i= 1:m-1
            for j = i+1

                if posmat(i,2)*posmat(j,2)<0
                  result(1,1) = posmat(i,1);
                  result(2,1) = posmat(j,1);  
                  runbound = [runbound;result];
                end
            i=i+1;

            end     
        end

        m=size(runbound,1);
        % Stores [runstart runend runvelocity]
        runpeak=[];
        tmptb=[];
        result=zeros(1,3);

        for i= 1:m-1
            for j=i+1
                [lower]=find(posmat(:,1) == runbound(i,1));
                [higher]=find(posmat(:,1) == runbound(j,1));
                results(1,1)=runbound(i,1);
                results(1,2)=runbound(j,1);
                tmptb = [posmat(lower:higher,2)];
                tmptb = abs(tmptb);
                results(1,3) = max(tmptb)*(sign(posmat(lower,2)));
                runpeak = [runpeak;results];
                i=i+1;
            end
        end

        % removes runs where velocity is less than cutoff
        velocutoff = 5

        [lower] = find(abs(runpeak(:,3))>=velocutoff);
        runpeak = runpeak(lower,:,:);

        % get the spikes for each run

        spkplot=[];
        tmptb=[];

        figure;
        % Generates a title from information provided
        ut = num2str(u);
        title(strcat('Spike raster plot for ',animalname,': Day ',dayt,' epoch ',epocht, ' Tetrode ', tett,'well ',ut));
        
        xlabel('cm');
        ylabel('trial');
        
        hold on;
        
        for i=1:size(runpeak,1)-1;
            [list]=find(spkmat(:,1) >= runpeak(i,1)-bnd & spkmat(:,1) <= runpeak(i,2)+bnd);
            tmpindex = zeros(size(list));
            tmpindex = tmpindex + i;
            tmptb=[spkmat(list,1)-runpeak(i,1) tmpindex spkmat(list,3) spkmat(list,2)];
            spkplot = [spkplot;tmptb];
            for k = 1:length(tmptb(:,1))
                line([tmptb(k,3) tmptb(k,3)], [i-1 i], [tmptb(k,4) tmptb(k,4)], 'LineWidth', 0.5, 'Color', 'k');
            end

            ylim([0 i])
            xlim([0 max(tmptb(:,3))+5]);
            %lim([0 i]);
            %line([well(1,1) well(1,1)], [0 i], 'LineWidth', 2, 'Color', 'r');
            %line([well(2,1) well(2,1)], [0 i], 'LineWidth', 2, 'Color', 'r');

            i=i+1;
        end
            % ----Saving----
            % Saves figure as pdf
            % First checks if a folder called Plot exists in the processed data folder,
            % if not, creates it.

            cd(strcat(datadir,animaldir,'/'));

            plotdir = dir('Plot');

            if (isempty(plotdir))
                %an a plot folder needs to be created
                !mkdir Plot
            end

            % change to that directory and saves the figure with file name
            % animal_day_epoch
            cd(strcat(datadir,animaldir,'/','Plot/'));
            figurename = strcat(animalname,'_ras_',dayt,'e',epocht,'t',tett,'w',ut);
            saveas(gcf, figurename, 'pdf');

            % Closes the figure
            close;
        u+1;
    end
end

cd(currentdir);





% PSTH

%figure;
%edges = [0:1:60];

%psth = histc(spkplot(:,3),edges);

%bar(edges,psth);
%xlim([0 60]);
%xlabel('cm');
%ylabel('spikes');







