function Plot_MU_veloc_rate(animaldir,day,binsize)

% animal dir - name of the directory containing processed data
% day - experimental day to be analysed
% epoch - experimental epoch to be analysed
% tet - tetrode number
%
% Uses the same start as Plotopenfieldrate.m

% get the current directory to go back to after the function is done
currentdir = pwd

for tet=1:8;
    for epoch = 2:2:6;
        epochstep=2;
        % ---- File loading ----
        % See if day number needs 0
        dsz = '';
           if (day < 10)
              dsz = '0';
           end

        tsz = '';
           if (tet < 10)
              tsz = '0';
           end


        % Specify data directory and load the file 
        animalname = animaldir(1:end-1);
        datadir = '/data14/jai/';

        % Converts the day and epoch into text
        dayt = num2str(day);
        epocht = num2str(epoch);
        tett = num2str(tet);

        % Loads the position data
        posfilename = strcat(datadir,animaldir,'/',animalname,'pos',dsz,dayt,'.mat');
        eval(['load ', posfilename]);

        % Loads the spike data

        % Specify data directory and load the file 
        animalname = animaldir(1:end-1);
        datadir = '/data14/jai/';
        dr = [datadir,animalname,'/',animalname,dsz,dayt,'/']
        cd(dr); 
        f = dir(strcat(tsz,tett,'*'));

        mufile = sprintf('%s/%s%s%s-%s%s_params.mat',f(1).name,animalname,dsz,dayt,tsz,tett);
        eval(['load ', mufile]);

        % Load posfile

        % Loads the position data
        posfilename = strcat(datadir,animaldir,'/',animalname,'pos',dsz,dayt,'.mat');
        eval(['load ', posfilename]);

        % make table with time filter properties
        %epochtb = [];
        %tmf = double(unique(clustdata.timefilters));
        %epochtb(:,1) = tmf(3:end,1);

        % plot mu epoch above threshold locations
 
        tma = clustdata.timefilterranges(epoch+1,1);
        tmb = clustdata.timefilterranges(epoch+1,2);
        [tmplist] = find(clustdata.params(:,1) >= tma  & clustdata.params(:,1) <= tmb );
        %tmptbl(:,1) = clustdata.params(tmplist,8);
        %tmptbl(:,2) = clustdata.params(tmplist,9);
        timepoints = clustdata.params(tmplist,1);

        % Get spiketimes
        [spikepos, posindex] = lookuptime3(timepoints/10000, pos{day}{epoch}.data(:,1),pos{day}{epoch}.data(:,2:4)');

        findgoodpoints = find((spikepos(1,:)~=0)&(spikepos(2,:)~=0));
        spikepos = spikepos(:,findgoodpoints)';
        timepoints = timepoints(findgoodpoints',1);
        posindex = posindex(findgoodpoints);


        % ---- Extract relevant data ----
        % Extracts spike time position data from spike file
        posx = spikepos(:,1);
        posy = spikepos(:,2);
        posindexfield = posindex;
        time = timepoints;

        if ~isempty(spikepos)
            spikes = [time posx posy posindexfield];

        else
            spikes = [0 0 -1];
        end

        % Extracts position data from pos file
        pospos = [pos{day}{epoch}.data(:,1) pos{day}{epoch}.data(:,2) pos{day}{epoch}.data(:,3)];

        timestep = 0.01;
        %timesetp = pospos(2,1) - pospos(1,1);

        goodspikes = [];
        goodspikes = spikes;
        %make a cell array, where each cell contains data for one trajectory.
        %inside each cell [binlocation occupancy spikecount firingrate]

        tmpposition = (pospos(:,[2 3]));
        tmpspikes = (goodspikes(:,[2 3]));

        % ---- Calculations ----
        % Calculate firing rate and spatial occupancy

        if ~isempty(tmpposition)
            minx = floor(min(tmpposition(:,1)));
            maxx = ceil(max(tmpposition(:,1)));
            binx = (minx:binsize:maxx);
            miny = floor(min(tmpposition(:,2)));
            maxy = ceil(max(tmpposition(:,2)));
            biny = (miny:binsize:maxy);

            [output.occupancy output.xticks output.yticks] = HIST2(tmpposition(:,1),tmpposition(:,2), binx, biny);

            nonzero = find(output.occupancy ~= 0);
            [output.spikes BX BY] = hist2(tmpspikes(:,1), tmpspikes(:,2), binx, biny);
            output.spikerate = zeros(size(output.spikes));
            output.spikerate(nonzero) = output.spikes(nonzero) ./(timestep* output.occupancy(nonzero) );

            % Calculate smoothed occupancy
            %g = gaussian2(std,(6*std));
            %output.smoothedspikerate = filter2(g,(output.spikerate)); % is this the right filter?
            %smoothedoccupancy = [];
            %smoothedoccupancy = zeros(size(output.spikes));
            %smoothedoccupancy = filter2(g, output.occupancy);
            %zero = find(smoothedoccupancy == 0);

            %output.smoothedspikerate(zero) = -1;
        end


        % ----Data extraction ----

        % Extract position data from day processed file

        xdata = pos{day}{epoch}.data(:,2);
        ydata = pos{day}{epoch}.data(:,3);
        veldata = pos{day}{epoch}.data(:,5);

        % Get the boundaries of data values
        xmin = min(xdata);
        xmax = max(xdata);

        ymin = min(ydata);
        ymax = max(ydata);


        %posdata = [xdata-xmin,ydata-ymin];

        % Make matrix with intervals of 1cm

        gridint = 1

        %xedges = [floor(xmin):gridint:ceil(xmax)];
        %yedges = [floor(ymin):gridint:ceil(ymax)];

        xedges = [floor(xmin):gridint:ceil(xmax)];
        yedges = [floor(ymin):gridint:ceil(ymax)];

        histbin = [zeros(size(xdata),3)];
        [dum,histbin(:,1)] = histc(xdata,xedges);
        [dum,histbin(:,2)] = histc(ydata,yedges);

        % Put velocity in 3rd column

        histbin(:,3) = veldata(:); 

        % Averages all the velocity values for each coordinate

        %accumdata = accumarray([histbin(:,2) histbin(:,1)],histbin(:,3),[size(xedges,1) size(yedges,1)], @mean);
        accumdata = accumarray([histbin(:,1) histbin(:,2)],histbin(:,3),[size(xedges,2) size(yedges,2)], @mean, NaN);

        % Match mean velocity and spike rate

        velmat = flipud(rot90(accumdata,1));
        ratemat = output.spikerate;

        % Makes a 2 column table velocity spike rate
        corrtbl=[]
        tmptbl=[]

        for i=1:(size(velmat,2)-1)
            tmptbl(:,1) = velmat(:,i);
            tmptbl(:,2) = ratemat(:,i);
            corrtbl = [corrtbl;tmptbl];

            i=i+1;
        end

        [n] = isnan(corrtbl(:,1));
        [n] = find(n ==0);
        corrtbl = corrtbl(n,:);

        R = corrcoef(corrtbl(:,1), corrtbl(:,2));
        R = R(1,2);

        figure;

        plot(corrtbl(:,1), corrtbl(:,2),'k.');

        h = lsline;
        set(h, 'LineWidth',1, 'Color', 'r');

        txt = strcat('coeff. ', num2str(R));
        ylabel('Hz');
        xlabel('cm/s');
        % Generates a title from information provided
        title(['MU spike rate vs Velocity ',animalname,': Day ',dayt,' tetrode ', tett,' epoch ',epocht]);

        ax = axis;

        l = text(ax(1,2)*0.6,ax(1,4)*0.8,txt,'Color','r');


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
        figurename = strcat(animalname,'_MUspr_v_',dayt,'t',tett,'e',epocht);
        saveas(gcf, figurename, 'pdf');

        % Closes the figure
        close;
        
        epoch = epoch+epochstep;
    end
    
    tet = tet+1;
end


    % ----Change directory----
    % Returns to directory at the start
    cd(currentdir);
 
 
 
 
    

    


