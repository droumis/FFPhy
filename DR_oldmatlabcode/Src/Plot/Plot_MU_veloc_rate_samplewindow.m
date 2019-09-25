function Plot_MU_veloc_rate(animaldir,day,timestep)

% animal dir - name of the directory containing processed data
% day - experimental day to be analysed
% 
% timestep - time interval between location
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

        % Extracts position data from pos file [time x y velocity]
        pospos = [pos{day}{epoch}.data(:,1)*10000 pos{day}{epoch}.data(:,2) pos{day}{epoch}.data(:,3) pos{day}{epoch}.data(:,5)];

        
        % Time bins
        n = 1:timestep+1:size(pospos,1);
        tedge = [pospos(n,1); max(pospos(:,1))];
        
        
        % get velocity for time points
        % Mtable lists time velocity spikerate
        Mtable = [];
        histbin = [];
        
        Mtable(:,1) = tedge;
        Mtable(:,3) = histc(clustdata.params(:,1),tedge);
        [dum, histbin(:,1)] = histc(pospos(:,1), tedge);
        
        Mtable(:,2) = accumarray(histbin(:,1), pospos(:,4), size(tedge), @mean, NaN);
        
        Mtable(:,3) = Mtable(:,3)*10000/(pospos(1+timestep,1)-pospos(1,1));
        
        
        [R,P]  = corrcoef(Mtable(:,2), Mtable(:,3),'rows','complete');
        R = R(1,2);
          
        figure; 
        
        plot(Mtable(:,2), Mtable(:,3),'k.');

        h = lsline;
        set(h, 'LineWidth',1, 'Color', 'r');

        txt = strcat('coeff. ', num2str(R));
        ylabel('Hz');
        xlabel('cm/s');
        % Generates a title from information provided
        tbin = num2str((pospos(1+timestep,1)-pospos(1,1))/10);
        title(['MU spike rate vs Velocity ',animalname,': Day ',dayt,' tetrode ', tett,' epoch ',epocht,' ', tbin,' ms bin' ]);

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
        figurename = strcat(animalname,'_MUspr_v',dayt,'t',tett,'e',epocht);
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
 
 

 
    

    


