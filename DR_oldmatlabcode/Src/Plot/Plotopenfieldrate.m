function [output] = openfieldratemap(animaldir,tet,day,epoch, binsize,std,spn,spr,sspr,ot)

% animal dir - name of the directory containing processed data
% day - experimental day to be analysed
% epoch - experimental epoch to be analysed
% tet - tetrode number
%
%[output] = openfieldratemap(spikes,pos, binsize, std)
%[output] = openfieldratemap(spikes,pos, binsize)
%
%binsize- the length of each spatial bin (default 1 cm)
%std - defines the shape of the 2d gaussian used to smooth spikerate. (default 1)
%spn - plot spike number
%spr - plot spike rate
%sspr - plot occupancy normalised spike rate
%ot - plot outliers, default 0
%Calculates the 2d occupancy normalized firing rate for the cell
%Plots the data

%
%The output is a structure with n matrices. The matrices are: occupancy, 
%bin vector x, bin vector y, bin spike count, occ normailized firing per bin, and smoothed
%occ normalized firing. 
%

warning('OFF','MATLAB:divideByZero');

% get the current directory to go back to after the function is done
currentdir = pwd

% Set default parameters
if (nargin < 5)
    binsize = 1;
end

if (nargin < 6)
    std = 1;
end

if (nargin < 7)
    spn = 1;
end

if (nargin < 8)
    spr = 1;
end

if (nargin < 9)
    spn = 1;
end

if (nargin < 10)
    ot = 0;
end


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

% Loads the position data
posfilename = strcat(datadir,animaldir,'/',animalname,'pos',dsz,dayt,'.mat');
eval(['load ', posfilename]);

% Loads the spike data
spkfilename = strcat(datadir,animaldir,'/',animalname,'spikes',dsz,dayt,'.mat');
eval(['load ', spkfilename]);


% ---- Extract relevant data ----
% Extracts spike time position data from spike file
posx = spikes{day}{epoch}{tet}{1}.data(:,2);
posy = spikes{day}{epoch}{tet}{1}.data(:,3);
posindexfield = spikes{day}{epoch}{tet}{1}.data(:,7);
time = spikes{day}{epoch}{tet}{1}.data(:,1);

if ~isempty(spikes)
    spikes = [time posx posy posindexfield];

else
    spikes = [0 0 -1];
end

% Extracts position data from pos file
pospos = [pos{day}{epoch}.data(:,1) pos{day}{epoch}.data(:,2) pos{day}{epoch}.data(:,3)];

timestep = pospos(2,1) - pospos(1,1);

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
    g = gaussian2(std,(6*std));
    output.smoothedspikerate = filter2(g,(output.spikerate)); % is this the right filter?
    smoothedoccupancy = [];
    smoothedoccupancy = zeros(size(output.spikes));
    smoothedoccupancy = filter2(g, output.occupancy);
    zero = find(smoothedoccupancy == 0);

    output.smoothedspikerate(zero) = -1;
end
% ----Plot figure----


% Plot the firing rate
% Defines the colormap, with the lowest value set to grey
if spr == 1

    nc = 1024;
    cmap = jet(nc);
    cmap(1,:) = 0.8;
    colormap(cmap);

    % max color scale value
    amax = max(output.spikerate(:));
    tett = num2str(tet);

    % Creates an image of the data

    imagedata = flipud(output.spikerate);

    % Firing rate
    image(imagedata,'CDataMapping', 'scaled');

    % Generates a title from information provided
    title(strcat('Occupancy normalised spike rate for ',animalname,': Day ',dayt,' epoch ',epocht, ' Tetrode ', tett));

    caxis([0 amax]);

    % Keeps the proportions of the data
    axis image;

    % Labels the sides as 'cm'
    xlabel('cm');
    ylabel('cm');

    % The tick marks are  apart
    set(gca, 'XTick', [0:5:maxx-minx]);
    set(gca, 'YTick', [0:5:maxy-miny]);

    % Adds color legend, tick names
    % Make the colorbar, specifying the location and the tick names

    scaleint = [10];

    datascale=[0:scaleint:amax];

    tickname = [0:scaleint:amax];

    c = colorbar('location','EastOutside','YTickLabel', tickname);

    % Adds text 'Legend' to color legend
    
    outliercutoff = 25;
    outliercuttofft = num2str(outliercutoff);
    
    if ot == 1;
        
        ylabel(c,strcat('Hz o:>',outliercuttofft, 'Hz'));

    else
        ylabel(c,strcat('Hz'));
        
    end
        
    % Changes the tick intervals of color legend at intervals defined by scalteint

    set(c, 'YTick',datascale);
    
    if ot ==1;
    
        hold on;
     
        ylabel(c,strcat('Hz o:>',outliercuttofft, 'Hz'));
    
        [outx,outy] = find(imagedata>outliercutoff);

        plot(outy,outx,'o','MarkerEdgeColor',[1 1 1]);
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
    figurename = strcat(animalname,'_spr_d',dayt,'e',epocht,'t',tett);
    saveas(gcf, figurename, 'pdf');

    % Closes the figure
    close;
end    

if sspr == 1

    % Smoothed firing rate

    figure;

    nc = 1024;
    cmap = jet(nc);
    cmap(1,:) = 0.8;
    colormap(cmap);

    imagedatas = flipud(output.smoothedspikerate);

    bmax = max(output.smoothedspikerate(:));
    bmin = min(output.smoothedspikerate(:));
    tett = num2str(tet);

    % Firing rate
    image(imagedatas,'CDataMapping', 'scaled');

    % Generates a title from information provided
    title(strcat('Smoothed occupancy normalised spike rate for ',animalname,': Day ',dayt,' epoch ',epocht, ' Tetrode ', tett));

    caxis([0 bmax]);

    % Keeps the proportions of the data
    axis image;

    % Labels the sides as 'cm'
    xlabel('cm');
    ylabel('cm');

    % The tick marks are  apart
    set(gca, 'XTick', [0:5:maxx-minx]);
    set(gca, 'YTick', [0:5:maxy-miny]);

    % Adds color legend, tick names
    % Make the colorbar, specifying the location and the tick names

    scaleint = [10];

    datascale=[0:scaleint:bmax];

    tickname = [0:scaleint:bmax];

    c = colorbar('location','EastOutside','YTickLabel', tickname);
    
    % Adds text 'Legend' to color legend
    outliercutoff = 25;
    outliercuttofft = num2str(outliercutoff);
    
    if ot == 1;
        
        ylabel(c,strcat('Hz o:>',outliercuttofft, 'Hz'));

    else
        ylabel(c,strcat('Hz'));
        
    end
        
    % Changes the tick intervals of color legend at intervals defined by scalteint

    set(c, 'YTick',datascale);
    
    if ot ==1;
    
        hold on;

        ylabel(c,strcat('Hz o:>',outliercuttofft, 'Hz'));
    
        [outx,outy] = find(imagedatas>outliercutoff);

        plot(outy,outx,'o','MarkerEdgeColor',[1 1 1]);
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
    figurename = strcat(animalname,'_sspr_d',dayt,'e',epocht,'t',tett);
    saveas(gcf, figurename, 'pdf');
    
    % Closes the figure
    close;
end

if spn == 1

    % Spikes

    figure;

    nc = 1024;
    cmap = jet(nc);
    cmap(1,:) = 0.8;
    colormap(cmap);

    imagedatasp = flipud(output.spikes);

    cmax = max(output.spikes(:));
    cmin = min(output.spikes(:));
    tett = num2str(tet);

    % Firing rate
    image(imagedatasp,'CDataMapping', 'scaled');

    % Generates a title from information provided
    title(strcat('Spikes ',animalname,': Day ',dayt,' epoch ',epocht, ' Tetrode ', tett));

    caxis([0 cmax]);

    % Keeps the proportions of the data
    axis image;

    % Labels the sides as 'cm'
    xlabel('cm');
    ylabel('cm');

    % The tick marks are  apart
    set(gca, 'XTick', [0:5:maxx-minx]);
    set(gca, 'YTick', [0:5:maxy-miny]);

    % Adds color legend, tick names
    % Make the colorbar, specifying the location and the tick names

    scaleint = [100];

    datascale=[0:scaleint:cmax];

    tickname = [0:scaleint:cmax];

    c = colorbar('location','EastOutside','YTickLabel', tickname);

    % Adds text 'Legend' to color legend
    
    outliercutoff = 25;
    outliercuttofft = num2str(outliercutoff);
    
    if ot == 1;
        
        ylabel(c,strcat('Events o:>',outliercuttofft));

    else
        ylabel(c,strcat('Events'));
        
    end
        
    % Changes the tick intervals of color legend at intervals defined by scalteint

    set(c, 'YTick',datascale);
    
    if ot ==1;
    
        hold on;    

        ylabel(c,strcat('Events o:>',outliercuttofft));
    
        [outx,outy] = find(imagedatasp>outliercutoff);

        plot(outy,outx,'o','MarkerEdgeColor',[1 1 1]);
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
    figurename = strcat(animalname,'_spn_d',dayt,'e',epocht,'t',tett);
    saveas(gcf, figurename, 'pdf');
    
    % Closes the figure
    close;
end

cd(currentdir);
    
warning('ON','MATLAB:divideByZero');
