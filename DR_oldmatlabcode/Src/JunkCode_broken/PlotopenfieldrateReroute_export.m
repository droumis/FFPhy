function [output] = PlotopenfieldrateReroute_export(animaldir,tet,cell,day,epoch, binsize,std,sspr,ot,excludetime)

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
%The output is a structure {1,day}{1,epoch}{1,tet}{1,cell}. The place field
%matrix and the normalised place field matrix are saved.



warning('OFF','MATLAB:divideByZero');

% get the current directory to go back to after the function is done
currentdir = pwd;

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
cellt=num2str(cell);

% Loads the position data
posfilename = strcat(datadir,animaldir,'/',animalname,'_',dsz,dayt,'.mat');
eval(['load ', posfilename]);

posfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'linpos',dsz,dayt,'.mat');
load(posfilename);


% Loads the spike data
spkfilename = strcat(datadir,animaldir,'/',animalname,'spikes',dsz,dayt,'.mat');
eval(['load ', spkfilename]);


% ---- Extract relevant data ----
% Extracts spike time position data from spike file
posx = spikes{1,day}{1,epoch}{1,tet}{1,cell}.data(:,2);
posy = spikes{1,day}{1,epoch}{1,tet}{1,cell}.data(:,3);
posindexfield = spikes{1,day}{1,epoch}{1,tet}{1,cell}.data(:,7);
time = spikes{1,day}{1,epoch}{1,tet}{1,cell}.data(:,1);



if ~isempty(spikes)
    spikes = [time posx posy posindexfield];
    
else
    spikes = [0 0 -1];
end

g = gaussian2(std,(3*std));


% Extracts position data from pos file
posposorig = [Data{day}{epoch}.Pos.correcteddata(:,1) Data{day}{epoch}.Pos.correcteddata(:,2) Data{day}{epoch}.Pos.correcteddata(:,3)];
timestep = posposorig(2,1) - posposorig(1,1);





% get included times
goodspikes = [];
goodspikes = spikes;
pospos = posposorig;
if ~isempty(excludetime)
    exclind = isExcluded(time,excludetime);  %find excluded time indices
    %pospos = posposorig(exclind==1,:);
    % get exclude index
    
    goodspikes = spikes(exclind==0,:);
    badspikes=spikes(exclind==1,:);
    exclind2 = isExcluded(posposorig(:,1),excludetime);  %find excluded time indices
    pospos = posposorig(exclind2==1,:);
    
end


% %Plot time vs spikes
%
% for inds=1:size(excludetime,1)
%     line([excludetime(inds,1) excludetime(inds,2)],[1 1]);
%     hold on;
% end
% hold on;
% plot(goodspikes(:,1),2,'.g');
% hold on;
% plot(badspikes(:,1),2.5,'.r');
%
% ylim([0 3]);

%make a cell array, where each cell contains data for one trajectory.
%inside each cell [binlocation occupancy spikecount firingrate]

tmpposition = (pospos(:,[2 3]));
tmpspikes = (goodspikes(:,[2 3]));

% ---- Calculations ----
% Calculate firing rate and spatial occupancy
minx = -10;
maxx = 400*0.5;
binx = (minx:binsize:maxx);
miny = -10;
maxy = 300*0.5;
biny = (miny:binsize:maxy);

% get total occupancy for all times in order to get indices of non-occupied positions
[occupos occuposx occuposy] = HIST2(pospos(:,2),pospos(:,3), binx, biny);

smoothedoccupancy = [];
smoothedoccupancy = ones(size(occupos));
smoothedoccupancy = filter2(g, occupos);

%smoothedoccupancy = occupos;


%smoothedoccupancy = smoothedoccupancy(1:66,1:81); % for some strange reason the columns change to 82 in some epochs
occuposzeros = find(smoothedoccupancy == 0);

if ~isempty(tmpposition)
    %     minx = floor(min(tmpposition(:,1)));
    %     maxx = ceil(max(tmpposition(:,1)));
    %     binx = (minx:binsize:maxx);
    %     miny = floor(min(tmpposition(:,2)));
    %     maxy = ceil(max(tmpposition(:,2)));
    %     biny = (miny:binsize:maxy);
    
    
    
    %[output.occupancy output.xticks output.yticks] = HIST2(tmpposition(:,1),tmpposition(:,2), binx, biny);
    [output.occupancy output.xticks output.yticks] = HIST2(pospos(:,2),pospos(:,3), binx, biny);
    %output.occupancy=output.occupancy(1:size(smoothedoccupancy,1)-1,1:size(smoothedoccupancy,2));
    
    
    nonzero = find(output.occupancy ~= 0);
    
    [output.spikes BX BY] = hist2(tmpspikes(:,1), tmpspikes(:,2), binx, biny);
    
    output.spikerate = zeros(size(output.spikes));
    
    output.spikerate(nonzero) = output.spikes(nonzero) ./(timestep* output.occupancy(nonzero) );
    
    % Calculate smoothed occupancy
    
    output.smoothedspikerate = filter2(g,(output.spikerate)); % is this the right filter?
    %output.smoothedspikerate = output.spikerate;
    
    output.smoothedspikerate(ind2sub(size(output.smoothedspikerate),occuposzeros)) = -1;
else
    
    output.spikerate=zeros(maxx,maxy);
    output.smoothedspikerate=zeros(maxx,maxy);
    output.spikes=zeros(maxx,maxy);
    
end
% ----Plot figure----


% Plot the firing rate
% Defines the colormap, with the lowest value set to grey


if sspr == 1
    
    % Smoothed firing rate
    
    figure;
    
    nc = 1024;
    cmap = jet(nc);
    cmap(1,:) = 1;
    colormap(cmap);
    
    imagedatas = flipud(output.smoothedspikerate);
    
%     bmax = max(output.smoothedspikerate(:));
%     bmin = min(output.smoothedspikerate(:));
%     tett = num2str(tet);
%     
%     % Firing rate
%     %image(imagedatas,'CDataMapping', 'scaled');
%     imagesc(imagedatas);
%     
%     % plot circle round rewarded wells
%     % get rewarded wells
%     Wells=Data{1,day}{1,epoch}.Wellinfo.rewardedwells;
%     % get reward well positions
%     Wellpos=linpos{1,day}{1,epoch}.wellSegmentInfo.wellCoord;
%     [output.occupancy output.xticks output.yticks]=HIST2(Wellpos(:,1),Wellpos(:,2), binx, biny);
%     
%     % plot circles around wells
%     hold on;
%     output.occupancy = flipud(output.occupancy);
%     [wella, wellb]=find(output.occupancy==1);
%     Wellpos=[wella,wellb];
%     
%     for wind=1:size(Wellpos,1);
%         plot(Wellpos(wind,2),Wellpos(wind,1),'ow','MarkerSize',15);
%     end
%     
%     
%     
%     
%     % Generates a title from information provided
%     title(sprintf('Smoothed occupancy normalised spike rate for %s: Day %s Epoch %s Tetrode %s cell %s \n %s spikes',...
%         animalname, dayt, epocht, tett,cellt,num2str(size(goodspikes,1))));
%     
%     caxis([bmin bmax]);
%     
%     % Keeps the proportions of the data
%     axis image;
%     
%     % Labels the sides as 'cm'
%     xlabel('cm');
%     ylabel('cm');
%     
%     % The tick marks are  apart
%     set(gca, 'XTick', [0:5:maxx-minx]);
%     set(gca, 'YTick', [0:5:maxy-miny]);
%     
%     % Adds color legend, tick names
%     % Make the colorbar, specifying the location and the tick names
%     
%     scaleint = [10];
%     
%     datascale=[0:scaleint:bmax];
%     
%     tickname = [0:scaleint:bmax];
%     
%     c = colorbar('location','EastOutside','YTickLabel', tickname);
%     
%     % Adds text 'Legend' to color legend
%     outliercutoff = 25;
%     outliercuttofft = num2str(outliercutoff);
%     
%     if ot == 1;
%         
%         ylabel(c,strcat('Hz o:>',outliercuttofft, 'Hz'));
%         
%     else
%         ylabel(c,strcat('Hz'));
%         
%     end
%     
%     % Changes the tick intervals of color legend at intervals defined by scalteint
%     
%     set(c, 'YTick',datascale);
%     
%     if ot ==1;
%         
%         hold on;
%         
%         ylabel(c,strcat('Hz o:>',outliercuttofft, 'Hz'));
%         
%         [outx,outy] = find(imagedatas>outliercutoff);
%         
%         plot(outy,outx,'o','MarkerEdgeColor',[1 1 1]);
%     end
    
    % save placefield file
    
    filename=strcat(datadir,animaldir,'/',animalname,'plf',dsz,dayt,'.mat');
    if(exist(filename,'file'));
        load(filename);
    end
    
    imagedatas_norm=imagedatas./max(imagedatas(:));
    imagedatas_norm(imagedatas_norm<0)=-1;
    
    celldata{1,day}{1,epoch}{1,tet}{1,cell}.plf=imagedatas;
    celldata{1,day}{1,epoch}{1,tet}{1,cell}.plf_norm=imagedatas_norm;
    
    save(filename, 'celldata');

    
    close;
end



cd(currentdir);

warning('ON','MATLAB:divideByZero');
