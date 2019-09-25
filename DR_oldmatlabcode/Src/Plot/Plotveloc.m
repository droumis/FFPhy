function plotveloc(animaldir,day,epoch)
% Plots the velocity of the animal during the specified epoch
% Saves the figure as a pdf in Plot directory
% animal dir - name of the directory containing processed data
% day - experimental day to be analysed
% epoch - experimental epoch to be analysed

% get the current directory to go back to after the function is done
currentdir = pwd

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
filename = strcat(datadir,animaldir,'/',animalname,'pos',dsz,dayt,'.mat');

% Load the mat file
eval(['load ', filename]);

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

% ----Plot figure----
% Defines the colormap, with the lowest value set to grey

nc = 1024;
cmap = jet(nc);
cmap(1,:) = 0.8;
colormap(cmap);

% max color scale value
amax = 25

% Creates an image of the data

image(rot90(accumdata,1),'CDataMapping', 'scaled');

% Generates a title from information provided
title(strcat('Mean velocity for ',animalname,': Day ',dayt,' epoch ',epocht));

caxis([0 amax]);

% Keeps the proportions of the data
axis image;

% Labels the sides as 'cm'
xlabel('cm');
ylabel('cm');

% The tick marks are 2 cm apart
set(gca, 'XTick', [0:5:xmax-xmin]);
set(gca, 'YTick', [0:5:ymax-ymin]);

% Adds color legend, tick names
% Make the colorbar, specifying the location and the tick names

scaleint = [10];

datascale=[0:2:amax];

tickname = [0:2:amax];

c = colorbar('location','EastOutside','YTickLabel', tickname);

% Adds text 'Legend' to color legend

ylabel(c,'cm/s o:outlier');

% Changes the tick intervals of color legend at intervals defined by scalteint

set(c, 'YTick',datascale);

hold on;

[outx,outy] = find(rot90(accumdata,1)>25);

plot(outy,outx,'o','MarkerEdgeColor',[1 1 1]);

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
figurename = strcat(animalname,'_velocity_d',dayt,'e',epocht);
saveas(gcf, figurename, 'pdf');

% Closes the figure
close;

% ----Change directory----
% Returns to directory at the start
cd(currentdir);
 













