function plottraj(animaldir,day,epoch)
% Plots the trajectory of the animal during the specified epoch
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

% Rescale data to zero

xmin = min(xdata);
xmax = max(xdata);

ymin = min(ydata);
ymax = max(ydata);

posdata = [xdata-xmin,ydata-ymin];


% ----Plot figure----
% Creates a new figure

figure;

% Plots data in black line
p = plot(posdata(:,1),posdata(:,2),'k-');

% Generates a title from information provided
title(strcat('Trajectory for ',animalname,': Day ',dayt,' epoch ',epocht));

% Keeps the size ratio according to the data  
axis image;

% Labels the sides as 'cm'
xlabel('cm');
ylabel('cm');

% The tick marks are 2 cm apart
set(gca, 'XTick', [0:5:max(posdata(:,1))]);
set(gca, 'YTick', [0:5:max(posdata(:,2))]);

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
figurename = strcat(animalname,'_traj_d',dayt,'e',epocht);
saveas(gcf, figurename, 'pdf');

% Closes the figure
close;

% ----Change directory----
% Returns to directory at the start
cd(currentdir);














