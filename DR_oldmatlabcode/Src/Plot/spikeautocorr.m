function spikeautocorr(animaldir,day,epoch,tet,gauss)

% Plots autocorrelogram of data
% animal dir - name of the directory containing processed data
% day - experimental day to be analysed
% epoch - experimental epoch to be analysed
% tet - tetrode number
% gauss - convolve data with gauss function or not


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

% Loads the spike data
spkfilename = strcat(datadir,animaldir,'/',animalname,'spikes',dsz,dayt,'.mat');
eval(['load ', spkfilename]);


% ---- Extract relevant data ----
% Extracts spike time position data from spike file multiply by 1000 to get
% to 1ms scale
spktime = 1000*spikes{day}{epoch}{tet}{1}.data(:,1);

% Gaussian convolution

sig = 1;
x = [-5:sig:5]';
k = exp(-(x/sig).^2/2)/(sig*sqrt(2*pi));
z = conv(spktime, k);

% autocorrelation

lag = 500;
[cor,lag] = xcorr(z,z,'coeff');

%stem(lag,cor);

plot(cor);
 
disp('l');




