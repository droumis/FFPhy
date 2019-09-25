

dir = '/data25/sjadhav/RippleInterruption/spikevideo';
cd(dir);

load REfspikes07
load REfpos07

data = video_reader('07_020611.mpeg','07_020611.mpegoffset','day_date.postimestamp');

day = 7;
ep = 4; epoch = ep;

ss = spikes{day}{ep};

clear neuron;

k = 0;
cellorder = [3 1 4 2];%  Tet12Ceel2, Tet9Cell1, Tet12Cell4, Tet9Cell4
cellorder = [2 4 1 3];
for i = 1:length(ss)
  if isempty(ss{i})
    continue;
  end
  for j = 1:length(ss{i})
    if isempty(ss{i}{j})
      continue;
    end
    if isempty(ss{i}{j}.data)
      continue;
    end
    
    % Shantanu - pick my cells
    if ((i==9) && (j==1 || j==4)) || ((i==12) && (j==2 || j==4))
        i;
        j;
        k = k + 1;
        neuron(cellorder(k)).tetrode = i;
        neuron(cellorder(k)).index = j;
        neuron(cellorder(k)).data = ss{i}{j}.data;
        neuron(cellorder(k)).spikewidth = ss{i}{j}.spikewidth;
        neuron(cellorder(k)).nspikes = size(neuron(k).data,1);
    end
  end
end

eegtets = [5 9 10 12];
% Load eeg files
cnt=0;
for tet=eegtets
    cnt=cnt+1;
    eegfile = sprintf('%seeg%02d-%01d-%02d.mat', 'REf', day,epoch,tet);
    load(eegfile);
    eval(['eegdata(',num2str(cnt),').data = eeg{day}{epoch}{tet}.data;'])
    if cnt==1
        teeg = geteegtimes(eeg{day}{epoch}{tet});
        eegstart = eeg{day}{epoch}{tet}.starttime;
        eegsamprate = round(eeg{day}{epoch}{tet}.samprate);
    end
end


NN = k;

colors = jet(NN);

% t1 = double(str2ts('0:36:48'))/10000; 
% t2 = double(str2ts('1:38:00'))/10000; 

t1=4244.5; t2=4252;
%t1=4698; t2=4707;

%t1=4334; t2=4342; %t1=4335; t2=4342;
%t1=4455; t2=4462;
%t1=4511; t2=4525; %t1=4520; t2=4525;
%t1=4571; t2=4580;

%t1=4832; t2=4840;

MAXSPIKES = 1000;
SCALE = 4;

% tetrode 12
k = 0;
Ks = [];
for n = 1:NN
  if (neuron(n).nspikes > MAXSPIKES)
    continue;
  end
  ww = neuron(n).data(:,1) > t1 & neuron(n).data(:,1) < t2;
  if sum(ww) > 0
    k = k + 1;
    Ks = [Ks n];
  end
end




% Ks = [8 22 34 39 53];

colors = jet(length(Ks));
colors(1,:) = [1 0.5 0];
%colors = {'k','m','g','c','r'};

cmpp = pos{day}{ep}.cmperpixel;

vid_start = find(double(data.timestamps) < t1*10000, 1, 'last');
vid_end = find(double(data.timestamps) > t2*10000, 1, 'first');

winsize = mean(diff(double(data.timestamps(vid_start:vid_end))));

picture = data.get_picture(vid_start);
figure(1);



% for frame_idx  = vid_end
ii = 0;
for frame_idx  = vid_start : vid_end 
  ii = ii + 1;
  hold off
  %picture = flipdim(data.get_picture(frame_idx),2);
  picture = data.get_picture(frame_idx);
  
  %picture(:,75:175,:) = 0;
  picture(:,1:105,:) = 0; % Shantanu W track
  picture = imresize(picture,SCALE);
  h = imshow(picture);
  hold on

  tf1 = (double(data.timestamps(frame_idx)) - winsize/2)/10000;
  tf2 = (double(data.timestamps(frame_idx)) + winsize/2)/10000;
  k = 0;
  for n = Ks
    ww = neuron(n).data(:,1) > t1 & neuron(n).data(:,1) < tf2;
    k = k + 1;
    if sum(ww) > 0
      plot(SCALE*neuron(n).data(ww,2)/cmpp,SCALE*neuron(n).data(ww,3)/cmpp,'o', ...
        'markersize', 12, ...
        'markerfacecolor',colors(k,:),'markeredgecolor','none');
    end
  end

  % blank space for spikes: xspace 75 - 175, yspace = 80 - 200
  % Shantanu: xspace 225-325  yspace 80-300
  height = 10; height2 = 15; spacing2 = 20;
  spacing = 15;
  %height = 15;
  %spacing = 20;
  for k = 1:length(Ks)
    baseline = 100 + k*spacing; 
    th = text(SCALE*5,SCALE * (baseline + height*0.5) ,sprintf('Neuron %d',k));
    set(th,'color',colors(k,:),'fontsize',30,...
      'fontname', 'trebuchet ms', 'fontweight', 'demi', ...
      'verticalalignment','middle');
  
  
    baseline = 5 + k*spacing2; 
    th = text(SCALE*5,SCALE * (baseline + height2*0.5) ,sprintf('EEG %d',k));
    set(th,'color',colors(k,:),'fontsize',30,...
      'fontname', 'trebuchet ms', 'fontweight', 'demi', ...
      'verticalalignment','middle');
%     th = text(SCALE*80,SCALE * (baseline + height*0.5) ,sprintf('neuron %d',k));
%     set(th,'color',colors(k,:),'fontsize',18,...
%       'fontname', 'trebuchet ms', 'fontweight', 'demi', ...
%       'verticalalignment','middle');
  end

%   line(SCALE*[170 170],SCALE*[100 + spacing - height, 100 + (length(Ks))*15 + 2*height], ...
%     'linewidth', 0.5, ...
%     'color', [0.75 0.75 0.75]);

  line(SCALE*[104 104],SCALE*[100 + spacing - height, 100 + (length(Ks))*15 + 2*height], ...
    'linewidth', 0.5, 'color', [0.75 0.75 0.75]);

  line(SCALE*[104 104],SCALE*[5 + spacing - height, 100 + (length(Ks))*15 + 2*height], ...
    'linewidth', 0.5, 'color', [0.75 0.75 0.75]);

  winsize = 3.6;
  k = 0;
  for n = Ks
    ww = neuron(n).data(:,1) > max(tf2-winsize,t1) & neuron(n).data(:,1) < tf2;
    k = k + 1;
    if sum(ww) > 0
      tt = neuron(n).data(ww,1);
      tt = (tt - (tf2 - winsize))/winsize;
      %tt = tt * 75 + 20 + 75;
      tt = tt * 75 + 25;
      baseline = 100 + k*spacing;
      spikeTrain(SCALE* tt, SCALE* baseline, SCALE* height,'color',colors(k,:),'linewidth',2);
    end
  end
  
  % EEG data
  for k = 1:length(eegtets)
      eidx = teeg > max(tf2-winsize,t1) & teeg < tf2;
      eet = teeg(eidx);
      eed = eegdata(k).data(eidx);
      
      eet = (eet - (tf2 - winsize))/winsize;
      eet = eet .* 75 + 25;
      dscale = max(eed)-min(eed); pscale = height2*9/10;
      eed = eed.*(pscale/dscale);     
    
      baseline = 12 + k*spacing2; 
      plot(SCALE* eet, SCALE*(baseline+eed),'color',colors(k,:),'linewidth',1);
     
  end

  %axis(SCALE*[75  250 80  (80 + 0.75*0.75*(295-75))]);  % assume title is 25%
  % (160x120 = 4x3 -> 214x120 = 16x9)

  %axis(SCALE*[0  320 0  240]);  
  % 320X240
  
  ff2 = getframe(gca);
  imwrite(ff2.cdata(1:end-1,1:end-1,:),sprintf('video2/spkeeg_%03d.png',ii),'png');

end

figure(2)
clf
subplot(6,1,1)
ww = pos{day}{ep}.data(:,1) > t1 & pos{day}{ep}.data(:,1) < t2;
plot(pos{day}{ep}.data(ww,1), pos{day}{ep}.data(ww,5),'r');
subplot(6,1,[2:5])
for i = 1:length(Ks)
  ww = neuron(Ks(i)).data(:,1) > t1 & neuron(Ks(i)).data(:,1) < t2;
  plot(neuron(Ks(i)).data(ww,1),i,'s', ...
      'markerfacecolor',colors(i,:),'markeredgecolor','none');
  hold on
end


%return


% ---------------Done with image files ----------------------------


load 07_020611-12.mat % spike waves
waves = double(waves);
timestamps = double(timestamps);
idx = find(timestamps > t1 & timestamps < t2);

%Eeg tetrode for background noise
%EEg or ripple file
tetrode = 9;
%EEG file
load(sprintf('REfeeg%02d-%d-%02d.mat',day,ep,tetrode));
ee = eeg{day}{ep}{tetrode};
estart = floor((t1 - ee.starttime) * ee.samprate - ee.samprate * winsize/10000 / 2);
elength = round((t2 - t1 + winsize/10000) * ee.samprate);
%EEG noise in background
background = ee.data(estart + [0:elength-1]);
background = interp(background,20);


%If ripple file
load(sprintf('REfripple%02d-%d-%02d.mat',day,ep,tetrode));
% ee = ripple{day}{ep}{tetrode};
% estart = floor((t1 - ee.starttime) * ee.samprate - ee.samprate * winsize/10000 / 2);
% elength = round((t2 - t1 + winsize/10000) * ee.samprate);
% % EEG noise in background
% background = ee.data(estart + [0:elength-1],1);
% background = interp(background,20);



tstart = (t1*10000 - winsize)*3;
%MultiUnit = all spikes, added to background 
for i = 1:length(idx)
  w = waves(:,1,idx(i));
  widx = timestamps(idx(i))*3 - tstart - 8;
  background(widx + [0:39]) = background(widx + [0:39]) + w;
end

for k = 1:length(Ks)
  
%   tetdir = dir(sprintf('%02d-*',neuron(Ks(k)).tetrode));
%   tetdir = tetdir([tetdir.isdir]==1).name;

  load(fullfile(sprintf('%02d_020611-%02d.mat',day,neuron(Ks(k)).tetrode)));
  waves = double(waves);
  timestamps = double(timestamps);
  ww = find(neuron(Ks(k)).data(:,1) > t1 & neuron(Ks(k)).data(:,1) < t2);
  idx = lookup(neuron(Ks(k)).data(:,1)*10000, timestamps);
  mw = mean(waves(:,:,idx),3);
  [y,bestchannel] = max(mean(abs(mw),1));
  % This for loop adds spike sounds to background sound
  for i = 1:length(ww)
    widx = round(neuron(Ks(k)).data(ww(i)) * 30000 - tstart) - 8;
    background(widx + [0:39]) = background(widx + [0:39]) + ...
      waves(:,bestchannel,idx(ww(i)))*(100/y) *20;
  end
end

bb = resample(background,147,100); % resample from 30 kHz to 44.1 kHz
wavwrite(bb / max(abs(1.2*bb)),44100,'video2/REf0704.wav');



% ----------------------------- End ----------------------





% To make video file, use ffmpeg as follows
%(1) either of the two lines to get the video from the png files

% ii = 0;
% for frame_idx  = vid_start : vid_end
%     ii = ii + 1;
%     %ffmpeg -r 29.97 -i fr_%03d.png REf0704spike1.mp4
%     ffmpeg -r 29.97 -i fr_%03d.png -vcodec mpeg2video REf0704spike1.mpeg
% end
% 
% %(2) to get sound
% ffmpeg -i REf0704spike1.mpeg -i REf0704spike1.wav REf0704spike1av.mpeg


% OR video in matlab - Big file
% (1) Video
%aviobj = avifile('REf0704spike1matlab.avi','fps',29.97);
%for i=1:226, f=sprintf('fr_%03d.png',i); A=imread(f); aviobj = addframe(aviobj,A);end
%aviobj=close(aviobj);

%(2) to get sound
% ffmpeg -i testmatlab.avi -i REf0704spike1.wav REf0704spike1matlabav.mpeg
