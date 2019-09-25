
load ../../Ten/tenspikes05
load ../../Ten/tenpos05

% data = video_reader('ten05.mpeg','ten05/ten05.mpegoffset','ten05/ten05.cpupostimestamp');

day = 5;
ep = 4;

ss = spikes{day}{ep};

clear neuron;

k = 0;
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
    k = k + 1;
    neuron(k).tetrode = i;
    neuron(k).index = j;
    neuron(k).data = ss{i}{j}.data;
    neuron(k).spikewidth = ss{i}{j}.spikewidth;
    neuron(k).nspikes = size(neuron(k).data,1);
  end
end

NN = k;

colors = jet(NN);

t1 = double(str2ts('1:14:48'))/10000; 
t2 = double(str2ts('1:15:00'))/10000; 

MAXSPIKES = 1000;
SCALE = 4;


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
% t1 = double(str2ts('1:47:35'))/10000; t2 = double(str2ts('1:47:48'))/10000; Ks = [11 23 30 32 34 47 52];
% t1 = double(str2ts('1:14:46'))/10000; t2 = double(str2ts('1:15:00'))/10000; 
Ks = [8 22 34 39 53];

colors = jet(length(Ks));

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
  picture = flipdim(data.get_picture(frame_idx),2);
  picture(:,75:175,:) = 0;
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
  height = 10;
  spacing = 15;
  for k = 1:length(Ks)
    baseline = 100 + k*spacing;
    th = text(SCALE*80,SCALE * (baseline + height*0.5) ,sprintf('neuron %d',k));
    set(th,'color',colors(k,:),'fontsize',18,...
      'fontname', 'trebuchet ms', 'fontweight', 'demi', ...
      'verticalalignment','middle');
  end

  line(SCALE*[170 170],SCALE*[100 + spacing - height, 100 + (length(Ks))*15 + 2*height], ...
    'linewidth', 0.5, ...
    'color', [0.75 0.75 0.75]);

  winsize = 6;
  k = 0;
  for n = Ks
    ww = neuron(n).data(:,1) > max(tf2-winsize,t1) & neuron(n).data(:,1) < tf2;
    k = k + 1;
    if sum(ww) > 0
      tt = neuron(n).data(ww,1);
      tt = (tt - (tf2 - winsize))/winsize;
      tt = tt * 75 + 20 + 75;
      baseline = 100 + k*spacing;
      spikeTrain(SCALE* tt, SCALE* baseline, SCALE* height,'color',colors(k,:),'linewidth',2);
    end
  end

  axis(SCALE*[75  295  80  (80 + 0.75*0.75*(295-75))]);  % assume title is 25%

  % (160x120 = 4x3 -> 214x120 = 16x9)
  ff2 = getframe(gca);
  imwrite(ff2.cdata(1:end-1,1:end-1,:),sprintf('video/fr_%03d.png',ii),'png');

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

tetrode = 9;
load 09-118/ten05-09.mat
waves = double(waves);
timestamps = double(timestamps);
idx = find(timestamps > t1 & timestamps < t2);
load(sprintf('../../Ten/EEG/teneeg%02d-%d-%02d.mat',day,ep,tetrode));
ee = eeg{day}{ep}{tetrode};
estart = floor((t1 - ee.starttime) * ee.samprate - ee.samprate * winsize/10000 / 2);
elength = round((t2 - t1 + winsize/10000) * ee.samprate);

background = ee.data(estart + [0:elength-1]);
background = interp(background,20);

tstart = (t1*10000 - winsize)*3;
for i = 1:length(idx)
  w = waves(:,1,idx(i));
  widx = timestamps(idx(i))*3 - tstart - 8;
  background(widx + [0:39]) = background(widx + [0:39]) + w;
end

for k = 1:length(Ks)
  tetdir = dir(sprintf('%02d-*',neuron(Ks(k)).tetrode));
  tetdir = tetdir([tetdir.isdir]==1).name;
  load(fullfile(tetdir,sprintf('ten%02d-%02d.mat',day,neuron(Ks(k)).tetrode)));
  waves = double(waves);
  timestamps = double(timestamps);
  ww = find(neuron(Ks(k)).data(:,1) > t1 & neuron(Ks(k)).data(:,1) < t2);
  idx = lookup(neuron(Ks(k)).data(:,1)*10000, timestamps);
  mw = mean(waves(:,:,idx),3);
  [y,bestchannel] = max(mean(abs(mw),1));
  for i = 1:length(ww)
    widx = round(neuron(Ks(k)).data(ww(i)) * 30000 - tstart) - 8;
    background(widx + [0:39]) = background(widx + [0:39]) + ...
      waves(:,bestchannel,idx(ww(i)))*(100/y) *20;
  end
end

bb = resample(background,147,100); % resample from 30 kHz to 44.1 kHz
wavwrite(bb / max(abs(1.2*bb)),44100,'video/test.wav');

