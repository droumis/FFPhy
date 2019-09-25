

%dir = '/data25/sjadhav/HPExpt/spikevideo/HPa_d2e4_7760to7790';
%cd(dir);

load HPaspikes02
load HPapos02

data = video_reader('02_020212.mpeg','02_020212.mpegoffset','day_date.postimestamp');

% The number of timesamps is one less than number of frames. So
% timestamps(end is spurious). Remove it, and either add a new
% timestamps(1) or a new timestamps (end)

% MANUAL EDITS NOT ALLOWED TO "video_reader" class
%data.timestamps(end)=[];
%data.timestamps = [data.timestamps(1)-33 ;data.timestamps]; % 33 ms frame time

% REMEMBER - LAST TIME STAMP IS SPURIOUS. So make a different structure
usets = double(data.timestamps);
usets(end) = [];
usets = [usets(1)-33; usets];

day = 2;
ep = 4;

ss = spikes{day}{ep};

clear neuron;

k = 0;
%cellorder = [3 1 4 2];%  Tet12Ceel2, Tet9Cell1, Tet12Cell4, Tet9Cell4
%cellorder = [2 4 1 3];

% Video 1
 t1=7660; t2=7690;
% Video 2
% t1=7230; t2=7248;
% Video 3
% t1 = 7314; t2 = 7333; 
% Video 4
% t1=7701; t2 = 7720;
% Video 5
%t1 = 7635; t2 = 7652;

% WRONGcellorder = [1 8 4 3 5 6 9 7 2 10]; % Video1

cellorder = [1 9 4 2 5 6 8 7 3 10]; % Video1
% cellorder = [10 2 7 9 6 5 3 4 8 1]; % Video1 - Cell 10 to cell 1
% cellorder = [3 9 1 7 8 6 5 4 2]; % Video2
% cellorder = [3 8 1 4 7 6 5 2]; % Video3
% cellorder = [2 7 1 6 5 4 3]; % Video 4
% cellorder = [3 8 1 6 7 5 4 2]; % Video 5

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
        % Video1
            if ( i==1 && (j==2 || j==5 || j == 6)) || ( i==4 && (j==1 || j==3 || j==4 || j==5 || j==6)) ...
                    || ( i==7 && (j==1)) || ( i==12 && (j==2))
        % Video2
%         if ( i==1 && (j==1 || j==2 || j==3 || j == 6)) || ( i==4 && (j==1 || j==3 || j==4 || j==5)) ...
%                 || ( i==12 && (j==2))
        % Video 3
%         if ( i==1 && (j==1 || j==2 || j==3 || j == 5 || j == 6)) || ( i==4 && (j==3 || j==4)) ...
%                 || ( i==8 && (j==1))
         % Video4
%         if ( i==1 && (j==1 || j==2 || j==3)) || ( i==4 && (j==1 || j==3 || j==4 || j==5))
        % Video 5
%           if ( i==1 && (j==1 || j==2 || j==3 || j == 6)) || ( i==4 && (j==1 || j==3 || j==4)) ...
%                  || ( i==12 && (j==2))
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

NN = k;

colors = jet(NN);

% t1 = double(str2ts('0:36:48'))/10000;
% t2 = double(str2ts('1:38:00'))/10000;



MAXSPIKES = 1000;
SCALE = 4;

% tetrode 12
k = 0;
Ks = [];
for n = 1:NN
    %if (neuron(n).nspikes > MAXSPIKES)
    %  continue;
    %end
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

vid_start = find(usets < t1*10000, 1, 'last');
vid_end = find(usets > t2*10000, 1, 'first');

winsize = mean(diff(usets(vid_start:vid_end)));

picture = data.get_picture(vid_start);
figure(1);



% %for frame_idx  = vid_end
ii = 0;
for frame_idx  = vid_start : vid_end
    ii = ii + 1;
    hold off
    
    picture = data.get_picture(frame_idx); % Get Picture as is
    picture(:,1:90,:) = 0; % Shantanu W track % Black area on left
    %Now crop off sides of picture to focus W-track
    picture = picture(1:180,1:220,:);
    
    % For flipping picture - DONT. Spikes get screwed Up
    %   picture = flipdim(data.get_picture(frame_idx),1); % Flip Up Down
    %   picture = flipdim(picture,2);  %Flip Left Right - Dont
    %   %picture(:,75:175,:) = 0;
    %   picture(:,1:100,:) = 0; % Shantanu W track % Black area on left
    %   % Now crop off sides of picture to focus W-track
    %   picture = picture(70:220,1:240,:);
    
    
    
    picture = imresize(picture,SCALE);
    h = imshow(picture);
    hold on
    
    tf1 = (usets(frame_idx) - winsize/2)/10000;
    tf2 = (usets(frame_idx) + winsize/2)/10000;
    k = 0;
    for n = Ks
        ww = neuron(n).data(:,1) > t1 & neuron(n).data(:,1) < tf2;
        k = k + 1;
        if sum(ww) > 0
            plot(SCALE*neuron(n).data(ww,2)/cmpp,SCALE*neuron(n).data(ww,3)/cmpp,'o', ...
                'markersize', 12, ...
                'markerfacecolor',colors(k,:),'markeredgecolor','none');
            %         % Adjust x-y posn If flipped, etc
            %       plot(SCALE*((neuron(n).data(ww,2)/cmpp)-75),SCALE*((neuron(n).data(ww,3)/cmpp)-125),'o', ...
            %         'markersize', 12, ...
            %         'markerfacecolor',colors(k,:),'markeredgecolor','none');
        end
    end
    
    % blank space for spikes: xspace 75 - 175, yspace = 80 - 200
    % Shantanu: xspace 225-325  yspace 80-300
    
    height = 8;
    spacing = 12;
    %height = 15;
    %spacing = 20;
    
    %shift=100;
    shift=0;
    for k = 1:length(Ks)
        baseline = shift + k*spacing;
        th = text(SCALE*5,SCALE * (baseline + height*0.5) ,sprintf('CELL %d',k));
        
        set(th,'color',colors(k,:),'fontsize',30,...
            'fontname', 'trebuchet ms', 'fontweight', 'bold', ...
            'verticalalignment','middle');
        %     th = text(SCALE*80,SCALE * (baseline + height*0.5) ,sprintf('neuron %d',k));
        %     set(th,'color',colors(k,:),'fontsize',18,...
        %       'fontname', 'trebuchet ms', 'fontweight', 'demi', ...
        %       'verticalalignment','middle');
    end
    
    %   line(SCALE*[170 170],SCALE*[100 + spacing - height, 100 + (length(Ks))*15 + 2*height], ...
    %     'linewidth', 0.5, ...
    %     'color', [0.75 0.75 0.75]);
    
    line(SCALE*[89 89],SCALE*[shift + spacing - height, shift + (length(Ks))*15 + 2*height], ...
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
            %tt = tt * 75 + 25;
            tt = tt * 75 + 10;
            baseline = shift + k*spacing;
            spikeTrain(SCALE* tt, SCALE* baseline, SCALE* height,'color',colors(k,:),'linewidth',2);
        end
    end
    
    %axis(SCALE*[75  250 80  (80 + 0.75*0.75*(295-75))]);  % assume title is 25%
    % (160x120 = 4x3 -> 214x120 = 16x9)
    
    %axis(SCALE*[0  320 0  240]);
    % 320X240
    
    ff2 = getframe(gca);
    imwrite(ff2.cdata(1:end-1,1:end-1,:),sprintf('video1b/fr_%03d.png',ii),'png');
    
    %if ii==50, keyboard; end
    
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


load 02_020212-01.mat % spike waves
waves = double(waves);
timestamps = double(timestamps);
idx = find(timestamps > t1 & timestamps < t2);

%Eeg tetrode for background noise
%EEg or ripple file
tetrode = 1;
%EEG file
load(sprintf('HPaeeg%02d-%d-%02d.mat',day,ep,tetrode));
ee = eeg{day}{ep}{tetrode};
estart = floor((t1 - ee.starttime) * ee.samprate - ee.samprate * winsize/10000 / 2);
elength = round((t2 - t1 + winsize/10000) * ee.samprate);
%EEG noise in background
background = ee.data(estart + [0:elength-1]);
background = interp(background,20);


%If ripple file
% load(sprintf('REfripple%02d-%d-%02d.mat',day,ep,tetrode));
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
    
    load(fullfile(sprintf('%02d_020212-%02d.mat',day,neuron(Ks(k)).tetrode)));
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
wavwrite(bb / max(abs(1.2*bb)),44100,'video1b/HPa0204.wav');



% ----------------------------- End ----------------------





% To make video file, use ffmpeg as follows
%(1) either of the two lines to get the video from the png files

% ii = 0;
% for frame_idx  = vid_start : vid_end
%     ii = ii + 1;
%     ffmpeg -r 29.97 -i fr_%03d.png REf0704spike1.mp4
%     ffmpeg -r 29.97 -i fr_%03d.png -vcodec mpeg2video REf0704spike1.mpeg
% end
%
% %(2) to get sound
% ffmpeg -i REf0704spike1.mpeg -i REf0704spike1.wav REf0704spike1av.mpeg


% OR video in matlab - Big file
% (1) Video

%   aviobj = avifile('Hpa0204spikematlab1.avi','fps',29.97);
%   for i=1:901, f=sprintf('fr_%03d.png',i); A=imread(f); aviobj = addframe(aviobj,A);end
%   aviobj=close(aviobj);

%(2) to get sound
% ffmpeg -i Hpa0204spikematlab2.avi -i HPa0204.wav HPa0204spike1matlabav.avi
% For UBUNTU machine, might need to specify audio codec
% ffmpeg -i Hpa0204spikematlab2.avi -i HPa0204.wav -acodec mp2 HPa0204spike1matlabav.avi
