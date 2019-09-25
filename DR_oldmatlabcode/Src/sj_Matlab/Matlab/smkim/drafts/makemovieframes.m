% makemovieframes.m
%
%   by smk, last updated 5 February 2008
%
% This script is designed to generate a sequence of movie frames which show
% the locations at which neurons spike, overlaid on top of video frames of
% behavior in an environment.
%
% This script will NOT work straight-out-of-the-box with your data. Read the
% code and modify to suit your needs.
%
%
% smk's Quicktime Pro key: V78F-68N2-DJUT-GDHD-BXLT
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs to feed into the script: you will need to modify these

% sprintf-style string pattern for loading the individual video frames
% note that the video frames must be in an image format that can be read
% by matlab's imread function
video_frames_filename_pattern = '/home/smkim/foo/foo%05d.png';

% sprintf-style string pattern for exporting png output
output_filename_pattern = '/home/smkim/foo/out%05d.png';

% a cell array of spike times: spiketimes{i} is a vector of spiketimes for
% the ith neuron
load('spiketimes_day1_run1.mat'); 
temp = spiketimes; 
clear('spiketimes');
spiketimes{1} = temp{10}.data;
spiketimes{2} = temp{9}.data;
spiketimes{3} = temp{6}.data;
spiketimes{4} = temp{8}.data;
clear('temp');

% define the colors for plotting spikes of different neurons: 
% these must be 3-element RGB vectors, not matlab predefined character codes!
colors{1} = [1 1 0]; 
colors{2} = [0 1 0];
colors{3} = [1 0 1];
colors{4} = [0 1 1];
if numel(colors) ~= numel(spiketimes)
    error('cells in spiketimes must match cells in colors');
end

% load the raw position data
% an array of position data: posdata(:,1) contains timestamps, 
% posdata(:,2) contains x coordinates, posdata(:,3) contains y coordinates
load('rawpos.mat');
temp = rawpos{1}{2}.data;
clear('rawpos');
posdata = temp;

% define the time interval of the movie. we need to align the posdata, the 
% video frames, and the spike times.
% range of desired rows of posdata
posdata_indices = 2729:3417;
% corresponding start time, in seconds
start_time = posdata(posdata_indices(1),1);
% the number of the first video frame (to be combined with 
% video_frames_filename_pattern)
video_frame_num = 1;

% define miscellaneous properties of the figure
% read in the first frame to get width and height (remember that the image
% function draws images in orientation of dimension_1 rows x dimension_2 cols
video_frame = imread(sprintf(video_frames_filename_pattern,video_frame_num));
video_frame_height = size(video_frame,1);
video_frame_width = size(video_frame,2);
% length of time interval over which spike trains are drawn
spike_train_window = 3; 
% pixel resolution of the spike train rasters
pixels_per_sec = 150;
% number of pixels to space apart different elements of the figure; you can 
% tinker with this value to see its effect on the figure layout
horizontal_space = 50;
vertical_space = 40;
% font size; all fonts are scaled relative to this
font_size = 16;
% figure caption in the upper left
figure_caption = 'CA1';

% transform posdata coordinates to match the orientation of the video frame
posdata(:,2) = video_frame_width - temp(:,3);
posdata(:,3) = video_frame_height - temp(:,2);
clear('temp');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw most of the first frame
h_figure = figure;
set(h_figure,'Renderer','opengl','Color','k','Position',[0 0 900 400]);
h_axes = axes('XDir','normal','YDir','normal',...
    'XLimMode','manual','YLimMode','manual','ZLimMode','manual',...
    'XLim',[(-2*horizontal_space-pixels_per_sec*spike_train_window) ...
    video_frame_width],'YLim',[0 video_frame_height],...
    'DataAspectRatio',[1 1 1],'Color','k','Units','pixels', ...
    'Parent',h_figure);
hold on; % lock the axes properties

% draw the video frame
h_image = image([1 video_frame_width],[1 video_frame_height], ...
    video_frame,'Parent',h_axes);

% draw the figure caption
text(-horizontal_space-pixels_per_sec*spike_train_window, ...
    video_frame_height,figure_caption,'HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',1.25*font_size,'Color','w');

% draw a vertical line corresponding to the current time in the spike trains
line([-horizontal_space -horizontal_space], ...
    [0 video_frame_height],'Color','w','LineWidth',2);
% and draw a matching text object which reports the current clock time
h_clock = text(-horizontal_space,video_frame_height, ...
    sprintf('%04.2f ',start_time),'VerticalAlignment','top', ...
    'HorizontalAlignment','right','FontSize',font_size,'Color','w');

% draw a time scale bar for the spike trains
line([(-horizontal_space-pixels_per_sec*spike_train_window) ...
    (-horizontal_space-(spike_train_window-0.5)*pixels_per_sec)],[0 0], ...
    'LineWidth',4,'Color','w');
text(-horizontal_space-pixels_per_sec*spike_train_window,10,'500 ms', ...
    'HorizontalAlignment','left', ...
    'FontSize',font_size,'Color','w');

% draw text labels for each neuron's spike train
for i = 1:length(spiketimes)
    text(-horizontal_space-pixels_per_sec*spike_train_window, ...
        video_frame_height-vertical_space*(i+1), ...
        ['neuron ' num2str(i)],'Color',colors{i}, ...
        'HorizontalAlignment','right','VerticalAlignment','middle', ...
        'FontSize',font_size);
end

% minimum size of spike bubbles
minimum_bubble_size = 5;
% maximum size of spike bubbles (i.e., size when first spawned)
maximum_bubble_size = 40;
% decay rate of bubbles, as a fraction per frame
bubble_decay_rate = 0.9;
% cell array containing vectors of handles to the bubbles representing 
% locations at which spikes occurred
for unit = 1:length(spiketimes)
    % h_spike_locations{i} contains a vector of handles to the marker bubbles
    % for the ith neuron
    h_spike_locations{unit} = [];
end

% how tall to draw spike tickmarks?
tickmark_halfheight = 8;
% cell array containing vectors of handles to the tickmarks representing 
% spike times in the time series raster plot on the left
for unit = 1:length(spiketimes)
    % h_spike_tickmarks{i} contains a vector of handles to the spike tickmarks
    % for the ith neuron
    h_spike_tickmarks{unit} = [];

    % draw faded tickmarks to represent impending spikes in the near future
    temp = find( (spiketimes{unit} > posdata(1,1)) ...
        & (spiketimes{unit} <= posdata(2,1)) );
    for j = temp(:)'
        xcoord = floor( -horizontal_space ...
            - pixels_per_sec * (posdata(1,1) - spiketimes{unit}(j)) );
        ycoord = video_frame_height - vertical_space*(unit+1);
        h_spike_tickmarks{unit}(end+1) = line( ...
            'XData',[xcoord xcoord],'YData',ycoord + ...
            [-tickmark_halfheight tickmark_halfheight], ...
            'Color',0.2*colors{unit},'LineWidth',2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now draw the remaining frames in sequence, with spike overlays
for posidx = posdata_indices

    % load and display video frame
    video_frame = imread( sprintf(video_frames_filename_pattern, ...
        video_frame_num) );
    set(h_image,'CData',video_frame);
    
    % update the clock time
    set(h_clock,'String',sprintf('%04.2f',posdata(posidx,1)));
    
    % iterate through neurons and update stuff
    for unit = 1:length(spiketimes)

        % shrink previously existing spike marker bubbles
        for j = 1:numel(h_spike_locations{unit})
            if get(h_spike_locations{unit}(j),'MarkerSize') > minimum_bubble_size
                set(h_spike_locations{unit}(j),'MarkerSize', ...
                    bubble_decay_rate*get(h_spike_locations{unit}(j),'MarkerSize'));
            end
        end
        % spawn new spike marker bubbles
        temp = find( (spiketimes{unit} > posdata(posidx,1)) ...
            & (spiketimes{unit} <= posdata(posidx+1,1)) );
        for j = temp(:)'
            interp_factor = (spiketimes{unit}(j) - posdata(posidx,1)) ...
                / (posdata(posidx+1,1) - posdata(posidx,1));
            xcoord = interp_factor*posdata(posidx,2) + (1-interp_factor)*posdata(posidx+1,2);
            ycoord = interp_factor*posdata(posidx,3) + (1-interp_factor)*posdata(posidx+1,3);
            h_spike_locations{unit}(end+1) = line('XData',[xcoord],'YData',[ycoord], ...
                'Color',colors{unit},'Linestyle','none','Marker','o','MarkerSize',maximum_bubble_size);
        end

        % update spike rasters
        % tickmarks march across the window at a rate of pixels_per_sec
        % first, erase all the old tick marks
        for j = 1:numel(h_spike_tickmarks{unit})
            delete(h_spike_tickmarks{unit}(j));
        end
        h_spike_tickmarks{unit} = [];
        % draw faded tickmarks to represent impending spikes in the near future
        temp = find( (spiketimes{unit} > posdata(posidx,1)) ...
            & (spiketimes{unit} <= posdata(posidx+1,1)) );
        for j = temp(:)'
            xcoord = floor( -pixels_per_sec * (posdata(posidx,1) - spiketimes{unit}(j)) ...
                - horizontal_space );
            ycoord = video_frame_height - vertical_space*(unit+1);
            h_spike_tickmarks{unit}(end+1) = line( ...
                'XData',[xcoord xcoord],'YData',ycoord+[-tickmark_halfheight tickmark_halfheight], ...
                'Color',0.2*colors{unit},'LineWidth',2);
        end
        % draw tickmarks to represent spikes which occurred in the recent past
        temp = find( (spiketimes{unit} > posdata(posidx,1) - spike_train_window) ...
            & (spiketimes{unit} <= posdata(posidx,1)) );
        for j = temp(:)'
            xcoord = floor( -pixels_per_sec * (posdata(posidx,1) - spiketimes{unit}(j)) ...
                - horizontal_space );
            ycoord = video_frame_height - vertical_space*(unit+1);
            h_spike_tickmarks{unit}(end+1) = line( ...
                'XData',[xcoord xcoord],'YData',ycoord+[-tickmark_halfheight tickmark_halfheight], ...
                'Color',colors{unit},'LineWidth',2);
        end
    end

    % save the drawn frame
    outfilename = sprintf(output_filename_pattern,video_frame_num);
    imobj = getframe(h_figure);
    %imwrite(imobj.cdata,outfilename,'png');
    pause(0.01);

    % very important!: advance to next video frame
    video_frame_num = video_frame_num + 1;
end


