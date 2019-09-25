function amplitude_cutoff = check_tt(filename,session)
%CHECK_TT View amplitude distribution of events in .tt file and interactively set amplitude cutoff
%
%   AMPLITUDE_CUTOFF = CHECK_TT(TT_FILENAME,SESSION) reads the
%   threshold-triggered snippets of waveforms samples from the packed binary
%   data format in TT_FILENAME and plots the max amplitude of the waveform
%   snippet as recorded on each channel of the tetrode. The user is given a
%   keyboard prompt to determine a desired amplitude cutoff for censoring
%   high-amplitude noise events. This function is intended to be used as a
%   complement to READ_SPIKE.
%
%   TT_FILENAME must specify an NSpike *.tt file of the type outputted by
%   nspike_extract -spike and must include the path to the file.
%
%   SESSION must be a struct array of session info all matching the same
%   subject/day (because a single .tt file contains data from one subject/day).
%
%   The return value AMPLITUDE_CUTOFF is a row vector of positive real values.
%
%   See also READ_SPIKE, written by smk.
%
%Depends on:
%   IS_SESSION (written by smk)
%   READ_BINARY_RECORDS (written by smk)
%   HISTCN (written by smk)
%
%Written by smk, 2009 June 21.
%

if (exist('is_session') ~= 2)
  error('CHECK_TT depends on m-file IS_SESSION (written by smk)');
end
if (exist('read_binary_records') ~= 2)
  error('CHECK_TT depends on m-file READ_BINARY_RECORDS (written by smk)');
end
if (exist('histcn') ~= 2)
  error('CHECK_TT depends on m-file HISTCN (written by smk)');
end

% Number of channels per tetrode
NCHANNELS = 4;
% Number of samples per window in the original tt file
NPRE_THRESH_POINTS = 8;
NPOST_THRESH_POINTS = 32;
NPOINTS_PER_SPIKE = NPRE_THRESH_POINTS + NPOST_THRESH_POINTS;

% Check that SESSION is a valid session struct array
if ~is_session(session)
  error('SESSION does not appear to be a valid sessions data struct array');
end
% Check whether all elements of session have matching subject and day
if (numel(unique({session(:).subject})) > 1)
  error('SESSION argument has elements whose subject fields differ');
end
if (numel(unique([session(:).day])) > 1)
  error('SESSION argument has elements whose day fields differ');
end

% Verify that filename specifies a real file and includes 
if ~ischar(filename)
  error('filename must be a string');
elseif (exist(filename) ~= 2)
  error('filename %s does not refer to a valid file on search path',filename);
elseif ~isdir(fileparts(filename))
  error('filename %s does not include path',filename)
end

% read data from tt file using READ_BINARY_RECORDS
recformat = cell2struct({ ...
    'timestamp', 'uint32', 1; ...
    'samples'  , 'int16' , NCHANNELS*NPOINTS_PER_SPIKE }, ...
    {'name','type','count'},2);
eoh = '%%ENDHEADER\n';
try
  % (note that this isn't a fully-populated spike data struct array; it will fail
  % to validate with IS_SPIKE)
  spike = read_binary_records(filename,eoh,recformat,Inf);
  num_events = numel(spike.timestamp);
  % samples field has num_events rows and (NCHANNELS*NPOINTS_PER_SPIKE) columns;
  % each row contains samples in interlaced order, e.g. 
  %
  % (channel1,sample1), (channel2,sample1), ... (channelNCHANNELS,sample1), 
  % (channel1,sample2), (channel2,sample2), ... (channelNCHANNELS,sample2), 
  % (channel1,sample3), ...
  % (channel1,sampleNPOINTS_PER_SPIKE), (channel2,sampleNPOINTS_PER_SPIKE), ...
  %
  % We want to convert this to a 3-dimensional array of size 
  % NPOINTS_PER_SPIKE * NCHANNELS * num_events. Also, invert sign to undo the -1
  % gain in NSpike (so that, as expected, positive values mean positive recorded
  % extracellular voltage)
  spike.samples = (-1) * permute(reshape(spike.samples', ...
      [NCHANNELS NPOINTS_PER_SPIKE num_events]),[2 1 3]);
catch
  error('Could not read waveforms data from file %s',filename);
end

% Check whether there are records
if (num_events == 0)
  warning('file %s contains no trigger events',filename);
  amplitude_cutoff = [0 0 0 0];
  return;
end

% Initialize amplitude_cutoff to be Inf
amplitude_cutoff = Inf*ones([1 NCHANNELS 1]);

% Initialize multi-axes figure
fig = figure();

% Plot max-amplitude scatterplots for interactive setting of amplitude_cutoff
MAX_SAMPLE = single(max(abs(spike.samples(:))));
maxval = single(permute(max(abs(spike.samples),[],1),[3 2 1]));
for c = 1:NCHANNELS
  if ~isfinite(amplitude_cutoff(c))
    edges{c} = [linspace(0,MAX_SAMPLE,101)'; Inf];
  else
    edges{c} = [linspace(0,amplitude_cutoff(c),101)'; Inf];
  end
end
colormap(gray);
[row, col] = ndgrid(1:(NCHANNELS-1),2:NCHANNELS);
pair = find(triu(ones(NCHANNELS-1)));
row = row(pair);
col = col(pair);
for i = 1:numel(pair)
  if pair(i)
    h(i) = subplot(NCHANNELS,NCHANNELS-1,pair(i));
    pixels = sparse(histcn([maxval(:,row(i)), ...
        maxval(:,col(i))], ...
        edges{row(i)},edges{col(i)}) > 0);
    % note that pixels is transposed when passed to imagesc because of the way
    % that imagesc orients the data for display
    imagesc([0 edges{row(i)}(end-1)],[0 edges{col(i)}(end-1)],pixels', ...
        'Parent',h(i));
    set(h(i),'XLim',[0 edges{row(i)}(end-1)], ...
        'YLim',[0 edges{col(i)}(end-1)], ...
        'XTick',[0 edges{row(i)}(end-1)], ...
        'YTick',[0 edges{col(i)}(end-1)], ...
        'DataAspectRatio',[1 1 1],'XDir','normal','YDir','normal');
    xlabel(h(i),sprintf('channel %d',row(i)));
    ylabel(h(i),sprintf('channel %d',col(i)));
  end
end

% Plot histogram of event density over time, with recording epoch intervals
% shaded
t_begin = min(session(1).timerange(1),spike.timestamp(1));
t_end = max(session(end).timerange(2),spike.timestamp(end));
bins = linspace(single(t_begin),single(t_end),200);
counts = hist(single(spike.timestamp),bins,'Color','k');
h(numel(pair)+1) = subplot(NCHANNELS,NCHANNELS-1, ...
    (NCHANNELS-1)^2 + [1:NCHANNELS-1]);
set(gca,'XLim',[t_begin t_end],'YLim',[0 max(counts)]);
patch([vertcat(session(:).timerange)'; ...
    flipud(vertcat(session(:).timerange)')],...
    repmat([0; 0; max(get(gca,'YLim')); max(get(gca,'YLim'))], ...
    [1,numel(session)]),'r','EdgeColor','none');
hold on;
bar(bins,counts,1.0,'FaceColor','k','EdgeColor','k');
xlabel('timestamp');
ylabel('event counts');

% Interactively set maximum amplitude cutoff
while 1
  % Keep track of how many events will be censored with this cutoff
  num_rejected = nnz(any(max(abs(spike.samples),[],1) > ...
      repmat(amplitude_cutoff,[1 1 num_events]),2));
  disp(sprintf('Amplitude cutoff is %s',mat2str(amplitude_cutoff)));
  disp(sprintf('%d of %d events will be censored with this cutoff', ...
      num_rejected,num_events));
  confirm = input('Are you satisfied with this cutoff? [n]/y/q: ','s');
  if strcmp(confirm,'y')
    delete(fig);
    break
  elseif strcmp(confirm,'q');
    disp('Amplitude cutoff will be set to all zeros (reject all events)');
    delete(fig);
    amplitude_cutoff = zeros(1,NCHANNELS);
    return;
  else
    while 1
      try
        newval = input('Enter a new value for amplitude cutoff: ');
        if ( isnumeric(newval) && isreal(newval) && ...
            all(newval > 0) && ~any(isnan(newval)) && ...
            (isscalar(newval) || (isvector(newval) && ...
            (numel(newval) == NCHANNELS))) )
          break;
        end
      catch
        disp('That was not a valid value. Please try again: ');
      end
    end
    if isscalar(newval)
      amplitude_cutoff = newval*ones(1,NCHANNELS);
    else
      amplitude_cutoff = newval; 
    end
    for c = 1:NCHANNELS
      if ~isfinite(amplitude_cutoff(c))
        edges{c} = [linspace(0,MAX_SAMPLE,101)'; Inf];
      else
        edges{c} = [linspace(0,amplitude_cutoff(c),101)'; Inf];
      end
    end
    for i = 1:numel(pair)
      if pair(i)
        h(i) = subplot(NCHANNELS,NCHANNELS-1,pair(i),'replace');
        pixels = sparse(histcn([maxval(:,row(i)), ...
            maxval(:,col(i))], ...
            edges{row(i)},edges{col(i)}) > 0);
        % note that pixels is transposed when passed to imagesc because of the way
        % that imagesc orients the data for display
        imagesc([0 edges{row(i)}(end-1)], ...
            [0 edges{col(i)}(end-1)], ...
            pixels','Parent',h(i));
        set(h(i),'XLim',[0 edges{row(i)}(end-1)], ...
            'YLim',[0 edges{col(i)}(end-1)], ...
            'XTick',[0 edges{row(i)}(end-1)], ...
            'YTick',[0 edges{col(i)}(end-1)], ...
            'DataAspectRatio',[1 1 1],'XDir','normal','YDir','normal');
        xlabel(h(i),sprintf('channel %d',row(i)));
        ylabel(h(i),sprintf('channel %d',col(i)));
      end
    end
  end
end

