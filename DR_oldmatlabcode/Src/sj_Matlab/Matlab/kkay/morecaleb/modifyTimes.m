function [names, ranges] = correctTimes(times_mat_filename, epoch, label, beginning, ending)
% can be called as: function [names, ranges] = correctTimes(times_mat_filename, epoch, label, beginning, ending)
% or
% function [names, ranges] = correctTimes(epoch, label, beginning, ending)
% 
% Empty vectors/strings for label, beginning, and ending leave current values.

if nargin == 4
  ending = beginning;
  beginning = label;
  label = epoch;
  epoch = times_mat_filename;
  times_mat_filename = 'times.mat';

  if ~isnumeric(epoch) || ~ischar(label) 
    epoch
    label
    error('Improper syntax - use "help correctTimes" for help.');
  end
end

load(times_mat_filename);

if epoch > size(ranges,1)-1
  error('Epoch not found in times.mat file.');
end

tokens = regexp(names{epoch+1},...
  '(\d+)\s+(.*)\s+(\d+):(\d+):*(\d+)\.*(\d*)\s*-\s*(\d+):(\d+):*(\d+)\.*(\d*)','tokens');

current_label = tokens{1}{2};
if ~isempty(label)
  current_label = label;
end


if isempty(tokens{1}{6})
  beginning_str = sprintf('%d:%02d:%02d',str2num(tokens{1}{3}),str2num(tokens{1}{4}),str2num(tokens{1}{5}));
else
  beginning_str = sprintf('%d:%02d:%02d.%d',str2num(tokens{1}{3}),str2num(tokens{1}{4}),str2num(tokens{1}{5}),str2num(tokens{1}{6}));
end
beginning_ts = ranges(epoch+1,1);

if isempty(tokens{1}{10})
  ending_str = sprintf('%d:%02d:%02d',str2num(tokens{1}{7}),str2num(tokens{1}{8}),str2num(tokens{1}{9}));
else
  ending_str = sprintf('%d:%02d:%02d.%d',str2num(tokens{1}{7}),str2num(tokens{1}{8}),str2num(tokens{1}{9}),str2num(tokens{1}{10}));
end
ending_ts = ranges(epoch+1,2);


if ~isempty(beginning) && ischar(beginning)
  beginning_str = beginning;
  beginning_ts = str2ts(beginning);
  if beginning_ts < ranges(epoch+1,1)
    warning('Specified beginning is before entry in times.mat file.');
  end
elseif ~isempty(beginning) && isnumeric(beginning)
  beginning_ts = beginning;
  if beginning_ts < ranges(epoch+1,1)
    warning('Specified beginning is before entry in times.mat file.');
  end
  beginning_str = timestampToString(beginning_ts);
end

if ~isempty(ending) && ischar(ending)
  ending_str = ending;
  ending_ts = str2ts(ending);
  if ending_ts > ranges(epoch+1,2)
    warning('Specified ending is after entry in times.mat file.');
  end
elseif ~isempty(ending) && isnumeric(ending)
  ending_ts = ending;
  if ending_ts > ranges(epoch+1,2)
    warning('Specified ending is after entry in times.mat file.');
  end
  ending_str = timestampToString(ending_ts);
end

ranges(epoch+1,:) = [beginning_ts ending_ts];
names{epoch+1} = sprintf('%d  %s %s-%s', epoch+1, current_label, beginning_str, ending_str);

save(times_mat_filename,'names','ranges');


function timestring = timestampToString(timestamp)

% DO NOT TAMPER WITH THIS CONSTANT CONVERSION FACTOR!
TS_PER_SEC = 1e4;

if ~isnumeric(timestamp)
  error('TIMESTAMP argument must be a number');
end

% perform calculations in double precision, because MATLAB behaves strangely
% when doing integer arithmetic (e.g. does not always truncate in the expected
% ways)
timestamp = double(timestamp);
xxxx = rem(timestamp,TS_PER_SEC);
ss = rem(floor(timestamp/TS_PER_SEC),60);
mm = rem(floor(timestamp/TS_PER_SEC/60),60);
hhh = floor(timestamp/TS_PER_SEC/3600);

timestring = sprintf('%d:%02d:%02d.%04d',hhh,mm,ss,xxxx);

