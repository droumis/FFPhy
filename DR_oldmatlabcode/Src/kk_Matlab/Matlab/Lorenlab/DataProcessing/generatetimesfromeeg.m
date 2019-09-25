function [names, ranges] = generatetimesfromeeg(mindelay)
% function generatetimesfromeeg(mindelay)
% Script to create a generic times.mat file from eeg files.
% If mindelay is specified, combine all periods separated by less
%   than mindelay.

if nargin > 0
   mindelay = mindelay * 10; % convert to 0.1 ms ticks
end

d = dir('*.eeg');
if isempty(d)
  error('No EEG files found.');
end
epochs = epoch_extract(d(1).name,0);

n_epochs = size(epochs,1)

if (nargin > 0) & (n_epochs > 1)
  i = 1;
  while 1
     if (epochs(i+1,1) - epochs(i,2)) < mindelay
        epochs(i,2) = epochs(i+1,2);
        epochs(i+1,:) = [];
     else
        i = i + 1;
     end
     if (i > size(epochs,1)-1)
        break
     end
  end
   n_epochs = size(epochs,1)
end

for n = 2:(n_epochs+1)
  names{n} = sprintf('%d  epoch %d %s-%s', n, n-1,...
    cell2mat(timetrans(epochs(n-1,1),10000,1)), ...
    cell2mat(timetrans(epochs(n-1,2),10000,1)));
  ranges(n,:) = epochs(n-1,:);
end

names{1} = '1» All points';
ranges(1,:) = [0 ranges(n,2)+1000];
for n = size(ranges,1)+1:32
   names{n} = num2str(n);
end

if nargout == 0
  save times names ranges
end
