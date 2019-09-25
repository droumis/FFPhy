function [out] = gethighthetatimes(animaldir,animalprefix, epochs, tetlist, varargin)
% out = gethighthetatimes(animaldir,animalprefix,epochs, tetlist, options)
%
%    animaldir and animal prefix are strings indicating the base director for
%    the animal's data and the prefix for the data files
%
%    epochs is an Nx2 list of days and epochs
%
%    tetlist is a list of tetrodes to use or an empty matrix if the
%    'cellfilter' option is used.
%
% options are
%  'cellfilter', 'cellfilterstring'
%         specifies a cell filter to select the tetrodes to use for
%         high theta filtering
%  'minthresh', minthresh 
%         specifies a minimum threshold in stdev units for a valid 
%        high theta event  (default 0)
%
% Produces a cell structure with a time field and an nhighthetatimes field which
% indicates the number of electrodes with a highthetaple at each time point
%
% Examples:
% gethighthetatimes('/data/name/Fre', 'fre', epochs, 1)
% gethighthetatimes('/data/name/Fre', 'fre', epochs, [], 'cellfilter', '(isequal($area, ''CA1''))')

% assign the options
cellfilter = '';
minenergy = 0;
minthresh = 0;

[otherArgs] = procOptions(varargin);

if nargin < 4
   tetlist = [];
end

%check to see if a cell filter is specified
if (~isempty(cellfilter))
   % this will cause us to ignore tetlist
   cellinfo = loaddatastruct(animaldir, animalprefix, 'cellinfo');
end

loaddays = unique(epochs(:,1));
hightheta = loaddatastruct(animaldir, animalprefix, 'hightheta', loaddays);
for i = 1:size(epochs,1)
   % if cellfilter is set, apply it to this day and epoch
   if (~isempty(cellfilter))
      tetlist =  evaluatefilter(cellinfo{epochs(i,1)}{epochs(i,2)}, ...
            cellfilter); 
      % get rid of the cell indeces and extract only the tetrode numbers 
      tetlist = unique(tetlist(:,1))';
   end

   if isempty(tetlist)
      tetlist = 1:length(hightheta{epochs(i,1)}{epochs(i,2)});
      for j = length(tetlist):-1:1
         if isempty(hightheta{epochs(i,1)}{epochs(i,2)}{j})
            tetlist(j) = [];
         end
      end
   end
   % go through the tetlist and construct an an array where each element 
   % represents the number of active tetrodes for each 1 ms timestep.
   try
      ht = hightheta{epochs(i,1)}{epochs(i,2)}{tetlist(1)};
   catch
      keyboard
   end

   thtimes = ht.timerange(1):0.01:ht.timerange(end);
   nhightheta = zeros(size(thtimes));
   for t = 1:length(tetlist)
      tmphightheta = hightheta{epochs(i,1)}{epochs(i,2)}{tetlist(t)};
      % apply the minthresh threhsold
      htvalid = find(tmphightheta.maxthresh > minthresh);
      httimes = [tmphightheta.starttime(htvalid) tmphightheta.endtime(htvalid)];
      % create another parallel vector with bordering times for zeros
      nhttimes = [(httimes(:,1) - 0.00001) (httimes(:,2) + 0.00001)];
      httimes = reshape(httimes', length(httimes(:)), 1); 
      httimes(:,2) = 1;
      nhttimes = [ht.timerange(1) ; reshape(nhttimes', ...
         length(nhttimes(:)), 1) ; ht.timerange(2)];
      nhttimes(:,2) = 0;
      % create a new list with all of the times in it
      tlist = sortrows([httimes ; nhttimes]);
      % use interp to create a set of ones and zeros for each time
      % and add to nhightheta to get a cumulative count of the number of
      % highthetatimes per timestep
      try
         nhightheta = nhightheta + interp1(tlist(:,1), tlist(:,2), thtimes, 'nearest');
      catch
         keyboard
      end
   end
   out{epochs(i,1)}{epochs(i,2)}.time = thtimes';
   clear times;
   out{epochs(i,1)}{epochs(i,2)}.nhightheta = nhightheta';
end
