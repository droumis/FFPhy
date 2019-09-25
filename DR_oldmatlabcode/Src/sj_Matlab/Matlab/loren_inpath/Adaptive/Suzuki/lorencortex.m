function c = cortex (source_filename, possible_clusters)

% @cortex\cortex.m
% the constructor for the @cortex class.
%
% Syntax 1: COBJ = CORTEX (PATH\FILENAME, CLUSTERRANGE);
%
% The first argument is the name of the CORTEX data file;
% CLUSTERRANGE is a vector containing the range of 'possible' clusters;
%
% Syntax 2: COBJ = CORTEX (PATH\FILENAME);
%
% Loads to memory a pre-processed object, previously created with the 
% method SAVE
%
% Syntax 3: COBJ = CORTEX
%
% Declares an empty @cortex object
%
% See also: SPIKEMATRIX, TRIALFIRINGRATE, MEANFIRINGRATE, DENSITY
% Last modified: 23 Dec 98

switch nargin
case 0
   c = emptyobject;
   disp ('Empty @cortex object created');
case 1
   f = openfile (source_filename, 'r');
   c = loadump (f);
   fclose (f);
   disp (['File ''', source_filename, ''', preprocessed from ''', c.name, ''', loaded and ready.']);
case 2
   f = openfile (source_filename, 'r');
   % strip the path from the file name and make sure
   % that the remaining string is not longer than 30 chars
   pos = max(find(source_filename == '\' | source_filename == '/')) + 1;
   if isempty (pos)
      name = source_filename;
   else
      name = source_filename (pos:length(source_filename));
   end
   if length (name) > 30
      name = name (1:30);
   end
   c = loadraw (f, name, possible_clusters);
   fclose (f);
   disp (['File ''', source_filename, ''' preprocessed and ready.']);
otherwise
   error ('wrong number of arguments');
end

%c=class(c, 'cortex');
return % exit function

%========================================
function c = emptyobject
%========================================

c.name = 'noname';
c.trials = 0;
c.cond_no           = [];
c.repeat_no         = [];
c.block_no          = [];
c.eog_rate          = [];
c.KHz_resolution    = [];
c.trial_type        = [];
c.exp_response      = [];
c.response          = [];
c.times             = [];
c.codes             = [];
c.cluster_bank      = [];
c.spike_start       = [];
c.spike_end         = [];
c.spikes            = [];
return % to main

%========================================
function c = loadump (f)
%========================================

c.name = deblank(char(fread (f,30,     'char'  ))');
t                   = fread (f, 1,     'uint16');
c.trials            = t;
c.cond_no           = fread (f, t,     'uint16');
c.repeat_no         = fread (f, t,     'uint16');
c.block_no          = fread (f, t,     'uint16');
c.eog_rate          = fread (f, t,     'uint8' );
c.KHz_resolution    = fread (f, t,     'uint8' );
c.trial_type        = fread (f, t,     'uint16');
c.exp_response      = fread (f, t,     'uint16');
c.response          = fread (f, t,     'uint16');
m                   = fread (f, 1,     'uint16');
c.times             = fread (f, [t,m], 'uint32' );
c.codes             = fread (f, [t,m], 'uint16');
u                   = fread (f, 1,     'uint16'); % the number of different clusters in this file
c.cluster_bank      = fread (f, u,     'uint16'); % the list of cluster codes
c.spike_start       = fread (f, [t,u], 'uint32' );
c.spike_end         = fread (f, [t,u], 'uint32' );

last_spike  = max (max (c.spike_end));

% this is necessary in case the file contains no spikes
if isempty (last_spike)
   last_spike = 0;
end

% all the time stamps for all the clusters
c.spikes    = uint32 (fread (f, last_spike, 'uint32' ));

return % to main

%========================================
function c = loadraw (f, name, possible_clusters)
%========================================

i               = 0;
trials          = 0;
max_codes       = 0;
total_spikes    = 0;
head            = zeros (14, 1);
cluster_bank    = [];

c.name = name;
disp(['Scanning cortex data file ', c.name, ': pass 1']);

while 1
   % read the first header entry and use it to test for eof
   % and for file integrity (this value should always be 26)
   header_len = fread (f, 1, 'uint16');
   if feof (f)
      break % end of pass 1, skip to the next stage
   elseif header_len ~= 26
      error ('header size mismatch; corrupted data file?');
   else
      trials = trials + 1;
   end
   
   %read the rest of the header
   head ( 2: 9) = fread (f, 8, 'uint16');
   head (10:11) = fread (f, 2, 'uint8' );
   head (12:14) = fread (f, 3, 'uint16');
   
   % check for file integrity
   if head (6) ~= head (7) * 2
      error ('timebuf vs. codebuf mismatch; corrupted data file?')'
   end
   
   %skip the time stamps but get the 'code' buffer
   fseek (f, head(6), 'cof');
   events_this_trial = head (7)/2;
   
   codes = fread(f, events_this_trial, 'uint16');
   
   %skip eog and epp buffers
   fseek (f, head(8)+head(9), 'cof');
   
   %look for clusters in this trial
   found_clusters = intersect (codes, possible_clusters);
   cluster_bank   = union (cluster_bank, found_clusters);
   spikes_this_trial = sum (ismember  (codes, found_clusters));
   
   total_spikes = total_spikes + spikes_this_trial;
   max_codes = max (max_codes, events_this_trial - spikes_this_trial);
end

num_of_clusters_this_file = length (cluster_bank);
msg = sprintf ('%d cluster(s) found; %d total spikes in %d trials',...
   num_of_clusters_this_file, total_spikes, trials);
disp (msg);

%now the various arrays and variables can be initialized
%for the second pass through the file

c.trials            = trials;
c.cond_no           = zeros(trials, 1);
c.repeat_no         = zeros(trials, 1);
c.block_no          = zeros(trials, 1);
c.eog_rate          = zeros(trials, 1);
c.KHz_resolution    = zeros(trials, 1);
c.trial_type        = zeros(trials, 1);
c.exp_response      = zeros(trials, 1);
c.response          = zeros(trials, 1);
c.times             = zeros(trials, max_codes);
c.codes             = zeros(trials, max_codes);
c.cluster_bank      = cluster_bank;
c.spike_start       = zeros (trials, num_of_clusters_this_file);
c.spike_end         = zeros (trials, num_of_clusters_this_file);
c.spikes            = uint32 (zeros (total_spikes, 1));

disp(['Scanning cortex data file ', c.name, ': pass 2']);
frewind (f);
current_spike_pos = 1;

% trial-by-trial loop
for i=1:trials
   
   head ( 1: 9) = fread (f, 9, 'uint16');
   head (10:11) = fread (f, 2, 'uint8' );
   head (12:14) = fread (f, 3, 'uint16');
   
   c.cond_no(i)           = head (2) + 1;
   c.repeat_no(i)         = head (3) + 1;
   c.block_no(i)          = head (4) + 1;
   c.eog_rate(i)          = head (10);
   c.KHz_resolution(i)    = head (11);
   c.trial_type(i)        = head (12);
   c.exp_response(i)      = head (13);
   c.response(i)          = head (14);
   
   total_num_of_events  = head (7) /2;
   
   % read in time-stamps:
   all_times = fread (f, total_num_of_events, 'uint32');
   
   % read in behavioral codes and spikes
   all_codes = fread (f, total_num_of_events, 'uint16');
   non_spikes = ones (total_num_of_events, 1);
   
   for j=1:num_of_clusters_this_file
      
      %returns the indices of the current cluster
      temp_cl = find (all_codes == cluster_bank(j));
      if ~isempty (temp_cl)
         %removes this cluster from the non-spike events
         non_spikes (temp_cl) = 0;
         
         %returns an array containing the time stamps of the corresponding spike
         temp_sp = nonzeros (all_times (temp_cl));
         
         %updates the pointers that will be needed to retrieve the relevant
         %spikes from the HUGE spikes array
         c.spike_start (i, j) = current_spike_pos;
         current_spike_pos    = c.spike_start (i, j) + length (temp_sp);
         c.spike_end   (i, j) = current_spike_pos - 1; 
         c.spikes (c.spike_start (i,j) : c.spike_end (i, j)) = temp_sp;
      end
   end
   
   %all that is left after removing spikes are behavioral events
   behav_indices                          = find (non_spikes);
   num_of_behav_events                    = length (behav_indices);
   
   if num_of_behav_events
      c.times (i, 1:num_of_behav_events) = all_times (behav_indices)';
      c.codes (i, 1:num_of_behav_events) = all_codes (behav_indices)';
   end
   % skip eog and epp buffers   
   fseek (f, head (8) + head (9), 'cof');
   
end % of pass 2
return % to main
