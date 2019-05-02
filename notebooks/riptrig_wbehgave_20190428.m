
%{


load spikes, rips, perform

get spikes in rip as incr process +-.5 s

label rip spikes with vec of ripnum or ripstarttime

perform change vs rip-spike 0lag xcorr

%}
%% choose ldr, seg, and brsh

animal = 'D13';
andef = animaldef(animal);
day = 6;
ep = 2;
tetA = 9;
tetB = 4;
clustA = 9;
clustB = 6;
%% load the pins
spikes = loaddatastruct(andef{2}, animal, 'spikes', day);
rips = loaddatastruct(andef{2}, animal, 'ca1rippleskons', day);
perform = load([andef{2} animal 'ca1rippleskons07.mat']);
cellinfo = loaddatastruct(andef{2}, animal, 'cellinfo', day);
ripeeg = loadeegstruct(andef{2}, animal, 'ripple', day, ep, [tetA tetB]);
eeg = loadeegstruct(andef{2}, animal, 'eeg', day, ep, [tetA tetB]);

%%  abbreviate
rips_start_end = [rips{day}{ep}{1}.starttime rips{day}{ep}{1}.endtime];

spikesA_list = spikes{day}{ep}{tetA}{clustA};
spikesB_list = spikes{day}{ep}{tetB}{clustB};
clustA_info = cellinfo{day}{ep}{tetA}{clustA};
clustB_info = cellinfo{day}{ep}{tetB}{clustB};

tetAripmag = ripeeg{day}{ep}{tetA}.data(:,3); %field 3 is envelope_magnitude
tetBripmag = ripeeg{day}{ep}{tetB}.data(:,3); %field 3 is envelope_magnitude
srate = eeg{day}{ep}{tetA}.samprate;

%% filter mask on sp for rip st end
ripspikesA = spikesA_list.data(logical(~isExcluded(spikesA_list.data, rips_start_end)));
ripspikesB = spikesB_list.data(logical(~isExcluded(spikesB_list.data, rips_start_end)));

%% reconstitute eeg times
eegtime = eeg{day}{ep}{tetA}.starttime:1/srate:eeg{day}{ep}{tetA}.endtime;

%% bin rip spikes to lfp time
% this is basical creating an increment process from spike times

ripspikesAeegtime = histcounts(ripspikesA, eegtime);
ripspikesBeegtime = histcounts(ripspikesB, eegtime);

%% xcorr of the entire spiketrainA v B..

zeroXc = xcorr(ripspikesAeegtime, ripspikesBeegtime, 0, 'coeff');
[xc, xclags] = xcorr(ripspikesAeegtime, ripspikesBeegtime, 750, 'coeff');

%% get xcorr for each event
% what if there aren't sufficent spikes for a given unit?? i.e. what if it
% wasn't valid for a certain period of time.. 
% this is why i wanted to combine spikes into marks for each tetrode. 
for rip = 1:size(rips_start_end, 1)
    ripspikesA
end

%% plotxcorr

% use the xcorr from the two spike trains to 

plot(xclags, xc)


%% 
%{
the reason why i wanted to plot the rips and rip trig spikes before was to
see if it was even reasonable to use single units for rip trig xc.. i.e.
there needs to be valid clusters across the whole epoch period for the xc
to even make sense.. otherwise it's just 00s.. 

D10 seems to have a couple reasonable rip trig su spikes.. but not all of
them and not very many compared to the num of tets.. so maybe just go
straight to mu?
%}














































