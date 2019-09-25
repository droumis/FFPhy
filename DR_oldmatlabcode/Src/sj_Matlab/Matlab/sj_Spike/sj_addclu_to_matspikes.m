function sj_addclu_to_matspikes(datafile,paramfile, clufile,sortmeth,converttogspikes)

% 
% sj_addclu_to_matspikes('02_022810-8','02_022810-8_params',02_022810.clu.8,'Kkwik',1);
% 
% Take clustering info from a .clu file (either from Kkwik or gsort) an Original "matspikes" + "params" file 
% and make a "matclust" AND a "gspikes" file structure 
%
% sortmeth = 'Kkwik' or 'G'
%
% if filetype=0, datafile = gspikefile: just add a spikes.hierarchy.assigns
% structure to get gspikefile_sort and convert to get a new matclust file
%
% if filetype=1 (DEFAULT), you start with a "param.mat" file and a "matspikes" file
%
%
% Shantanu, 04/21/2010

if nargin<3,
    sortmeth='Kkwik';
end

if nargin<4,
    converttogspikes=1;  % Default: make gspikefile also with sorted assignments
end

% Get Name for Output file (sort files)
gspike_sortfile = [datafile '_' sortmeth 'sort_gspike'];
matclustfile = ['matclust_' paramfile];

% Load data - Get nspikes
load(datafile);
nspikes = length(timestamps);

% Load clu file
clu = load(clufile);
nclu = clu(1);
assigns = clu(2:end);

% Check
if length(assigns)~=nspikes,
    disp('Error: No. of assignments does not equal number of spikes');
    return
end

% Load parameter file

load(paramfile);



%%%%% Save gspikefile: Easy

if converttogspikes==1
    
    spikes.waveforms_ch1=squeeze(waves(:,1,:));
    spikes.waveforms_ch2=squeeze(waves(:,2,:));
    spikes.waveforms_ch3=squeeze(waves(:,3,:));
    spikes.waveforms_ch4=squeeze(waves(:,4,:));
    spikes.Fs=Fs;
    spikes.swtimes=double(timestamps)./UnitsPerSec; % in sec
    spikes.fstimes=spikes.swtimes*1000;    % in ms
    spikes.ftimes=int16(spikes.fstimes);
    spikes.spiketimes=spikes.swtimes;
    spikes.threshT=9;
    spikes.threshV=[-Inf min(thresh)];
    spikes.waveforms=[spikes.waveforms_ch1, spikes.waveforms_ch2, spikes.waveforms_ch3, spikes.waveforms_ch4];
    
    spikes.hierarchy.assigns = assigns;
    spikes.hierarchy.nclusters = nclu;
    
    save(gspike_sortfile,'spikes');
end











