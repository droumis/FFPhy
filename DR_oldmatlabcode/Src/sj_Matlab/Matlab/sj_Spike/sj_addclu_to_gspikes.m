function sj_addclu_to_gspikes(datafile,clufile,sortmeth,converttomatclust)
% 
% sj_addclu_to_gspikes('spiketetcut_all_thr-100_reduced2_tetNF','spiketetcut_all_thr-100_reduced2_tetNF.clu.4','Kkwik');
% 
% Take clustering info from a .clu file (either from Kkwik or gsort) an Original "gspikes" file 
% Optional: make a "matspikes" and a "matclust" file structure 
%
% sortmeth = 'Kkwik' or 'G'
% datafile = gspikefile: just add a spikes.hierarchy.assigns
% structure to get gspikefile_sort and convert to get a new matclust file
%
%
% Shantanu, 04/21/2010

if nargin<3,
    sortmeth='Kkwik';
end

if nargin<4,
    converttomatclust=0;
end


% Get Name for Output file (sort files)
gspike_sortfile = [datafile '_' sortmeth 'sort_gspike'];

% Load data - Get nspikes
load(datafile);
nspikes = length(spikes.fstimes);

% Load clu file
clu = load(clufile);
nclu = clu(1);
assigns = clu(2:end);

% Check
if length(assigns)~=nspikes,
    disp('Error: No. of assignments does not equal number of spikes');
    return
end

%%% Output gspikefile_sort
% Put assignments in structure
spikes.hierarchy.assigns = assigns;
spikes.hierarchy.nclusters = nclu;
save(gspike_sortfile,'spikes');


%%% Output "matspikes", "params" file and "matclust" file
if converttomatclust==1,
    
    % 1st file
    timestamps=uint32(spikes.fstimes*10);
    waves=int16(zeros(size(spikes.waveforms_ch1,2),nch,size(spikes.waveforms_ch1,1)));
    waves(:,1,:)=int16(1000*spikes.waveforms_ch1');
    waves(:,2,:)=int16(1000*spikes.waveforms_ch2');
    waves(:,3,:)=int16(1000*spikes.waveforms_ch3');
    waves(:,4,:)=int16(1000*spikes.waveforms_ch4');
    
    save(matspike_file,'timestamps','waves');
    
    %% Have to check if param file exists or not. Instead, Just use
    % sj_gspikes_to_fetm BEFORE sorting
    % Converting to matclust is complicated. Instead: load param file and
    % import cluster assignments
    
end
