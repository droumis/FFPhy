function [spikes,aspkfile]=sss_splitspkfileaftersort_1(filename, savename)

%sss_splitspkfileaftersort_1('joinedsortfile_tetNF-dejit-energy32subsortf','MIC18-G7R3-08-24-07');
%%Split files after sorting (files were joined with sss_joinspikefile_cont5

%% Assume that no file will run longer than 10000secs/ safer/ use
%% 10^5+swtimes, and 10^8*fstimes for each file and push

%% For joining spk files in multiple directories before sorting. By
%% default, you have to adjust spike times to prevent overlap of spktimes
%% during sorting
%% POST Sort, deconstruct the spikesagain based on fstimes, and 

%% sss_splitspkfileaftersort_1('joinedsortfile-energy32subsortf');

%%%%%%%%%%%%  NO NEED TO SAVE WAVEFORMS - SAVE SPACE
if nargin<2, savename=filename; end

load(filename); ospikes=spikes; spikes=[];
reftimes = ospikes.swtimes;

%%%%%%%%%%%%% Get Multipliers from swtimes
ntimes=mod(reftimes, 10^5);
multipliers = (reftimes-mod(reftimes, 10^5))/10^5;
unqmul=unique(multipliers);

for n=1:length(unqmul)    
    mul=unqmul(n);
    curridxs = (find(multipliers==mul));
    spikes.Fs=ospikes.Fs;
    spikes.threshT=ospikes.threshT;
    spikes.threshV=ospikes.threshV;

    spikes.swtimes = ntimes(curridxs);  
    spikes.fstimes = ntimes(curridxs)*1000;
    spikes.ftimes=ospikes.ftimes(curridxs);
    spikes.spiketimes=ospikes.spiketimes(curridxs);
    spikes.hierarchy.assigns=ospikes.hierarchy.assigns(curridxs);
    
    if isfield(ospikes,'nsecs_1'),
        cmd=sprintf('spikes.nsecs = ospikes.nsecs_%d',mul); eval(cmd);
    end
    
    if isfield(ospikes,'nsecs_reduced'),
        cmd=sprintf('spikes.nsecs_reduced = ospikes.nsecs_reduced_%d',mul); eval(cmd);
        cmd=sprintf('spikes.nsecs_reduced_factor = ospikes.nsecs_reduced_factor_%d',mul); eval(cmd);
    end
    
    save([savename '-' num2str(mul)], 'spikes');
    
end





