


%{
Check to what extent my ripple detection and associated analyses are
impacted by the mountainsort masking of 'artifacts'
Anna shared her code used to asses her data:
-> trodes2ff_dr/annas_mask_debug.m  
Anna's appraoch:
- load segment of raw.mda and filt.mda
    - she is loading filt masked, unmasked from tmp dir
        - where is she getting the hardpath filenames from? maybe prv file?
        - i doubt i still have any tmp files.. so i'd need to recreate.. i
        guess one ntrode chan per day might be enough to get a better
        threshold..? maybe? or would each chan need their own threshold?
            - could open a raw mda in mountainview and recreate filt
            - could rerun a sort
- plot the raw, filtered unmasked, filtered masked, mask threshold
- check duration of loaded segment that was set to 0 with 
%}

% which ca1 ntrode chan with cells, ripples?
% where do i keep this info? 

raw = readmda('/media/droumis/DR_swapdata9/JZ4/mountainlab_output/20170424/nt28/raw.mda');                          
filt = readmda('/media/droumis/DR_swapdata9/JZ4/mountainsort_tmp/output_82339955944c4ec09e23f9941d8f3234cf37c278_timeseries_out.mda');
%%
chunk_len = 1e6;
raw_ch1 = raw(1,chunk_len:chunk_len*2);
filt_ch1 = filt(1,chunk_len:chunk_len*2);
clear filt raw
%%
points_masked = sum(filt_ch1(1,:)==0);
sec_masked = points_masked/3e4;
sec_data = size(filt_ch1,2)/3e4;
fprintf('%.02f s masked of %.02f s \n', sec_masked, sec_data);
fprintf('%.03f pct of data masked\n', (sec_masked/sec_data)*100);

% JZ4 day 3 nt 28 ch1: <1 % of data masked

%% Recreate the masking threshold
% square and chunk the data (drop any extra that doesn't fit in a chunk)
chunk_size = 1000;
numchunks = floor(length(filt_ch1)/chunk_size);
sq = filt_ch1.^2;
chunks = reshape(sq(1:chunk_size*numchunks),chunk_size,numchunks);
% sum and sqrt each chunk
rss = sqrt(sum(chunks));
% calculate threshold based on mean and sd of rss
meanrss = mean(rss);
thresh = mean(rss) + 3 * std(rss);
% how the fuck is this the threshold?? i dont understand

%% viz 
figure(1)
set(gcf);
% ax2 = subplot(3,1,2); hold on;
hold off
plot(raw_ch1(1:chunk_len),'k'); %title('unmasked filt')
hold on
plot(filt_ch1(1:chunk_len),'b'); %title('unmasked filt')
plot([1 chunk_len],[thresh thresh],'r:', 'linewidth', 4)
plot([1 chunk_len],[-thresh -thresh],'r:', 'linewidth', 4)
plot([1 chunk_len],[meanrss meanrss],'m:', 'linewidth', 4)
%%
% 1. Does the masking threshold correctly mask suprethreshold data?
% 2. Are SWR's incorrectly masked? By how much? is this consistent across 
%   ntrodes, data, animals?
% 3. What should the artifact mask threshold be instead? 



