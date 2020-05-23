


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

raw = readmda('');
filt_masked = readmda('');

points_masked = sum(filt_masked(1,:)==0);
sec_masked = points_masked/3e4;
sec_data = size(filt_masked,2)/3e4;
sprintf('%.02f s masked of %.02f s', sec_masked, sec_data);

%% Recreate the masking threshold
% square and chunk the data (drop any extra that doesn't fit in a chunk)
chunk_size = 1000;
numchunks = floor(length(filt_unmasked_ch1)/chunk_size);
sq = filt_unmasked_ch1.^2;
chunks = reshape(sq(1:chunk_size*numchunks),chunk_size,numchunks);
% sum and sqrt each chunk
rss = sqrt(sum(chunks));
% calculate threshold based on mean and sd of rss
thresh = mean(rss) + 10 * std(rss);

%% viz 
% 1. Does the masking threshold correctly mask suprethreshold data?
% 2. Are SWR's incorrectly masked? By how much? is this consistent across 
%   ntrodes, data, animals?
% 3. What should the artifact mask threshold be instead? 

filt = readmda('');




