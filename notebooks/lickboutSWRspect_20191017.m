


%{
Is the riptriglfp TFpower different between 
lickbout and nonlickbout trials?
licbout correct and lickbout error trials?

- complex morlet wavelet method
- multitaper method

- compare as two conditions in a designmat 

%}


%% load riptrig lickbout lfp
% swr times, lickbout times, ca1 lfp
% gather the riptrig lfp data into blocks 

%% create wavelets focused around the ripple range

%% get fft of wavelets

%% concat all the data and compute fft

%% for each wavelet fft, conv in freq space with data fft

%% do ifft on result and trim the neg freqs and 1/2+1 wavelet length

%% reshape the data back to tensor form 

%% get magnitude of real component as power and compute mean over rips 

%% log transform mean power 

%% plot