
%% first analysis: power spectrum of each tetrode, each epoch (one day)

% input is a struct containing tetrode fields
% eeg{something}{day}{tetrode}.<field>

params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];

%% note spectrums is a 3D output matrix (tetrode,spectrum,epoch)
%% meanspectrum is a 2D output matrix (tetrode,meanspectrum)

params.Fs = eeg{1}{1}{1}.samprate;
spectra=cell(1,7,14);

for d=1
    for e=1:7                                             % 7 epochs / day
        for t=1:14  
            [spectra{d,e,t},f,Serr]=mtspectrumc(double(eeg{d}{e}{t}.data(1:450000)),params);   %% only analyzing 11 minutes / epoch
        end
    end
end

spectra=cell2mat(spectra);               % coefficents x epochs x tetrodes
meandayspectra = squeeze(mean(spectra,2));  % coefficients x tetrodes


%% ratiospectra (to day mean)

ratiospectra=nan(size(spectra,1),7,14); % initialize
for t=1:14
    ratiospectra(:,:,t) = bsxfun(@rdivide,spectra(:,:,t),meandayspectra(:,t));
end

smoothwindow=size(spectra,1)/20;
smoothkernel = gausskernel(smoothwindow/4,smoothwindow/6);   % Gaussian

for t=1:14

    figure
    
for e=[7]          % sleep?
    hold on;
    semilogy(f,conv(ratiospectra(:,e,t),smoothkernel,'same'),'color',jet(1)/8);
end

for e=[2]
    semilogy(f,conv(ratiospectra(:,e,t),smoothkernel,'same'),'g');

for e=3:6          % run?
    hold on;
    semilogy(f,conv(ratiospectra(:,e,t),smoothkernel,'same'),'r');
end
end
end








