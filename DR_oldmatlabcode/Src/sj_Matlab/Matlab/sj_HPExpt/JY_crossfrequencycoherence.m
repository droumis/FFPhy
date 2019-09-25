function out = JY_crossfrequencycoherence(index, excludetimes, eeg,meaneegspectrograms, varargin)

%Define sampling rate in Hz
frequency_amplitude = [20:2:100];

% parse the options
params = [];
params.Fs = 1500;
params.fpass = [1 100];
params.err = [1 0.05];
window =[1 0.5];
params.tapers=[3 5];

%Set options
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'frequency_amplitude'
                frequency_amplitude = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% if using referenced eeg signal
e = eeg{index(1)}{index(2)}{index(4)}.data;
estart = eeg{index(1)}{index(2)}{index(4)}.starttime;
clear eeg

% if using eeg with reference added back
e = meaneegspectrograms{index(1)}{index(2)}{index(4)}.data;
estart = meaneegspectrograms{index(1)}{index(2)}{index(4)}.starttime;
clear meaneegspectrograms

m = [];
%Compute coherence between the raw signal and each gamma filtered trace
for g = frequency_amplitude(1:end)
    %Compute amplitude time series for given gamma frequency
    filterstring = sprintf('/home/mcarr/Src/Matlab/Filters/cfcampfilt%s.mat',num2str(g));
    eval(['load ', filterstring]);
    amplitude = abs(hilbert(filtfilt(cfcampfilt,1,e)));
           
    %Compute coherence for each gamma frequency
    [C,phi,S12,S1,S2,t,f]=cohgramc(amplitude,e,window, params);

    %apply excludetimes
    if isempty(m)
        omit = ~isExcluded(t+estart, excludetimes);
        m = zeros(length(frequency_amplitude),length(f));
    end
    coherence = mean(C(omit==0,:));
    m(find(frequency_amplitude==g),:) = coherence;
end


  %imagesc(t-win(1),f,S(:,:,ii)');
                     imagesc(f,frequency_amplitude,m,[0.4 0.5]);
                    
                    %imagesc(tmt-win(1),fmt,epocheegmean',[0, 20]);
                        axis xy;
                    %title(string);
                    colorbar

% 


out.cfc =m;
out.frequency = f;
out.amp_freq=frequency_amplitude;
out.index = index;
end