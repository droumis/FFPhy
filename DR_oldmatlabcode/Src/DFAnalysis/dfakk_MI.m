function out = dfakk_MI(index, excludetimes, eeg, varargin)

%Define sampling rate in Hz

% Define vector of frequencies whose amplitudes are modulated
srate = 1500;
%PhaseFreqVector=0:2:50;
%AmpFreqVector=10:5:300;
PhaseFreqVector=0:.5:50;
AmpFreqVector=10:5:300;


PhaseFreq_BandWidth=4;
AmpFreq_BandWidth=10;

nbin = 18;
    position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
    winsize = 2*pi/nbin;
    for j=1:nbin 
        position(j) = -pi+(j-1)*winsize; 
    end

%Set options
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'nbin'
                nbins = varargin{option+1};            
            case 'frequency_amplitude'
                frequency_amplitude = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end


e = eeg{index(1)}{index(2)}{index(3)}.data;
timevec = geteegtimes(eeg{index(1)}{index(2)}{index(3)});
data_length = length(e);
clear eeg
m = [];

%% Filter and Hilbert transform
'CPU filtering'
tic
Comodulogram=single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
AmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);

for ii=1:length(AmpFreqVector)
    Af1 = AmpFreqVector(ii);
    Af2=Af1+AmpFreq_BandWidth;
    AmpFreq=eegfilt(e',srate,Af1,Af2); % just filtering
    AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
end

for jj=1:length(PhaseFreqVector)
    Pf1 = PhaseFreqVector(jj);
    Pf2 = Pf1 + PhaseFreq_BandWidth;
    PhaseFreq=eegfilt(e',srate,Pf1,Pf2); % this is just filtering 
    PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
end
toc


%% Apply excludetimes to filtered eeg

includedindices = ~isExcluded(timevec,excludetimes);
AmpFreqTransformed = AmpFreqTransformed(:,includedindices);
PhaseFreqTransformed = PhaseFreqTransformed(:,includedindices);


%% Do comodulation calculation
'Comodulation loop'

counter1=0;
for ii=1:length(PhaseFreqVector)
counter1=counter1+1;

    Pf1 = PhaseFreqVector(ii);
    Pf2 = Pf1+PhaseFreq_BandWidth;
    
    counter2=0;
    for jj=1:length(AmpFreqVector)
    counter2=counter2+1;
    
        Af1 = AmpFreqVector(jj);
        Af2 = Af1+AmpFreq_BandWidth;
        [MI,MeanAmp]=ModIndex_v2(PhaseFreqTransformed(ii, :), AmpFreqTransformed(jj, :), position);
        Comodulogram(counter1,counter2)=MI;
    end
end
toc

figure
contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulogram',30,'lines','none')
set(gca,'fontsize',14)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
colorbar


out.index = index;
out.comodulogram =Comodulogram;
out.frequency_amplitude = AmpFreqVector;
out.frequency_phase = PhaseFreqVector;
out.phasefreq_bandwidth = PhaseFreq_BandWidth;
out.ampfreq_bandwidth = AmpFreq_BandWidth;

end