function out = dfa_XfreqCoupling(index, dummy, eeg, varargin)

statespec_descript = '';
phasenbin = 18; % default # of circular phase bins for each phase frequency: 18
selection_phasetet = 0; %ntrode selected as the phase reference; 0 if same ntrode
PhaseFreqVector = [];
AmpFreqVector = [];
PhaseFreq_BandWidth = [];   
AmpFreq_BandWidth = [];

% % Set options
if (~isempty(varargin))
    assign(varargin{:});
end

% identify animal directory
animalinfo = animaldef(animal);
animaldir = animalinfo{2};
animalprefix = animalinfo{3};
task = loaddatastruct(animaldir,animalprefix,'task');
tetinfo = loaddatastruct(animaldir,animalprefix,'tetinfo');

day = index(1);
ep = index(2);
tet = index(3);

disp(['XFC: ' num2str([day ep tet]) ])

position=zeros(1,phasenbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/phasenbin;
for j=1:phasenbin
    position(j) = -pi+(j-1)*winsize;
end

eeg_amp = eeg{index(1)}{index(2)}{index(3)}.data;
samprate = round(eeg{index(1)}{index(2)}{index(3)}.samprate);
timevec = geteegtimes(eeg{index(1)}{index(2)}{index(3)});
data_length = length(eeg_amp);

%% check if there another ntrode selected as the phase reference, otherwise use same ntrode
if selection_phasetet > 0
    % load phase eeg
    % first identify reference retrode
    regnum = selection_phasetet;
    phasetet = task{day}{ep}.STA_reftet(1,regnum);
    % load eeg data
    edata = loadeegstruct(animaldir,animalprefix,'eeg',day,ep,phasetet);
    eeg_phs = edata{day}{ep}{phasetet}.data;
    timevec2 = geteegtimes(edata{day}{ep}{phasetet});
    if length(timevec) ~= length(timevec2)   ||   ~all(timevec == timevec2)
        % interpolate to match the trace ends
        eeg_phs = interp1(timevec2,eeg_phs,timevec,'nearest');
    end
    clear timevec2 edata
else
    phasetet = tet;
    eeg_phs = eeg_amp;
end
clear eeg
%% Filter and Hilbert transform
tic
AmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);
disp('filtering for amplitude envelope')
parfor iAfreq=1:length(AmpFreqVector)
    Af1 = AmpFreqVector(iAfreq);
    Af2=Af1+AmpFreq_BandWidth;
    AmpFreq=eegfilt(eeg_amp',samprate,Af1,Af2); % just filtering
    AmpFreqTransformed(iAfreq, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
        disp(sprintf('Amp filtering %d - %d Hz (%d of %d)',Af1, Af2, iAfreq, length(AmpFreqVector)))
%     disp('.');
end
disp('filtering for phase time series')
parfor iPfreq=1:length(PhaseFreqVector)
    Pf1 = PhaseFreqVector(iPfreq);
    Pf2 = Pf1 + PhaseFreq_BandWidth;
    PhaseFreq=eegfilt(eeg_phs',samprate,Pf1,Pf2); % this is just filtering
    PhaseFreqTransformed(iPfreq, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
        disp(sprintf('Phase filtering %d - %d Hz (%d of %d)',Pf1, Pf2, iPfreq, length(PhaseFreqVector)))
%     disp('.');
end
out.toc = toc; %save timer

%% Iterate over states, first getting included periods then calculating XFC

excludeperiods = cell(1,length(statespec));
duration = nan(1,length(statespec));
lenstatespec = length(statespec);
for istate = 1:lenstatespec
    XFCall{istate} = single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
    % TIMEFILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timefilter = statespec{istate};
    for itimefilt = 1:length(timefilter)
        output = feval(timefilter{itimefilt}{1}, animaldir, animalprefix, [day ep], timefilter{itimefilt}{3:end});
        output2 = evaluatefilter2(output, timefilter{itimefilt}{2});
        tmpexcludeperiods = getExcludePeriods(output2{day}{ep}(:,1),output2{day}{ep}(:,2)) ;
%         if isempty(excludeperiods)
%             excludeperiods{istate} = tmpexcludeperiods;
%         else
            excludeperiods{istate} = combineExcludePeriods(excludeperiods{istate}, tmpexcludeperiods);
%         end
        clear tmpexcludeperiods
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply excludetimes to filtered eeg
    includedindices = ~isExcluded(timevec,excludeperiods{istate});
    duration(istate) = sum(includedindices)/samprate;
    disp([num2str(day) ' ' num2str(ep) ' ' num2str(tet) ' state ' num2str(istate) ' : ' num2str(duration(istate)) ' seconds'])
    % time-filter the filtered data
    AmpFreqTransformed_filtered = AmpFreqTransformed(:,includedindices);
    PhaseFreqTransformed_filtered = PhaseFreqTransformed(:,includedindices);
    
    % Do comodulation calculation
    counter1=0;
    lenPhaseFreqVec = length(PhaseFreqVector);
    lenAmpFreqVec = length(AmpFreqVector);
    for iPfreq=1:lenPhaseFreqVec
%         counter1=counter1+1;
        Pf1 = PhaseFreqVector(iPfreq);
        Pf2 = Pf1+PhaseFreq_BandWidth;
%         counter2=0;
        for iAfreq=1:lenAmpFreqVec
%             counter2=counter2+1;
            Af1 = AmpFreqVector(iAfreq);
            Af2 = Af1+AmpFreq_BandWidth;
            [XFC,MeanAmp]=ModIndex_v2(PhaseFreqTransformed_filtered(iPfreq, :), AmpFreqTransformed_filtered(iAfreq, :), position) ;
            XFCall{istate}(iPfreq,iAfreq) = XFC;
%             disp([num2str(day) ' ' num2str(ep) ' ' num2str(tet) ' state ' num2str(istate) ' XFC score: ' num2str(XFC)])
        end
        disp(sprintf('state %d of %d, phase %d of %d complete', istate,lenstatespec, iPfreq, lenPhaseFreqVec));
    end
end

out.date = date;
out.index = index ;
out.phasetet = phasetet;
out.statespec = statespec;
out.statespec_descript = statespec_descript;
out.duration = duration ;
out.XFC = XFCall;
out.frequency_amplitude = AmpFreqVector ;
out.frequency_phase = PhaseFreqVector ;
out.phasefreq_bandwidth = PhaseFreq_BandWidth ;
out.ampfreq_bandwidth = AmpFreq_BandWidth ;
end