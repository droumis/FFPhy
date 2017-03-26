function out = dfakk_MI(index, dummy, eeg, target_det, animal, statespec, statespec_descript, varargin)

if ~isempty(target_det)
    if ~rowfind(index,target_det)
        out.dummy = [];
        return
    end
end

d = index(1);
ep = index(2);
tet = index(3);

disp(['MI: ' num2str([d ep tet]) ])

% default # of circular phase bins for each phase frequency: 18
nbin = 18;
    position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
    winsize = 2*pi/nbin;
    for j=1:nbin 
        position(j) = -pi+(j-1)*winsize; 
    end

selection_phasetet = 0;
stareftet_only = 0;

% Set options
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'nbin'
                nbins = varargin{option+1};            
            case 'frequency_amplitude'
                frequency_amplitude = varargin{option+1};
            case 'PhaseFreqVector'
                PhaseFreqVector = varargin{option+1};             
            case 'AmpFreqVector'
                AmpFreqVector = varargin{option+1};   
            case 'PhaseFreq_BandWidth'
                PhaseFreq_BandWidth = varargin{option+1}; 
            case 'AmpFreq_BandWidth'
                AmpFreq_BandWidth = varargin{option+1};      
            case 'selection_phasetet'
                selection_phasetet = varargin{option+1};    
            case 'stareftet_only'
                stareftet_only = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% identify animal directory
animalinfo = animaldef(animal);
animaldir = animalinfo{2};
animalprefix = animalinfo{3};
    cd([animaldir 'MIDATA_new/']);
    task = loaddatastruct(animaldir,animalprefix,'task');
    tetinfo = loaddatastruct(animaldir,animalprefix,'tetinfo');
    
% check here if you've already calculated MI w/ the set of states for this adet (an-day-ep-tet)
    % to do this, load the file (if any) and see if the statespec_descript string matches
    
    filename = dir(sprintf('*MI-%d-%d-%d.mat',d,ep,tet));
    concat_flag = 0;
    if ~isempty(filename)
        load([filename.name],'out')
        for ee = 1:length(out)
            if strcmp(statespec_descript,out(ee).statespec_descript)
                clear out
                out.dummy = [];
                return
            end
        end
        concat_flag = 1;  % indicates that we're going to concatenate a dataentry for this adet
    end
    
if stareftet_only
%     if ~ismember(tet,task{d}{ep}.STA_reftet(1,3:end))  % CA1 or higher
%         disp('not STA ref tet, skipping')
%         out.dummy = [];
%         return
%     end
    if ~ismember(tet,task{d}{ep}.STA_reftet(:,7))  % EC
        disp('not EC tet, skipping')
        out.dummy = [];
        return
    end
end
    
e_amp = eeg{index(1)}{index(2)}{index(3)}.data;
    srate = round(eeg{index(1)}{index(2)}{index(3)}.samprate);
timevec = geteegtimes(eeg{index(1)}{index(2)}{index(3)});
data_length = length(e_amp);

if selection_phasetet > 0
    % load phase eeg
        % first identify reference retrode
    regnum = selection_phasetet;
    phasetet = task{d}{ep}.STA_reftet(1,regnum);
        % load eeg data
    edata = loadeegstruct(animaldir,animalprefix,'eeg',d,ep,phasetet);
        e_phs = edata{d}{ep}{phasetet}.data;
        timevec2 = geteegtimes(edata{d}{ep}{phasetet});
    if length(timevec) ~= length(timevec2)   ||   ~all(timevec == timevec2)
        % interpolate to match the trace ends
        e_phs = interp1(timevec2,e_phs,timevec,'nearest');
    end
    clear timevec2 edata 
else
    phasetet = tet;
    e_phs = e_amp;
end


clear eeg


%% Filter and Hilbert transform
disp('CPU filtering')
tic
AmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);

for ii=1:length(AmpFreqVector)
    Af1 = AmpFreqVector(ii);
    Af2=Af1+AmpFreq_BandWidth;
    AmpFreq=eegfilt(e_amp',srate,Af1,Af2); % just filtering
    AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
end

for jj=1:length(PhaseFreqVector)
    Pf1 = PhaseFreqVector(jj);
    Pf2 = Pf1 + PhaseFreq_BandWidth;
    PhaseFreq=eegfilt(e_phs',srate,Pf1,Pf2); % this is just filtering 
    PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
end
toc



%% Iterate over states, first getting included periods then calculating MI

excludeperiods = cell(1,length(statespec));
duration = nan(1,length(statespec));
COMOD = cell(1,length(statespec));

for state = 1:length(statespec)
    
    COMOD{state} = single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
    
    % TIMEFILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timefilter = statespec{state};
    for ttt = 1:length(timefilter)
        output = feval(timefilter{ttt}{1}, animaldir, animalprefix, [d ep], timefilter{ttt}{3:end});
        output2 = evaluatefilter2(output, timefilter{ttt}{2});
        tmpexcludeperiods = getExcludePeriods(output2{d}{ep}(:,1),output2{d}{ep}(:,2)) ;
        if isempty(excludeperiods)
            excludeperiods{state} = tmpexcludeperiods;
        else
            excludeperiods{state} = combineExcludePeriods(excludeperiods{state}, tmpexcludeperiods);
        end
        clear tmpexcludeperiods
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % apply excludetimes to filtered eeg
    includedindices = ~isExcluded(timevec,excludeperiods{state});
    duration(state) = sum(includedindices)/srate;
    
    disp([num2str(d) ' ' num2str(ep) ' ' num2str(tet) ' state ' num2str(state) ' : ' num2str(duration(state))])
    
    % time-filter the filtered data
    AmpFreqTransformed_filtered = AmpFreqTransformed(:,includedindices);
    PhaseFreqTransformed_filtered = PhaseFreqTransformed(:,includedindices);
    
    % Do comodulation calculation
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
            [MI,MeanAmp]=ModIndex_v2(PhaseFreqTransformed_filtered(ii, :), AmpFreqTransformed_filtered(jj, :), position) ;
            COMOD{state}(counter1,counter2) = MI;
        end
    end
end

% plot if checking
if 0
    state_toplot = 1;
    figure
    contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,COMOD{state_toplot}',30,'lines','none')
    caxis([0 .001])
    set(gca,'fontsize',14)
    ylabel('Amplitude Frequency (Hz)')
    xlabel('Phase Frequency (Hz)')
    title([num2str(d) ' ' num2str(ep) ' ' num2str(tet)])
    colorbar
end

A.date = date;
A.index = index ;
A.phasetet = phasetet;
A.statespec = statespec;
A.statespec_descript = statespec_descript;
A.duration = duration ;
A.COMOD = COMOD;
A.frequency_amplitude = AmpFreqVector ;
A.frequency_phase = PhaseFreqVector ;
A.phasefreq_bandwidth = PhaseFreq_BandWidth ;
A.ampfreq_bandwidth = AmpFreq_BandWidth ;

if concat_flag
    out = [out  A];
else
    out = A;
end

% save output in the MIDATA directory
cd([animaldir 'MIDATA_new/']);
filestring = sprintf('%sMI-%d-%d-%d',animalprefix,d,ep,tet);
save(filestring,'out','-v7.3')

% don't actually return anything to the filter framework..
clear out
out.dummy = [];

end