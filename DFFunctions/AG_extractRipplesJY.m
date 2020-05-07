function AG_extractRipplesJY(animaldir, animalprefix, eventname, eventeegname, days, tetfilters, min_suprathresh_duration, nstd, varargin)
%function extractconsensus(directoryname, fileprefix, eventname, day, tetrode, min_suprathresh_duration, nstd, varargin)

% based on jai's version of kk_extractconsensus2
% improves upon older extractconsensus in these two ways:
%
%  % 1. Unlike before, normalise and smooths individual tetrodes eeg --before-- taking the mean across sites,
%           instead of smoothing after summing
%           Per jai, the median is used instead of the mean
%
%  % 2. instead of an arbitrary threshold, uses immobility times during all run sessions of that day to estimate gaussian noise of the consensus ripple band
%           and let that determine an appropriate threshold
%
%animaldir - example '/data99/user/animaldatafolder/', a folder
%                containing processed matlab data for the animal
%
%fileprefix	- folder name where the day's data is stored
%day		- the day(s) to process
%min_suprathresh_duration - the time (in sec) which the signal
%       must remain above threshold to be counted as as ripple; default 0.015 s.

%nstd		- the number of standard dev that ripple must be from mean to
%			be detextractconsensusected. Start with 2.
%           - if empty, nstd is calculated using noise estimation from data
%
%options	'stdev', stdev   sets the size of the standard deviation used to
%				allow direct comparison across epochs
%       	'baseline', b   sets the size of the baseline used to
%				allow direct comparison across epochs
%           'maxpeakval, m	- ripples with maximal peaks above this value
%				are exluded.  Use this avoid detecting noise
%				events. Default 1000
%           'calcbaselineestimate', calulcates and saves estimate of baseline (changes dramatically with velocity)
%           'immobripsonly', only output rip events that occur during immobility (same cutoff as for noise estimate)
%
% Outputs:
%ripples 	- structue with various fields, including the following which
%			describe each ripple.
%	starttime - time of beginning of ripple
%	endtime	  - time of end of ripple
%	maxthresh - the largest threshold in stdev units at which this ripple
%			would still be detected.

disp(['processing ' animalprefix])

stdev = 0;
baseline = 0;
savetrace = 0;
maxpeakval = 2000;
% define the standard deviation for the Gaussian smoother which we
% apply before thresholding (this reduces sensitivity to spurious
% flucutations in the ripple envelope)
smoothing_width = 0.004; % 4 ms
%smoothing_width = 0.008; % 8 ms
numTF = length(tetfilters);
calcbaselineestimate = 0;
immobripsonly = 0;
epfilter = 'isequal($environment,''run'')';
immobilefilter = {'kk_get2dstate', '($immobility == 1)','immobility_velocity',4,'immobility_buffer',0};

% process varargin and overwrite default values 
if (~isempty(varargin))
    assign(varargin{:});
end

tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo');
task = loaddatastruct(animaldir, animalprefix, 'task');

runeps_all = evaluatefilter(task,epfilter);   


for d = days
    
    try  %check to see if a ca1ripplescons already has been created for this day
        ca1ripplesconsJY = loaddatastruct(animaldir, animalprefix, 'ca1ripplesconsJY',d);
        disp('ca1ripplesconsJY detected; appending/overwriting')
    catch
        disp('no ca1ripplesconsJY detected; creating')
    end
    
    
    % ---part 1: estimate baseline
    
    baseline_TF = nan(1,numTF);
    stdev_TF = nan(1,numTF);
    thresh_TF = nan(1,numTF);
%     if calcbaselineestimate
%         noise_baseline_TF = nan(1,numTF);
%         noise_stdev_TF= nan(1,numTF);
%         noise_thresh_TF = nan(1,numTF);
%     end
    runeps = unique(runeps_all(runeps_all(:,1) == d,2))';
    % load eeg data
    
    
    for TF = 1:length(tetfilters)       
        tetfilter = tetfilters{TF};
        tetlist =  evaluatefilter(tetinfo{d},tetfilter);
        tetlist = unique(tetlist(ismember(tetlist(:,1),runeps),2));
        if isempty(tetlist)
            disp(sprintf('no valid detecting tetrodes, d %d',d))
            continue
        end
        eventeeg = loadeegstruct(animaldir, animalprefix, eventeegname, d, runeps, tetlist');
        powervals = [];
        noisevals=[];
        noisevals_all=cell(max(runeps),1);
        
        for ep = runeps
            
            % use times vector from one of the tetrodes as a reference
            eegtimesvec_ref{ep} = geteegtimes(eventeeg{d}{ep}{tetlist(1)});
            epoch_starttime{ep} = eegtimesvec_ref{ep}(1);
            epoch_endtime{ep} = eegtimesvec_ref{ep}(end);
            %numsamples_refeeg = length(eventeeg{d}{ep}{tetlist(1)}.data(:,1));
            %if numsamples_refeeg ~= length(eegtimesvec_ref)
            %    keyboard
            %end            
            amptet = nan(length(tetlist),length(eegtimesvec_ref{ep}));  % initialize
            
            for tt = 1:length(tetlist)
                tet = tetlist(tt);
                % take the amplitude
                % data = double(eventeeg{d}{ep}{tet}.data(:,1))';
                % take the envelope
                data = double(eventeeg{d}{ep}{tet}.data(:,3))';
                
                % if there are brief periods of wild noise (filtered amplitude exceeds maxval),set values to 0 as a patch fix and report
                wildnoise = (abs(data) > maxpeakval);
                if any(wildnoise)
                    ms_of_wildnoise = round(1000*sum(wildnoise)/1500);
                    if ms_of_wildnoise > 5000
                        % if more than 5 s of wild noise, don't report consensus events in this case -- abandon the epoch
                        disp(sprintf('anim %s day%dep%d tet %d: wild noise periods totalling %d millisec (if 1500 Hz) -- (!!) exceeds 5 s',...
                            animalprefix,d,ep,tet,ms_of_wildnoise))
                        disp('(!!!!!) excluding epoch because of this -- consider eliminating tet as detector!')
                        keyboard
                        continue
                    else
                        disp(sprintf('anim %s day%dep%d tet %d: wild noise periods totalling %d millisec (if 1500 Hz) excluded in threshold calc',...
                            animalprefix,d,ep,tet,ms_of_wildnoise))
                        % simply set violating values to zero
                        data(wildnoise) = 0;
                    end
                end
                
                
                %if length(data) == numsamples_refeeg
                curr_amptet= data;
                %else
                % if mismatch in # of samples for this tetrode (different DSPs often start at slightly different times),
                % 'interpolate' (really just quick way of appending NaNs on the deficient end)
                %    eegtimesvec = geteegtimes(eventeeg{d}{ep}{tet});
                %    interpdata = interp1(eegtimesvec,data,eegtimesvec_ref,'nearest')';
                %   curr_amptet = interpdata;
                %end

                % smooth envelope
                Fs = eventeeg{d}{ep}{tet}.samprate;
                kernel = gaussian(smoothing_width*Fs, ceil(8*smoothing_width*Fs));
                curr_amptet = (smoothvect(curr_amptet, kernel));
                
                %normalalize each tet
                curr_std=nanstd(curr_amptet);
                curr_mean=nanmean(curr_amptet);
                curr_amptet=(curr_amptet-curr_mean)./curr_std;
                
                amptet(tt,:)=curr_amptet;   
            end
            
            % take the median to avoid one tetrode dominating power and save this for later
            powertrace{ep} = nanmedian(amptet,1);
            %powertrace_mean=mean(amptet,1);
            
            %get immobility periods and save these for later
            immoperiodstmp = kk_evaluatetimefilter(animaldir,animalprefix, {immobilefilter}, [d ep]);
            immoperiods{ep}=immoperiodstmp{d}{ep};
            immovec{ep} = logical(list2vec(immoperiods{ep},eegtimesvec_ref{ep}));
            validinds = immovec{ep} & ~isnan(powertrace{ep}');            
            powervals = [ powervals   powertrace{ep}(validinds) ];
            
%             if calcbaselineestimate
%                 % estimate continuous baseline level & remove large deviations
%                 powertrace_z=zscore(powertrace);
%                 powertrace_filt=powertrace_z<2;
%                 load('Hz0_5.mat');
%                 % filter with low pass filter to get baseline variation
%                 datafilt05 = filtfilt(filter_05, 1, powertrace(powertrace_filt));
%                 % interploate back to all times by filling the gaps where large deviations use to be
%                 powertrace_fill=interp1(eegtimesvec_ref(powertrace_filt),datafilt05,eegtimesvec_ref);
%                 noisevals=[noisevals powertrace_fill(validinds)];
%                 noisevals_all{ep}=powertrace_fill; 
%                 % calculate mean and sd of baseline estimation; noise threshold is 2SD above mean
%                 noise_baseline_TF(TF) = mean(noisevals);
%                 noise_stdev_TF(TF) = std(noisevals);
%                 noise_thresh_TF(TF) = noise_baseline_TF(TF) + 2* noise_stdev_TF(TF);
%             end
            
            totalbaselinedur = round( length(powervals) / Fs);
            disp(sprintf('TF #%d: %d s of baseline data',TF,totalbaselinedur));
                    
            % calculate mean and std of powertrace, to use in extractevents
            baseline_TF(TF) = mean(powervals);
            stdev_TF(TF) = std(powervals);

            % get threshold by estimating noise spread
            [thresh_TF(TF),riphistpdf(TF,:)] = jy_variableripthreshold_corecalculation({powervals},'animal_dir',animaldir,'animal_name',animalprefix,'day',d,'target_epochs',runeps);
        end
 
    % Now, calculate the actual power trace + consensus events
    
    clear eventscons ev
    clear eventtrace tr
    
    for ep = runeps

            mindur = round(min_suprathresh_duration * Fs);
            %powertrace = sum(amptet,1);
            %powertrace = median(amptet,1);

                baseline = baseline_TF(TF);
                stdev = stdev_TF(TF);
                thresh = thresh_TF(TF);
                
                %noise_baseline = noise_baseline_TF(TF);
                %noise_stdev = noise_stdev_TF(TF);
                %noise_thresh = noise_thresh_TF(TF);                
                %powertrace_fill_z=(powertrace_fill-noise_baseline)./noise_stdev;
                
            % extract the events if this is a valid trace
            if (thresh >-5) && any(find(powertrace{ep} < baseline))
                disp(sprintf('%s day %d epoch %d ',animalprefix,d,ep))
                
                % some odd error within the extractevents .mex makes it lock up sometimes... 
                % to avoid, duplicate the powertrace and then take the first half of the output
                % also need to make sure all values are positive so add an offset factor
                offsetfactor=10;
                tmpevent = extractevents([powertrace{ep}'; powertrace{ep}']+offsetfactor, thresh+offsetfactor, baseline+offsetfactor, 0, mindur, 0)';
                if mod(size(tmpevent,1),2)==0
                    lastind = size(tmpevent,1)/2;
                else
                    lastind = (size(tmpevent,1)-1)/2;
                end
                tmpevent=tmpevent(1:lastind,:);
                
                % only save events that occur during immobility
                startind = tmpevent(:,1);
                endind = tmpevent(:,2);
                immoevents = isExcluded(epoch_starttime{ep} + startind / Fs, immoperiods{ep}) & isExcluded(epoch_starttime{ep} + endind / Fs, immoperiods{ep});
                
                % Event output
                ev.eventname = eventname;
                ev.nstd = nstd;
                ev.min_suprathresh_duration = min_suprathresh_duration;
                ev.tetfilter = tetfilter;
                ev.tetlist = tetlist;
                ev.starttime = epoch_starttime{ep} + startind(immoevents) / Fs;
                ev.endtime = epoch_starttime{ep} + endind(immoevents) / Fs;
                ev.threshold=thresh;
                ev.maxthresh = (tmpevent(immoevents,9) - baseline-offsetfactor) / stdev;
                ev.baseline = baseline;
                ev.stdev = stdev;
                ev.powertrace = powertrace{ep};
                ev.riphistpdf = riphistpdf;
                %alternative maxthresh recalculated as multiple of thresh, not baseline+sd
                ev.maxthresh_abovenoise = (ev.maxthresh*stdev + baseline)/thresh;
                if 0
                    figure; hold on
                    otherrips=[ev.starttime ev.endtime]-epoch_starttime{ep};
                    ttt=otherrips;
                    nspikes   = length(ttt); % number of elemebts / spikes
                    for nn = 1:nspikes % for every spike
                        jbfill([(ttt(nn,1)) (ttt(nn,2))],...
                            [10 10],[-2 -2],...
                            [250 220 150]./255,[],[],1);
                        %line(([tt(nn) tt(nn)+t(nn)]-mean(tstartdata))/10000,[0 length(spikes)+1],'Color','b','LineWidth',0.1);
                    end
                    plot(eegtimesvec_ref{ep}-eegtimesvec_ref{ep}(1),powertrace{ep},'k')
                    hold on
                    %plot(eegtimesvec_ref -eegtimesvec_ref(1),noisevals_all{ep},'r')
                    % line([eegtimesvec_ref(1) eegtimesvec_ref(end)]-eegtimesvec_ref(1),[noise_thresh noise_thresh],'color','g')
                     line([eegtimesvec_ref{ep}(1) eegtimesvec_ref{ep}(end)]-eegtimesvec_ref{ep}(1),[thresh thresh],'color',[0.5 0.5 0.5])
                     line([eegtimesvec_ref{ep}(1) eegtimesvec_ref{ep}(end)]-eegtimesvec_ref{ep}(1),...
                         [baseline baseline],'color','c')
                   %line([eegtimesvec_ref(1) eegtimesvec_ref(end)]-eegtimesvec_ref(1),...
                        % [baseline baseline]-nstd*stdev,'color','c')
                     line([eegtimesvec_ref{ep}(1) eegtimesvec_ref{ep}(end)]-eegtimesvec_ref{ep}(1),...
                         [baseline baseline]+2*stdev,'color',[0.5 0.5 0.5],'linestyle','--')
                        line([eegtimesvec_ref{ep}(1) eegtimesvec_ref{ep}(end)]-eegtimesvec_ref{ep}(1),...
                         [baseline baseline]+3*stdev,'color',[0.5 0.5 0.5],'linestyle',':')
                      plot(eegtimesvec_ref{ep}(immovec{ep})-eegtimesvec_ref{ep}(1),-1*ones(sum(immovec{ep}),1),'o')%,'color',[1 0 0],'o'
                end

            else
                fprintf('%s day %d epoch %d no valid data for extraction\n',animalprefix,d,ep)
                % Event output
                ev.eventname = '';
                ev.nstd = [];
                ev.min_suprathresh_duration = [];
                ev.tetfilter = '';
                ev.tetlist = [];
                ev.starttime = [];
                ev.endtime = [];
                ev.maxthresh = [];
                ev.mean_runimmo = [];
                ev.stdev_runimmo = [];
                ev.threshold=[];
            end
            ev.timerange = [epoch_starttime{ep} epoch_endtime{ep}];
            ev.samprate = Fs;
            ev.baseline = baseline;
            ev.std = stdev;

            ca1ripplesconsJY{d}{ep}{TF} = ev;
            clear ev;
   
        end 
    end

    save(sprintf('%s%s%sconsJY%02d.mat', animaldir, animalprefix, eventname,d),'ca1ripplesconsJY');
    clear ripplescons immoperiods immovec eegtimesvec_ref;
    disp(sprintf('consensus saved for %s %s, day %d',animalprefix,eventname,d))   
end

end
