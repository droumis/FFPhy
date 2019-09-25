function extractsleep(animaldir, animalprefix, tetfilter, epochfilter, varargin)

%%% OBSOLETE

%	Reads in the ripple files from the specified day and tetrodes and
%	extracts all of the ripples from that tetrodes.

%   Reads in these files:
% delta
% theta
% position
% gamma (20-80 Hz)
% spindle

%animaldir - example '/data99/user/animaldatafolder/', a folder
%                containing processed matlab data for the animal
%
%fileprefix	- folder name where the day's data is stored
%

%day		- the day to process
%
%tetfilter	- the tetrode(s) to use to detect sleep (i.e. cortical tets)
%
%min_rem_duration
%		- the time (in seconds) which the signal
%       must remain above ratio to be counted as a rem period
%




% Outputs:
% sws 	- structure w/ various fields:
%	starttime - time of beginning of sws period
%	endtime	  - time of end of sws period
%	maxthresh - the largest threshold in stdev units at which this sws
%                   detected
% rem 	- structure w/ various fields:
%	starttime - time of beginning of sws period
%	endtime	  - time of end of sws period
%	maxthresh - the largest threshold in stdev units at which this sws
%                   detected

ctx_flag = 1;

smoothing_width = 0.1;
maxpeakval = 1000;
samethreshperday = 0;

nrem_flag = 0;
rem_flag = 0;
spindle_flag = 0;
sws_flag = 0;
bimodal_study = 0;
immobility = [10 1]; 


rem_thresh = 1.00;     % find using bimodality study
nrem_thresh = 0.60;
spindle_thresh = 1;  % following destexhe, SDs from baseline magnitude

mindur_spindle = 0.5;
mindur_rem = 0;
mindur_nrem = 0;



for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'sleep'
                sleep_flag = varargin{option+1};
            case 'rem'
                rem_flag = varargin{option+1};
            case 'nrem'
                nrem_flag = varargin{option+1};
            case 'spindle'
                spindle_flag = varargin{option+1};
            case 'sws'
                sws_flag = varargin{option+1};
            case 'cortex'                   % this gets MUA and rmsgamma from all ctx tetrodes first (takes a few min..)
                ctx_flag = varargin{option+1};
            case 'immobility'      % time transpired at after which sleep is counted
                time_immobile = varargin{option+1}(1);    % in seconds
                velocity_thresh = varargin{option+1}(2); 
            case 'smoothing_width'      % time transpired at after which sleep is counted
                smoothing_width = varargin{option+1};    % in seconds
            case 'maxpeakval'
                maxpeakval = varargin{option+1};
            case 'samethreshperday'
                samethreshperday = varargin{option+1};
            case 'bimodal_study'
                bimodal_study = varargin{option+1};
            case 'rem_thresh'
                rem_thresh = varargin{option+1};     % thetaratio over which rem is taken
            case 'nrem_thresh'
                nrem_thresh = varargin{option+1};    % thetaratio under which nrem is taken
            case 'spindle_thresh'
                spindle_threhs = varargin{option+1};
            case 'sws_thresh'
                sws_thresh = varargin{option+1};    % thetaratio under which nrem is taken
            case 'mindur_rem'
                mindur_rem = varargin{option+1};                
            case 'mindur_nrem'
                mindur_nrem = varargin{option+1};
            case 'mindur_sws'
                mindur_sws = varargin{option+1};
            case 'mindur_spindle'
                mindur_spindle = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end


% for bimodality study
if bimodal_study
figure
plot_counter = 1;
end
%





tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo');
task = loaddatastruct(animaldir,animalprefix,'task');
    dayep = evaluatefilter(task,epochfilter);
    days = unique(dayep(:,1))';


    

    
for day = days;
    
    % retrieve pos
    pos = loaddatastruct(animaldir, animalprefix, 'pos', day);
    
    % retrieve mua
    multi = loaddatastruct(animaldir, animalprefix, 'multi', day);
    
    % identify sleep epochs for this day
    epochs = dayep(find(dayep(:,1)==day),2)';
    
    for ep = epochs
        
        % First infer sleep periods via immobility parameters
        velocity = pos{day}{ep}.data(:,9);
        postimes = pos{day}{ep}.data(:,1);
            % generate [start end] list of immobile epochs
        immobile = vec2list(velocity < velocity_thresh,postimes);
            % truncate immobility periods in front by time_immobile to get
                % the inferred sleep periods
        immobile = [immobile(:,1)+time_immobile, immobile(:,2)];
        sleepperiods = immobile(find( (immobile(:,2)-immobile(:,1)) > 0),:);
            % next, grab tetrodes
        tetlist =  evaluatefilter(tetinfo{day}{ep},tetfilter)';        
        
        % move to the EEG directory
        cd([animaldir,'/EEG']);
        if ~isempty(tetfilter)
            tet = 1:1000;
        end
        
        
        % If ctx_flag specified, collect from all the CTX tetrodes:
               % i. MUA activity
               % ii. rmsgamma   (after Mukovski--Volgushev-2007)
        if ctx_flag
            
            % mua variables
            muakernel_length = .010;          % in sec
                % outputs
            ctxmua = [];
            muafiring = {};
            muafiring_allctx = [];
            
            % rmsgamma variables
            rmsframe = .005;
            smooth_duration = .010;       % in sec
                % outputs
            ctxrms_smooth = [];
            rmsgamma_smooth = {};
            rmsgamma_smooth_allctx = [];
            
            % general variables
            numctxtet = 0;
            dummyflag = 0;  
            

                
            for tet = 1:length(tetlist)
                if isfield(tetinfo{day}{ep}{tet},'area')
                    
                    % retrieve an eeg times vector to use as a reference
                    % for cortical mua/eeg -- doesn't matter which tetrode
                    if ~dummyflag
                        tmp1 = dir(sprintf('%seeg%02d-%d-%02d.mat',animalprefix, day,ep, tet));
                        tmp2 = dir(sprintf('%sgamma%02d-%d-%02d.mat',animalprefix, day,ep, tet));
                        load(tmp1.name);
                        load(tmp2.name);
                        ctxtimes = geteegtimes(eeg{day}{ep}{tet})';
                        Fs_eeg = floor(eeg{day}{ep}{tet}.samprate);
                        Fs_gamma = floor(gamma{day}{ep}{tet}.samprate);
                        ctxreftet = tet;
                        % error check the sampling rates
                        if (Fs_eeg ~= Fs_gamma)
                            error('sampling rates for filtered eegs need to match')
                        end
                        dummyflag = 1;
                    end
                    
                    % get mua
                    mua = multi{day}{ep}{tet}/10000;
                    if strcmp(tetinfo{day}{ep}{tet}.area,'ctx')
                        numctxtet = numctxtet + 1;
                        ctxmua = [ctxmua ; mua];
                    end
                        % bin the spikes
                    N = hist(mua,ctxtimes);
                        % smooth the rate
                    muakernel = gaussian(muakernel_length*Fs_eeg,8*muakernel_length*Fs_eeg);
                    muafiring{tet} = smoothvect(N,muakernel)' * Fs_eeg;
                    
                    % get rmsgamma
                    tmp2 = dir(sprintf('%sgamma%02d-%d-%02d.mat',animalprefix, day,ep, tet));
                    load(tmp2.name);
                    % check for different eeg length (different dsps)
                        % pad NaN if different length to be able to take
                        % mean later
                    
                    
                    if length(gamma{day}{ep}{tet}.data(:,1)) == length(ctxtimes)
                        gammadata = double(gamma{day}{ep}{tet}.data(:,1));
                    % if mismatch, pad difference w/ NaNs
                    else
                        tmp_gamma = [];
                        tmp_times = geteegtimes(gamma{day}{ep}{tet})';
                        % if ctxref eeg start is earlier
                            if ctxtimes(1) < tmp_times(1)
                                ind = lookup(tmp_times(1),ctxtimes);
                                tmp_gamma(1:(ind-1)) = NaN;
                                % and then if ctxref ends earlier
                                if ctxtimes(end) < tmp_times(end)
                                    midind = lookup(ctxtimes(end),tmp_times);
                                    tmp_gamma(ind:(ind+midind-1)) = gamma{day}{ep}{tet}.data(1:midind,1);
                                % and then if ctxref ends later
                                else
                                    tmp_gamma(ind:(ind+length(tmp_times)-1)) = gamma{day}{ep}{tet}.data(1:midind,1);
                                    % pad w/ NaN
                                    padlength = length(ctxtimes)-length(tmp_gamma);
                                    tmp_gamma = [ tmp_gamma  nan(1,padlength) ];
                                end
                            % if ctxref eeg start is later
                            else
                                startind = lookup(ctxtimes(1),tmp_times);
                                % and then if ctxref eeg ends earlier
                                if ctxtimes(end) < tmp_times(end)
                                    tmp_gamma = gamma{day}{ep}{tet}.data(startind:(startind + length(ctxtimes)-1),1);
                                % and then if ctxref eeg ends later
                                else
                                    endind = lookup(tmp_times(end),ctxrefeeg);
                                    tmp_gamma(1:endind) = gamma{day}{ep}{tet}.data(startind:(startind+endind-1),1);
                                    % pad w/ NaN
                                    tmp_gamma((endind+1):length(ctxtimes)) = NaN;
                                end
                            end
                        gammadata = double(tmp_gamma);
                    end
                    
                    % a. following Mukovski, calculate rms in 5 ms window of filtered eeg   
                    rmsgamma = smooth(gammadata,floor(rmsframe*Fs_eeg));
                    % b. next smooth the signal
                    rmsgamma_smooth{tet} = smooth(rmsgamma,floor(smooth_duration * Fs_eeg),'mean');
                    ctxrms_smooth = [ctxrms_smooth ; rmsgamma_smooth{tet}'];
                end
            end
            
        % combine mua from all ctx tetrodes
        N_all = hist(ctxmua,ctxtimes);
        muafiring_allctx = ( smoothvect(N_all,muakernel)' * Fs_eeg ) / numctxtet;
        
        % combine rmsgamma from all ctx tetrodes
        rmsgamma_smooth_allctx = mean(ctxrms_smooth,1);      
        
        end
        
        

        
        for tet = tetlist
            
            tmpflist1 = dir(sprintf('%stheta%02d-%d-%02d.mat', animalprefix, day,ep, tet));
            tmpflist2 = dir(sprintf('*delta%02d-%d-%02d.mat', day,ep, tet));
            tmpflist3 = dir(sprintf('*supratheta%02d-%d-%02d.mat', day,ep, tet));
            if sws_flag
                % 20-80 Hz
                tmpflist4 = dir(sprintf('%sgamma%02d-%d-%02d.mat',animalprefix, day,ep, tet));
                tmpflist5 = dir(sprintf('%seeg%02d-%d-%02d.mat',animalprefix, day,ep, tet));
            end
            if spindle_flag
            tmpflist6 = dir(sprintf('%sspindle%02d-%d-%02d.mat',animalprefix, day,ep, tet));
            end            

            
            % iterate through all epochs for this day and tet
                
                % load the eeg files
                load(tmpflist1.name);
                load(tmpflist2.name);
                load(tmpflist3.name);
                if sws_flag
                   load(tmpflist4.name);
                   load(tmpflist5.name); 
                end
                if spindle_flag
                    load(tmpflist6.name); 
                end
                
                % convert the envelope fields to double
                tenv = double(theta{day}{ep}{tet}.data(:,3));
                denv = double(delta{day}{ep}{tet}.data(:,3));
                stenv = double(supratheta{day}{ep}{tet}.data(:,3));
                
                % check that samprates are equivalent
                fs1 = theta{day}{ep}{tet}.samprate;
                fs2 = delta{day}{ep}{tet}.samprate;
                fs3 = supratheta{day}{ep}{tet}.samprate;
                if (fs1 ~= fs2) || (fs1 ~= fs3) || (fs2 ~= fs3)
                    error('sampling rates for filtered eegs need to match')
                end
                    
                % smooth the magnitude envelope
                samprate = fs1;
                kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
                tenv = smoothvect(tenv, kernel);
                denv = smoothvect(denv, kernel);
                stenv = smoothvect(stenv, kernel);
                
                % result
                thetaratio = tenv ./ (denv + stenv);
                deltaratio = denv ./ (tenv);
                times_filteeg = geteegtimes(theta{day}{ep}{tet})';
                

              % (Optional: study bimodality of theta/(delta + supratheta) ratio)
                if bimodal_study
                    perioddur = 1;
                    period_ind = floor(perioddur * fs1);
                    numwindows = floor((length(thetaratio)/fs1)/perioddur);
                    powers = [];
                    for p = 1:numwindows
                        powers = [powers ; mean(thetaratio(1+(p-1)*period_ind:p*period_ind))];
                    end
                    hold on
                    randcolor = [rand rand rand];
                    [N X] = hist(powers,0:.05:3);
                    if plot_counter > 60
                        figure
                        hold on
                        plot_counter = 1;
                    end
                    subplot(6,10,plot_counter)
                    h = bar(X,N);
                    xlim([0 3])
                    set(h(1),'facecolor',randcolor)
                    set(h(1),'edgecolor',randcolor)
                    patch=findobj(h,'Type','patch');
                    set(patch,'facealpha',.4)
                    set(patch,'edgealpha',.4)
                    titlestring=sprintf('%d %d %d',[day e tet]);
                    title(titlestring)
                    plot_counter = plot_counter + 1;
                    continue
                end
                
                % Find REM vs NREM periods within the inferred sleep periods.

                % initialize
                sleepvec = zeros(size(times_filteeg));
                
                % transform sleepperiods list into 0-1 vector
                for p=1:size(sleepperiods,1)
                    startind = lookup(sleepperiods(p,1),times_filteeg);
                    endind = lookup(sleepperiods(p,2),times_filteeg);
                    sleepvec(startind:endind) = 1;
                end
                
                % take intersection of sleep times w/ high (REM) vs. low (NREM) theta times
                rem_vec = sleepvec & (thetaratio >= rem_thresh) ;
                nrem_vec = sleepvec & (thetaratio < nrem_thresh) ;
                sws_vec = sleepvec & (deltaratio > sws_thresh) ;
            
                % transform into [start end] format
                remlist = vec2list(rem_vec,times_filteeg);
                nremlist = vec2list(nrem_vec,times_filteeg);
                swslist = vec2list(sws_vec,times_filteeg);
                
                % eliminate periods that are shorter than minimum durations
                if mindur_rem > 0
                    remlist(find((remlist(:,2)-remlist(:,1)) < mindur_rem),:) = [];
                end
                if mindur_nrem > 0
                    nremlist(find((nremlist(:,2)-nremlist(:,1)) < mindur_nrem),:) = [];
                end
                if mindur_sws > 0
                    swslist(find((swslist(:,2)-swslist(:,1)) < mindur_nrem),:) = [];
                end                
                
                %% Spindle analysis (after Destexhe--) %%%
                if spindle_flag
                    % retrieve spindle filtered eeg
                    spenv = double(spindle{day}{ep}{tet}.data(:,1));
                    % obtain "rms envelope" :
                    spbaseline = mean(spenv);
                    sprms = sqrt((spenv - spbaseline).^2);
                        sprms_std = std(sprms);
                        sprms_mean = mean(sprms);
                    % smooth w/ 100 ms Gaussian window"
                    spindlekernel = gaussian(fs1 * 0.1, 8 * fs1 * 0.1);
                    sprms_smooth = smoothvect(sprms,spindlekernel);
                    % identify where above spindle threshold (try 1 SD, after Destexhe) for at least minimum duration
                    spindle_vec = sprms_smooth-sprms_mean > spindle_thresh * sprms_std;
                    % take intersection with nrem
                    spindle_vec = spindle_vec & nrem_vec;
                    % convert to list
                    spindlelist = vec2list(spindle_vec,geteegtimes(spindle{day}{ep}{tet}));
                    % remove violations of minimum duration
                    if mindur_spindle > 0
                        spindlelist(find((spindlelist(:,2)-spindlelist(:,1)) < mindur_spindle),:) = [];
                    end
                    % lastly, iterate through each spindle and retrieve mean
                    % magnitude for each
                    spindle_mean_rms_zscore = [];
                    sptimes = geteegtimes(spindle{day}{ep}{tet});
                   for s = 1:size(spindlelist,1)
                       startind = lookup(spindlelist(s,1),sptimes);
                       endind = lookup(spindlelist(s,2),sptimes);
                       spindle_mean_rms_zscore(s,:) = (mean(sprms(startind:endind))-sprms_mean) / sprms_std;
                   end
                 end
                %%

                
            
           
                
                
                %% SWS analysis for UDS  
                      % NREM with low enough power threshold should be composed of SWS.

                if sws_flag

                if isfield(tetinfo{day}{ep}{tet},'area') && strcmp(tetinfo{day}{ep}{tet}.area,'ctx')
   
                    % set how long rms frame should be
                    frame = .005;  
                        frame_numindices = floor(frame * gamma{day}{ep}{tet}.samprate);
                    Fs = gamma{day}{ep}{tet}.samprate;
                    
                    % if spindles specified, then remove spindles from sws periods
                    if spindle_flag
                        sws_vec = sws_vec & ~spindle_vec;
                        swslist = vec2list(sws_vec,times_filteeg);
                    end
                   
                    % remove violations of minimum duration
                    if mindur_sws > 0
                        swslist(find((swslist(:,2)-swslist(:,1)) < mindur_sws),:) = [];
                    end
                    
                    % note that since filtered gamma eeg is sampled at higher Fs
                        % we thus need to reconstitute swslist into gamma times vector    
                    sws_vec2 = list2vec(swslist,geteegtimes(gamma{day}{ep}{tet}));
                        sws_vec2(sws_vec2==0) = NaN;
                        % a. following Mukovski, calculate rms in 5 ms window of filtered eeg                    
                    rmsgamma = smooth(double(gamma{day}{ep}{tet}.data(:,1)),frame_numindices);
                        % b. next smooth the signal
                    smooth_duration = .050;  % in sec
                        smooth_duration_ind = floor(smooth_duration * gamma{day}{ep}{tet}.samprate);
                    rmsgamma_smooth = smooth(rmsgamma,smooth_duration_ind,'mean');

                    % plot each sws period, centered in a window
                    
                    windowsize = 5 ;  % in sec
                    timevec = (1:(windowsize * 1500)) / 1500;
                    times_eeg = geteegtimes(eeg{day}{ep}{tet});
                    
                    if 0
                    for w = 1:size(swslist,1)
                        figure; hold on;
                        title(sprintf('%d %d %d',day,ep,tet));
                        midtime = (swslist(w,2)+swslist(w,1)) / 2;
                        starttime = midtime - windowsize / 2;
                        startind = lookup(starttime,times_eeg);
                        plotindices = startind:(startind+length(timevec)-1);
                        % eeg
                        plot(timevec,eeg{day}{ep}{tet}.data(plotindices),'Color',[.8 .8 .8]);              % all eeg
                            maxeeg = max(eeg{day}{ep}{tet}.data(plotindices));
                            mineeg = min(eeg{day}{ep}{tet}.data(plotindices));
                        plot(timevec,sws_vec2(plotindices).*eeg{day}{ep}{tet}.data(plotindices),'Color',[0 0 0]);                 % sws eeg in solid
                        % mua
                            maxmua = max(muafiring_allctx(plotindices));
                            minmua = min(muafiring_allctx(plotindices));
                            scalefactor1 = (maxeeg-mineeg)/(maxmua-minmua);
                        plot(timevec,scalefactor1 * muafiring_allctx(plotindices)  + mineeg,'Color',[255 192 203]/255) % all
                        plot(timevec,scalefactor1 * (muafiring_allctx(plotindices) .* sws_vec2(plotindices)) + mineeg,'Color',[255 105 180]/255) % solid                           
                        % rms
                            maxrg = max(rmsgamma(plotindices));
                            minrg = min(rmsgamma(plotindices));
                            scalefactor2 = (maxeeg-mineeg)/(maxrg-minrg);
                        plot(timevec,scalefactor2 * rmsgamma_smooth_allctx(plotindices) + 0.5 * mineeg,'Color',[.8 .8 1])              % all
                        plot(timevec,scalefactor2 * rmsgamma_smooth_allctx(plotindices)' .* sws_vec2(plotindices) + 0.5 * mineeg,'b')   % solid    
                        
                        pause
                        close all
                    end
                    end
                    
                    % look for bimodality in rmsgamma
                    signal = rmsgamma_smooth.*sws_vec2;
                        signal = signal(signal ~= 0);
                    figure
                    hist(signal,200)
                    
                    % look for bimodality in mua
                        % interpolate muafiring_allctx if sampled
                        % differently from current tetrode
                        muafiring_allctx2 = interp1(ctxtimes,muafiring_allctx,geteegtimes(gamma{day}{ep}{tet}))';
                    signal = muafiring_allctx2.*sws_vec2;
                    figure
                    hist(signal,500)
                    %%%%
                    
                    
                end
                end
                
                
                
                  
                
                
                % output rem periods if there are any
                if rem_flag
                    if ~isempty(remlist)
                        r.starttime = remlist(:,1);
                        r.endtime = remlist(:,2);
                    else
                        r.starttime = [];
                        r.endtime = [];
                    end
                    
                    r.timerange = [0 length(tenv)/fs1] + theta{day}{ep}{tet}.starttime;
                    r.samprate = theta{day}{ep}{tet}.samprate;
                    r.rem_thresh = rem_thresh;
                    r.baseline = 'relying on static threshold a la Buzsaki';
                    r.std = 'relying on static threshold a la Buzsaki';
                    r.mindur_rem = mindur_rem;
                    
                    % install output in proper place
                    rem{day}{ep}{tet} = r;
                    clear r;
                end
                           
                % output nrem periods if there are any
                if nrem_flag
                    if ~isempty(nremlist)
                        nr.starttime = nremlist(:,1);
                        nr.endtime = nremlist(:,2);
                    else
                        nr.starttime = [];
                        nr.endtime = [];
                    end
                    
                    nr.timerange = [0 length(tenv)/spindle{day}{ep}{tet}.samprate] + spindle{day}{ep}{tet}.starttime;
                    nr.samprate = spindle{day}{ep}{tet}.samprate;
                    nr.spindle_thresh = spindle_thresh;
                    nr.baseline = 'relying on absolute threshold a la Buzsaki';
                    nr.std = 'relying on absolute threshold a la Buzsaki';
                    nr.mindur_nrem = mindur_nrem;
                    
                    % install output in proper place
                    nrem{day}{ep}{tet} = nr;
                    clear nr;
                end
           
                % output spindle periods if there are any
                if spindle_flag
                    if ~isempty(spindlelist)
                        sp.starttime = spindlelist(:,1);
                        sp.endtime = spindlelist(:,2);
                        sp.spindle_mean_rms_zscore = spindle_mean_rms_zscore;
                    else
                        sp.starttime = [];
                        sp.endtime = [];
                        sp.spindle_mean_rms_zscore = [];
                    end
                    
                    sp.timerange = [0 length(spenv)/fs1] + spindle{day}{ep}{tet}.starttime;
                    sp.samprate = spindle{day}{ep}{tet}.samprate;
                    sp.spindle_rms_mean = sprms_mean;               % mean of rms (of filt eeg) over epoch
                    sp.spindle_rms_std = sprms_std;                % std of rms (of filt eeg) over epoch
                    sp.mindur_spindle = mindur_spindle;
                    sp.nrem_thresh = nrem_thresh;
                    
                    % install output in proper place
                    spindles{day}{ep}{tet} = sp;
                    clear sp;
                end
                
                
                % output nrem periods if there are any
                if sws_flag
                    if ~isempty(swslist)
                        sw.starttime = swslist(:,1);
                        sw.endtime = swslist(:,2);
                    else
                        sw.starttime = [];
                        sw.endtime = [];
                    end
                    
                    sw.timerange = [0 length(tenv)/fs1] + theta{day}{ep}{tet}.starttime;
                    sw.samprate = theta{day}{ep}{tet}.samprate;
                    sw.sws_thresh = sws_thresh;
                    sw.baseline = 'relying on absolute threshold a la Buzsaki';
                    sw.std = 'relying on absolute threshold a la Buzsaki';
                    sw.mindur_sws = mindur_sws;
                    
                    % install output in proper place
                    sws{day}{ep}{tet} = sw;
                    clear sw;
                end
                
                
        end
        
        % output sleep periods if there are any
        if sleep_flag
            if ~isempty(sleepperiods)
                slp.starttime = sleepperiods(:,1);
                slp.endtime = sleepperiods(:,2);
                slp.timerange = [pos{day}{ep}.data(1,1) pos{day}{ep}.data(end,1)];
            else
                slp.starttime = [];
                slp.endtime = [];
            end
            
            slp.time_immobile = time_immobile;    % in seconds
            slp.velocity_thresh = velocity_thresh;
            % install output in proper place
            sleep{day}{ep} = slp;
            clear slp;
        end
        
    end
    if ~bimodal_study
        if sleep_flag
                save(sprintf('%s/%ssleep%02d.mat', animaldir, animalprefix, day), 'sleep');
            clear sleep
        end        
        if rem_flag
            save(sprintf('%s/%srem%02d.mat', animaldir, animalprefix, day), 'rem');
            clear rem
        end
        if nrem_flag
            save(sprintf('%s/%snrem%02d.mat', animaldir, animalprefix, day), 'nrem');
            clear nrem
        end
        if sws_flag
            save(sprintf('%s/%ssws%02d.mat', animaldir, animalprefix, day), 'sws');
            clear sws
        end
        if spindle_flag
            save(sprintf('%s/%sspindles%02d.mat', animaldir, animalprefix, day), 'spindles');
            clear spindles
        end        
    end
end

end
