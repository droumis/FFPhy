% analyze: get trial by trial firing rates
%adapted for mars 9/18/02cc
getcellsmars;
bdur=300;   %set bdur here; cells prior to 8/14/02 set bdur to 300ms
cueencode=1000;
ndatasets = length(dataset);
maxcells = 2;
% Epsx is the learning rate for the trial time function
Epsx =  5 * [1 1 1 1];
% Epst is the learning rate for the interspike interval function
Epst = 1 * [1 1 1 1];

% CP.x is the set of control points for the trial time function
CP.x = [0 0:50:4000 1e6 1e6];
% CP.t is the set of control points for the interspike interval function (not used at the moment, as Epst is set to 0)
CP.t = [0.000 0.001:0.004:0.025 0.050:0.025:10 1e6 1e6] * 1e3;

for i = 1:ndatasets
%for i=1:1
    dataset{i}.CP.x = CP.x;
    dataset{i}.CP.t = CP.t;    
    dataset{i}.Epsx = Epsx;    
    dataset{i}.Epst = Epst;
    cobj = dataset{i}.cobj;
    
    %% get the list of condition numbers
    dataset{i}.conds = unique(cobj.cond_no);
    nconds = length(dataset{i}.conds);
    
    if nconds~=4
        %display initial values for dataset{i}.conds(conds) and the number
        %of trials for each value of conds(#trials)
        %it is advised to only remove the cond that appears to have only
        %been presented once, since there is often one condition that was not actually
        %presented to the monkey but still appears as trial 1
        fprintf('%s\t %s\n','cond','#trials');
        for condno=1:nconds
            fprintf('%s\t\t %s\n',num2str(dataset{i}.conds(condno)),num2str(length(find(dataset{i}.cobj.cond_no==dataset{i}.conds(condno)))));
        end
        %dataset{i}.conds and nconds are updated to reflect user input
        dataset{i}.conds=input('Enter a vector with conditions to be analyzed: ');
        %ie[29 30 31 32]
        nconds = length(dataset{i}.conds);
    end
    
    %% get the list of correct trials
    correct = (cobj.response == 0);
    
    %% set up a cell array for the trials for each condition number
    %this is the same as condtrial in analyze.m
    dataset{i}.trialnum = cell(nconds,1);
    
    % get a list of the trials that do NOT contain conditions specified in dataset{i}.conds
    badcondtrialnum = zeros(size(cobj.cond_no));
    for trialn = 1:cobj.trials
        badcondtrialnum(trialn) = isempty(find(dataset{i}.conds(:)==cobj.cond_no(trialn)));
    end     
    
    %% get the total number of trials
    ntrials = length(cobj.trial_type);
    
    badtrial = zeros(ntrials,1);
    
    %%get a list of the fixation times (beginning of baseline) for each trial
    fixationtimes = zeros(1, ntrials);
    alltrialnum=[];
    for trialn = 1:ntrials,
        %% get rid of 
        %% 1. trials that do not contain condition numbers specified by the
        %% user--the user must use discretion in choosing these values(see
        %% above badcondtrialnum)
        %% 2. Trials where cue is not initiated
        if (badcondtrialnum(trialn)| isempty(find(cobj.codes(trialn,:)== cueencode))) 
            fixationtimes(trialn) = 0;
            % set this trial to be a bad trial
            badtrial(trialn) = 1;
        else
            %find index in c.codes where cue start is encoded
            fixationind(trialn) = max(find(cobj.codes(trialn,:)==cueencode));
            fixationtimes(trialn) = cobj.times(trialn,fixationind(trialn))- bdur;
        end
    end;
    
    badtrial(find((cobj.response ~= 0) & (cobj.response ~= 6))) = 1;
    % get a list of the spike times for each cell
    ncells = length(cobj.cluster_bank);
    spiketimes = cell(ncells,1);
    for cellnum = 1:ncells  
        % each cell's spike times are stored in one element of the spiketimes cell array
        spiketimes{cellnum} = ones(size(cobj.spikes)) * -1e10;
        for trialn = 1:ntrials
            for j = cobj.spike_start(trialn,cellnum):cobj.spike_end(trialn,cellnum)
                if (j > 0)
                    spiketimes{cellnum}(j) = double(cobj.spikes(j)) - fixationtimes(trialn);
                end;
            end
        end
    end
    % create an array for the values for the control point heights with one row for each trial
    
    
    dataset{i}.cell = cell(ncells,1);
    for cellnum = 1:ncells  
        dataset{i}.cell{cellnum}.trialtimerate = zeros(ntrials, length(CP.x));
        for condno = 1:nconds
            targetcode = dataset{i}.conds(condno);
            % the testID is 1 for the trials currently being analyzed and 0 for the rest of the trials
            % the code appears in column one instead of column 2 for mars'
            % cortex files
            testID = (cobj.codes(:,2) == targetcode+200);
            trials = find(cobj.cond_no == targetcode); 
            % get rid of trials where  There is a response other than 0 or 6 (correct or incorrect target)
            
            testID(find(badtrial)) = 0;
            % make the list of the trials for this condition number
            dataset{i}.trialnum{condno} = find(testID);
            % get the largest number of timesteps for any trial
            ntimesteps = max(cobj.times') - fixationtimes;
            % run the adaptive estimation program to get the firing rates
            timestep = 1;
            stmp = spiketimes{cellnum};
            [thetainit thetahat KSStat] = adaptivesuzuki(spiketimes{cellnum}, cobj.spike_start(:,cellnum), cobj.spike_end(:,cellnum), timestep, ntimesteps, fixationtimes, testID, Epsx, Epst, CP.x, CP.t);
            %[thetainit thetahat KSStat] = adaptivesuzuki(stmp, cobj.spike_start(:,cellnum), cobj.spike_end(:,cellnum), timestep, ntimesteps, fixationtimes, testID, Epsx, Epst, CP.x, CP.t);
            
            % get the list of times corresponding to the start time for each trial
            trialtimes = cumsum(ntimesteps);
            % get the list of valid trials (the trials for this condition number)
            validtrial = dataset{i}.trialnum{condno};
            % get the number of valid trials
            nvalidtrials = length(validtrial);
            alltrialnum=[alltrialnum,validtrial'];
            for trialn = 1:nvalidtrials 
                % we want the value of the curve at the end of each valid trial, so we get the time
                % at the beginning of the next trial unless we're at the last trial, in which case we
                % specify a time after the end of the last trial (instantspline will treat that as 
                % indicating that you want the last timestep).
                if (validtrial(trialn) < ntrials) 
                    timeindex = trialtimes(validtrial(trialn) + 1);
                else
                    timeindex = max(trialtimes) + 100000;
                end
                [dataset{i}.cell{cellnum}.trialtimerate(validtrial(trialn),:) tvals] = ...
                    instantspline(CP.x, CP.t, thetainit, thetahat, timeindex);
            end;
        end;
    end;
    dataset{i}.alltrialnum=sort(unique(alltrialnum));
end;        
%getareamars;


