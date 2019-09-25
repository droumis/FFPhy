% analyze: get trial by trial firing rates

getcells;
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

%for i = 1:ndatasets
for i=1:1
    dataset{i}.CP.x = CP.x;
    dataset{i}.CP.t = CP.t;    
    dataset{i}.Epsx = Epsx;    
    dataset{i}.Epst = Epst;
    cobj = dataset{i}.cobj;
    % get the list of condition numbers
    dataset{i}.conds = unique(cobj.cond_no);
    nconds = length(dataset{i}.conds);
    % get the list of correct trials
    correct = (cobj.response == 0);
    % set up a cell array for the trials for each condition number
    dataset{i}.condtrial = cell(nconds,1);
   
    % get a list of the trials that were fixation only
    fixonly = zeros(size(cobj.cond_no));
    for trialnum = 1:cobj.trials
        fixonly(trialnum) = ~isempty(find(cobj.codes(trialnum,:) == 105));
    end     
    % get the total number of trials
    ntrials = length(cobj.trial_type);
    
    
    badtrial = zeros(ntrials,1);
    %get a list of the fixation times for each trial
    fixationtimes = zeros(1, ntrials);
    for trialnum = 1:ntrials,
        % get rid of 
        % 1. fixation only trials
        % 2. Trials where fixation doefixas not occur (codes == 8) (or is it 23?)
        if (fixonly(trialnum) | isempty(find(cobj.codes(trialnum,:)== 29))) 
            fixationtimes(trialnum) = 0;
            % set this trial to be a bad trial
            badtrial(trialnum) = 1;
        else
            fixationind(trialnum) = max(find(cobj.codes(trialnum,:)==29));
            fixationtimes(trialnum) = cobj.times(trialnum,fixationind(trialnum))- 300;
        end
    end;
    badtrial(find((cobj.response ~= 0) & (cobj.response ~= 6))) = 1;
    % get a list of the spike times for each cell
    ncells = length(cobj.cluster_bank);
    spiketimes = cell(ncells,1);
    for cellnum = 1:ncells  
        % each cell's spike times are stored in one element of the spiketimes cell array
        spiketimes{cellnum} = ones(size(cobj.spikes)) * -1e10;
        for trialnum = 1:ntrials
            for j = cobj.spike_start(trialnum,cellnum):cobj.spike_end(trialnum,cellnum)
                if (j > 0)
                    spiketimes{cellnum}(j) = double(cobj.spikes(j)) - fixationtimes(trialnum);
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
            testID = (cobj.codes(:,1) == targetcode + 200);
            trials = find(cobj.cond_no == targetcode); % for condition number 853
            % get rid of trials where  There is a response other than 0 or 6 (correct or incorrect target)
            
            testID(find(badtrial)) = 0;
            % make the list of the trials for this condition number
            dataset{i}.condtrial{condno} = find(testID);
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
            validtrial = dataset{i}.condtrial{condno};
            % get the number of valid trials
            nvalidtrials = length(validtrial);
            
            for trialnum = 1:nvalidtrials 
                % we want the value of the curve at the end of each valid trial, so we get the time
                % at the beginning of the next trial unless we're at the last trial, in which case we
                % specify a time after the end of the last trial (instantspline will treat that as 
                % indicating that you want the last timestep).
                if (validtrial(trialnum) < ntrials) 
                    timeindex = trialtimes(validtrial(trialnum) + 1);
                else
                    timeindex = max(trialtimes) + 100000;
                end
                [dataset{i}.cell{cellnum}.trialtimerate(validtrial(trialnum),:) tvals] = ...
                    instantspline(CP.x, CP.t, thetainit, thetahat, timeindex);
            end;
        end;
    end;
end;        



