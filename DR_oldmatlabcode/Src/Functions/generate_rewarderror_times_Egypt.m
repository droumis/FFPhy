% generate cell array of output and input trigger times on both rewarded
% and unrewarded trials


%Need in workspace:
%  rewardinfo and wellsdio from sj_findwellsfromdio1_Egypt (w-track) or
%  sj_findwellsfromdio1_Egypt_lin (linear track)

%% Manually enter these variables:

days = 7;                      %days to evaluate
epochs = [9];              %epochs to evaluate


%% For error times, select input times for incorrect trials (logic=0) 
errortimesall=cell(days(end),epochs(end));
errortimes=cell(days(end),epochs(end));
pseudotimes=cell(days(end),epochs(end));


for d=days
    for e=epochs
        inderr = find(rewardinfo{d}{e}(:,3) == 0);
        if ~isempty(inderr)
            errortimesall{d,e}=[errortimesall{d,e}; inderr rewardinfo{d}{e}(inderr,1) rewardinfo{d}{e}(inderr,2)];
            errortimes{d,e}=[errortimes{d,e}; errortimesall{d,e}(:,3)];
        end
        
    end
end

% (!!!) Evaluate up to here, then analyze errortimesall with position data to find
%exclusions: use Mari's plot_rewardtimes_position.m

%Find pseudotimes (last of reward output) after eliminating erroneous error
%times
for d=days
    for e=epochs
pseudotimes{d,e}=[pseudotimes{d,e}; errortimes{d,e}(:,1)+5000];
    end
end



%inderr = index of error trial in 'rewardinfo'
%rewardinfo{d}{e}(inderr,1) = error well (1=8,center ; 0=9 ; 2=10)
%rewardinfo{d}{e}(inderr,2) = input trigger time

%% Select input and output times for correct trials (logic=1)
rewardtimesall=cell(days(end),epochs(end));
testdiff=cell(days(end), epochs(end));
rewardtimes=cell(days(end), epochs(end));



for d=days
    for e=epochs
        ind = find(rewardinfo{d}{e}(:,3) == 1);
        if ~isempty(ind)
            rewardtimesall{d,e}=[rewardtimesall{d,e}; ind rewardinfo{d}{e}(ind,1) rewardinfo{d}{e}(ind,2) rewardinfo{d}{e}(ind,5)];
            testdiff{d,e} = [testdiff{d,e}; rewardtimesall{d,e}(:,3)-rewardtimesall{d,e}(:,4)];
            indright = find(testdiff{d,e}(:,1)>4900); %units are in 100s of microseconds
            if ~isempty(indright)    %Determine for correct trials how many outputs/inputs are not offset by 500ms, and exclude those   
                rewardtimes{d,e}=[rewardtimes{d,e}; rewardtimesall{d,e}(indright,3) rewardtimesall{d,e}(indright,4)];
                rewardtrialslost=length(ind)-length(indright)     %display trials excluded for each epoch
            else isempty(indright)
                rewardtimes{d,e}=[rewardtimes{d,e}; rewardinfo{d}{e}(ind,2) rewardinfo{d}{e}(ind,5)]; %for paradigms with no delay
            end
        end
    end
end

%ind = index of correct trial in 'rewardinfo'
%indright = index of output/input pairs offset by 500ms
%rewardinfo{d}{e}(ind,1) = correct well (1=8,center ; 0=9 ; 2=10)
%rewardinfo{d}{e}(ind,2) = output trigger time
%rewardinfo{d}{e}(ind,5) = input trigger time


disp('Generated errortimes, and rewardtimes with exclusions')

