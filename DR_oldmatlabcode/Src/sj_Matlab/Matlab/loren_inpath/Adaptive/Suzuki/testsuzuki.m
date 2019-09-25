% adaptest: Runs the adaptive estimation routine.
% You need to have set the following variables:
% datasetnum
% cellnum    (the number of the cell in the dataset)
% targetcode (the condition number + 200)
% Epsx       (the learning rate for the trial time axis)
% Epst       (the learning rate for the interspike interval axis)
% 

Epsx = 5 * [1 1 1 1];
Epst = 0 * [1 1 1 1];
pushd Data
cobj = lorencortex('A020703.4', 1);
%cobj = lorencortex('A050801.4', 1);
%cobj = lorencortex('C110600.4', 1);
popd
timestep = 1; % 1 ms
ntrials = length(cobj.trial_type);
cellnum = 1;
CP.t = [0.000 0.001:0.004:0.025 0.050:0.025:10 1e6 1e6] * 1e3;
%CP.x = [0 0:50e-3:7100e-3 1e6 1e6];
CP.x = [0 0:50:7100 1e6 1e6];
%CP.x = [0 0:10:7100 1e6 1e6];

%C110600.4
%targetcode = 267;
%targetcode = 270;
%targetcode = 853;
%targetcode = 854;
%targetcode = 855;
%A050801.1
%targetcode = 579;
%targetcode = 580;
%targetcode = 581;
%targetcode = 554;
%A020703.4
targetcode = 513;
trials = find(cobj.cond_no == (targetcode - 200)); 
fixonly = zeros(size(cobj.cond_no));
for i = 1:cobj.trials
    fixonly(i) = ~isempty(find(cobj.codes(i,:) == 105));
end
%resp = cobj.response(trial);

correct = (cobj.response == 0);


testID = (cobj.codes(:,1) == targetcode);

%testID = (cobj.codes(:,1) >= 267 & cobj.codes(:,1) <=270) | (cobj.codes(:,1) >= 853 & cobj.codes(:,1) <= 855);

%testID = (cobj.codes(:,1) >= 267 & cobj.codes(:,1) <=270) + 2*(cobj.codes(:,1) >= 853 & cobj.codes(:,1) <= 855);

spiketimes = zeros(size(cobj.spikes));
fixationtimes = zeros(1, ntrials);
for i = 1:ntrials,
%    if isempty(find(cobj.codes(i,:)==23)),
    % get rid of 
    % 1. fixation only trials
    % 2. Trials where fixation does not occur (codes == 8) (or is it 23?)
    %if (fixonly(i) | isempty(find(cobj.codes(i,:)==8))) 
    %if (fixonly(i) | isempty(find(cobj.codes(i,:)==23))) 
    if (fixonly(i) | isempty(find(cobj.codes(i,:)==29))) 
	    fixationtimes(i) = 0;
	    testID(i) = 0;
    else,
        %fixationind(i) = max(find(cobj.codes(i,:)==8));
        %fixationind(i) = max(find(cobj.codes(i,:)==23));
        fixationind(i) = max(find(cobj.codes(i,:)==29));
	     fixationtimes(i) = cobj.times(i,fixationind(i))-300;
    end
    for j = cobj.spike_start(i,cellnum):cobj.spike_end(i,cellnum)
       if (j > 0) 
      	  spiketimes(j) = double(cobj.spikes(j,cellnum)) - fixationtimes(i);
       end
    end;
end;

% get rid of trials where  There is a response other than 0 or 6 (correct or incorrect target)
badtrial = find((cobj.response ~= 0) & (cobj.response ~= 6));
testID(badtrial) = 0;
ntimesteps = max(cobj.times') - fixationtimes;


correct = (cobj.response == 0);

[thetainit thetahat KSStat] = adaptivesuzuki(spiketimes, cobj.spike_start, cobj.spike_end, timestep, ntimesteps, fixationtimes, testID, Epsx, Epst, CP.x, CP.t);

sconvert


%figure(1);
%figure(2);
%'move figure 2'
%pause
%lastfive = zeros(5,1);
%for i=1:589,
%  if (testID(i)),
%    [xvals tvals] = instantspline(CP.x, CP.t, thetainit, thetahat, 7100*i);
%    [i thetahat(1,(i+1))]
%    for j = 1:4
%       lastfive(j) = lastfive(j+1);
%    end
%    if (cobj.response(i) ~= 0 ) 
%      lastfive(5) = 1;
%      figure(1);
%      plot(CP.x(1:end-2), xvals(1:end-2));
%      %axis([0 2 0 45]);
%      set(gca, 'XLim', [0 2]);
%      set(gca, 'XTick', [.3 .8 1.5], 'XTickLabel', char('300', '800', '1500'));
%      title(sprintf('%2f%% correct', sum(lastfive)/5 * 100));
%      figure(2);
%      semilogx(CP.t(1:end-2), tvals(1:end-2));
%    else
%      lastfive(5) = 0;
%      figure(1);
%      plot(CP.x(1:end-2), xvals(1:end-2), 'r');
%      set(gca, 'XLim', [0 2]);
%      %axis([0 2 0 45]);
%      set(gca, 'XTick', [.3 .8 1.5], 'XTickLabel', char('300', '800', '1500'));
%      title(sprintf('%2f%% correct', sum(lastfive)/5 * 100));
%      figure(2);
%      semilogx(CP.t(1:end-2), tvals(1:end-2), 'r');
%    end;
%  end;
%end; 
