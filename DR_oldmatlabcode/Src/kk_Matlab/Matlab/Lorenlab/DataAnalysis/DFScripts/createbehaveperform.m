% get the list of 1's and 0's for correct and incorrect trials for each epoch
%-----------------------------------------------------

% set the background probability of correct
backprob = 1/2;

animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Frank'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = [];
epochfilter{1} = ['isequal($environment, ''TrackA'')'];
epochfilter{2} = ['isequal($environment, ''TrackB'')'];

iterator = 'singleepochanal';

timefilter = { {'getlinvelocity', ['((abs($velocity) > -1))']} };

f = createfilter('animal',animals,'epochs',epochfilter, 'excludetimefilter', timefilter, 'iterator', iterator);


f = setfilterfunction(f, 'calcproprewarded', {'linpos', 'task'});
f = runfilter(f);

% concatenate the days and get the estimated probability of a correct response
% Track A

for a = 1:length(f)
    behavperform(1).task = 'TrackA';
    behavperform(2).task = 'TrackB';
    for t = 1:length(f(a).output)
	behavperform(t).outreward = [];
	behavperform(t).inreward = [];
	behavperform(t).dayintrials = zeros(length(f(a).output{t}),2);
	behavperform(t).dayouttrials = zeros(length(f(a).output{t}),2);
	behavperform(t).dayepoch = zeros(length(f(a).output{t}),2);
	intrialind = 1;
	outtrialind = 1;
	for i = 1:length(f(a).output{t})
	    b = f(a).output{t}(i);
	    outrew = b.correct(find(b.inbound == 0));
	    inrew = b.correct(find(b.inbound == 1));
	    ntout = length(outrew) - 1;
	    ntin = length(inrew) - 1;
	    behavperform(t).outreward = [behavperform(t).outreward; outrew];
	    behavperform(t).inreward = [behavperform(t).inreward; inrew];
	    behavperform(t).dayouttrials(i,:) = ...
	    				[outtrialind (outtrialind + ntout)];
	    behavperform(t).dayintrials(i,:) = [intrialind (intrialind + ntin)];
	    behavperform(t).dayepoch(i,:) = f(a).epochs{t}(i,:);
	    outtrialind = outtrialind + ntout + 1;
	    intrialind = intrialind + ntin + 1;
	end
	% get the estimated probability of a correct run
	behavperform(t).outprobcorrect = getestprobcorrect(...
	    behavperform(t).outreward, backprob, 1); 
	behavperform(t).inprobcorrect = getestprobcorrect(...
	    behavperform(t).inreward, backprob, 1); 
    end
    eval(sprintf('save %s%sbehavperform.mat behavperform', f(a).animal{2}, ...
                  f(a).animal{3}));
    clear behavperform;
end


