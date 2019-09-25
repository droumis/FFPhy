% get the list of 1's and 0's for correct and incorrect trials for each epoch
%-----------------------------------------------------

% set the background probability of correct
backprob = 0.5; %1/3;

%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
animals = {'HPa'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
dayfilter = '1:8';
epochfilter = [];
epochfilter{1} = ['isequal($environment, ''wtr1'')'];
epochfilter{2} = ['isequal($environment, ''wtr2'')'];

iterator = 'singleepochanal';

%timefilter = { {'sj_getlinvelocity', ['((abs($velocity) > -1))']} }; % Equivalent to all times, can skip
%f = createfilter('animal',animals, 'days',dayfilter, 'epochs',epochfilter, 'excludetimefilter', timefilter, 'iterator', iterator);

f = createfilter('animal',animals, 'days',dayfilter, 'epochs',epochfilter, 'iterator', iterator);


f = setfilterfunction(f, 'DFAsj_calcproprewarded', {'linpos', 'task'});
f = runfilter(f);

% concatenate the days and get the estimated probability of a correct response
% Track A

for a = 1:length(f)
    behavperform(1).task = 'wtr1';
    behavperform(2).task = 'wtr2';
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
	    behavperform(t).dayintime{i} = b.time(find(b.inbound == 1),:);
	    behavperform(t).dayouttime{i} = b.time(find(b.inbound == 0),:);
	    behavperform(t).dayepoch(i,:) = f(a).epochs{t}(i,:);
	    outtrialind = outtrialind + ntout + 1;
	    intrialind = intrialind + ntin + 1;
	end
	% get the estimated probability of a correct run
	behavperform(t).outprobcorrect = getestprobcorrect(...
	    behavperform(t).outreward, backprob, 0); 
	behavperform(t).inprobcorrect = getestprobcorrect(...
	    behavperform(t).inreward, backprob, 0); 
    end
    
    
    
    eval(sprintf('save %s%sbehavperform.mat behavperform', f(a).animal{2}, ...
                  f(a).animal{3}));
    clear behavperform;
end


