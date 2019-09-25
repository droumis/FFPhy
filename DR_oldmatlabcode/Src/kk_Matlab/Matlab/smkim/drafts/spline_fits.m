		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% common code
clear;
load('aligned_to_release_sort_by_RT.mat');

% define timebins of interest
BINWIDTH = 0.01;
timebins = BINWIDTH*(-50:100)';
numbins = numel(timebins)-1;

% what is maximum range of firing rates?
MAXRATE = 25;

% define knots of the stimulus-dependent spline
knots = (-0.75:0.25:1.5)';
splinemat = [ ...
	-0.5 1.5 -1.5 0.5; ...
	1 -2.5 2 -0.5; ...
	-0.5 0 0.5 0; ...
	0 1 0 0 ...
	];
% create a look-up table of values to calculate points along 
% the stimulus-dependent spline
spline_LUT = zeros(numbins,numel(knots));
for j = 1:numbins
	% identify spline knots around each time bin
	[junk, knotidx] = histc(timebins(j),knots);
	tmp = (timebins(j) - knots(knotidx))./(knots(knotidx+1) - knots(knotidx));
	tmp = [tmp.^3 tmp.^2 tmp.^1 1] * splinemat;
	spline_LUT(j,(knotidx-1):(knotidx+2)) = tmp;
end

% keep track of which trials we used for training
testidx{1}{1} = [];
testidx{2}{1} = [];
testidx{2}{2} = [];
testidx{3}{1} = [];
testidx{3}{2} = [];
testidx{3}{3} = [];
% and keep track of which trial type this was
trialtype{1}{1} = [];
trialtype{2}{1} = [];
trialtype{2}{2} = [];
trialtype{3}{1} = [];
trialtype{3}{2} = [];
trialtype{3}{3} = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rewarded trials

% populate spikecounts, a vector of length numtrials*numbins
spikecounts_r = [];
% s and t these are two *linked vectors*
% for rewarded trials
s_r = [1; 2; 3];
t_r = [1; 2; 3];
color_r = [0 0 1];
% counter for the number of trials
numtrials_r = 0;
spikesbytrial_r = cell(0,1);

% fit rewarded trials
for i = 1:numel(s_r)
	for j = 1:numel(aligned.data{s_r(i)}{t_r(i)})
		% coin flip to see whether this gets included in the training set
		if rand(1) < 0.8
			numtrials_r = numtrials_r + 1;
			% spikesbytrial{j} holds all spike times for the jth trial
			% we need to save this for the time-rescaling test later
			spikesbytrial_r{numtrials_r} = aligned.data{s_r(i)}{t_r(i)}(j).spikes;

			% now scan through spikesbytrial{end} and bin counts
			for k = 1:numbins
				spikecounts_r(end+1,1) = nnz( ...
					spikesbytrial_r{numtrials_r} >= timebins(k) & ...
					spikesbytrial_r{numtrials_r} < timebins(k+1) );
			end
		else
			testidx{s_r(i)}{t_r(i)}(end+1) = j;
			% rewarded trials are coded as 1
			trialtype{s_r(i)}{t_r(i)}(end+1) = 1;
		end
	end
end

% populate stimulus inputs, a matrix of size [numel(spikecounts), numel(knots)]
inputs_r = repmat(spline_LUT,numtrials_r,1);

% estimate parameters
[b_r,dev_r,stats_r] = glmfit(inputs_r,spikecounts_r,'poisson');

%

% stimulus-dependent component alone
% construct bar graph of raw spike counts
figure;
set(gcf,'Color','w');
%
rate_hist = reshape(spikecounts_r,numbins,numtrials_r)';
rate_hist = mean(rate_hist,1)/BINWIDTH;
h = bar(1e3*timebins,[rate_hist 0],'histc');
set(h,'FaceColor',0.2*color_r,'FaceAlpha',0.8,'EdgeColor',0.2*color_r);
%
% overlay spline confidence intervals
%
[lambda_r_hat,lambda_r_hi,lambda_r_lo] = glmval(...
	b_r(1:(1+numel(knots))),spline_LUT,'log',stats_r,'confidence',0.99);
lambda_r_hat = lambda_r_hat/BINWIDTH;
lambda_r_hi = lambda_r_hi/BINWIDTH;
lambda_r_lo = lambda_r_lo/BINWIDTH;
patch(1e3*[timebins(1:end-1); timebins(end-1:-1:1)], ...
	[lambda_r_hat-lambda_r_lo; flipud(lambda_r_hat+lambda_r_hi)], ...
	color_r,'FaceAlpha',0.2,'EdgeColor','none');
% zero alignment line
line([0 0],[0 MAXRATE],'Color','m','LineWidth',3);
% draw spline fit
line(1e3*timebins(1:end-1),lambda_r_hat,'Color',color_r,'LineWidth',3); 
% mark knots
line(1e3*knots,exp(b_r(1)+b_r(2:end))/BINWIDTH,...
	'Color',color_r,'LineStyle','none','MarkerSize',15,'Marker','.')
% label plot
set(gca,'XLim',1e3*[knots(2) knots(end-2)],'YLim',[0 MAXRATE],'Box','on',...
	'FontSize',30,'LineWidth',3,'TickDir','out','TickLength',[0.02 0.02]);
xlabel('time after lever release (ms)');
ylabel('firing rate (spikes/s)');

%
% time-rescaling test on interspike intervals
z_r = [];
for i = 1:numtrials_r
	rescaled = [];
	idxs = find( ...
		spikesbytrial_r{i} >= timebins(1) & ...
		spikesbytrial_r{i} < timebins(end) );
	if idxs(1) == 1
		idxs(1) = [];
	end
	% identify pairs of consecutive spike times
	[junk, bins2] = histc(spikesbytrial_r{i}(idxs),timebins); bins2(end) = [];
	[junk, bins1] = histc(spikesbytrial_r{i}(idxs-1),timebins); bins1(end) = [];
	if ~isempty(bins1) && ~isempty(bins2) && all(size(bins1)==size(bins2))
		if bins1(1) == 0
			bins1(1) = [];
			bins2(1) = [];
		end
		for j = 1:numel(bins2)
			rescaled(end+1,1) = BINWIDTH*sum(lambda_r_hat(bins1(j):bins2(j)));
		end
	end
	z_r = [z_r; 1 - exp(-rescaled)];
end
z_r = sort(z_r);

% K-S plot on rescaled intervals
figure;
set(gcf,'Color','w');
modelz_r = ((1:numel(z_r))'-0.5)/numel(z_r);
line(z_r,modelz_r,'LineWidth',3,'Color',color_r);
line(modelz_r,modelz_r,'LineStyle','--','Color',color_r);
conf95 = 1.36/sqrt(numel(modelz_r));
patch([modelz_r; flipud(modelz_r)],...
	[modelz_r + conf95; flipud(modelz_r-conf95)],...
	color_r,'FaceAlpha',0.2,'EdgeColor','none');
set(gca,'DataAspectRatio',[1 1 1],'XLim',[0 1],'YLim',[0 1],'Box','on',...
	'FontSize',30,'LineWidth',3,'TickDir','out','TickLength',[0.02 0.02]);
xlabel('observed quantiles');
ylabel('model predicted quantiles');
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rewarded trials

% populate spikecounts, a vector of length numtrials*numbins
spikecounts_n = [];
% s and t these are two *linked vectors*
% for rewarded trials
s_n = [2; 3; 3];
t_n = [1; 1; 2];
color_n = [0 0 0];
% counter for the number of trials
numtrials_n = 0;
spikesbytrial_n = cell(0,1);

% fit rewarded trials
for i = 1:numel(s_n)
	for j = 1:numel(aligned.data{s_n(i)}{t_n(i)})
		% coin flip to see whether this gets included in the training set
		if rand(1) < 0.8
			numtrials_n = numtrials_n + 1;
			% spikesbytrial{j} holds all spike times for the jth trial
			% we need to save this for the time-rescaling test later
			spikesbytrial_n{numtrials_n} = aligned.data{s_n(i)}{t_n(i)}(j).spikes;

			% now scan through spikesbytrial{end} and bin counts
			for k = 1:numbins
				spikecounts_n(end+1,1) = nnz( ...
					spikesbytrial_n{numtrials_n} >= timebins(k) & ...
					spikesbytrial_n{numtrials_n} < timebins(k+1) );
			end
		else
			testidx{s_n(i)}{t_n(i)}(end+1) = j;
			% rewarded trials are coded as 1
			trialtype{s_n(i)}{t_n(i)}(end+1) = 1;
		end
	end
end

% populate stimulus inputs, a matrix of size [numel(spikecounts), numel(knots)]
inputs_n = repmat(spline_LUT,numtrials_n,1);

% estimate parameters
[b_n,dev_n,stats_n] = glmfit(inputs_n,spikecounts_n,'poisson');

%

% stimulus-dependent component alone
% construct bar graph of raw spike counts
figure;
set(gcf,'Color','w');
%
rate_hist = reshape(spikecounts_n,numbins,numtrials_n)';
rate_hist = mean(rate_hist,1)/BINWIDTH;
h = bar(1e3*timebins,[rate_hist 0],'histc');
set(h,'FaceColor',0.2*color_n,'FaceAlpha',0.8,'EdgeColor',0.2*color_n);
%
% overlay spline confidence intervals
%
[lambda_n_hat,lambda_n_hi,lambda_n_lo] = glmval(...
	b_n(1:(1+numel(knots))),spline_LUT,'log',stats_n,'confidence',0.99);
lambda_n_hat = lambda_n_hat/BINWIDTH;
lambda_n_hi = lambda_n_hi/BINWIDTH;
lambda_n_lo = lambda_n_lo/BINWIDTH;
patch(1e3*[timebins(1:end-1); timebins(end-1:-1:1)], ...
	[lambda_n_hat-lambda_n_lo; flipud(lambda_n_hat+lambda_n_hi)], ...
	color_n,'FaceAlpha',0.2,'EdgeColor','none');
% zero alignment line
line([0 0],[0 MAXRATE],'Color','m','LineWidth',3);
% draw spline fit
line(1e3*timebins(1:end-1),lambda_n_hat,'Color',color_n,'LineWidth',3); 
% mark knots
line(1e3*knots,exp(b_n(1)+b_n(2:end))/BINWIDTH,...
	'Color',color_n,'LineStyle','none','MarkerSize',15,'Marker','.')
% label plot
set(gca,'XLim',1e3*[knots(2) knots(end-2)],'YLim',[0 MAXRATE],'Box','on',...
	'FontSize',30,'LineWidth',3,'TickDir','out','TickLength',[0.02 0.02]);
xlabel('time after lever release (ms)');
ylabel('firing rate (spikes/s)');

%
% time-rescaling test on interspike intervals
z_n = [];
for i = 1:numtrials_n
	rescaled = [];
	idxs = find( ...
		spikesbytrial_n{i} >= timebins(1) & ...
		spikesbytrial_n{i} < timebins(end) );
	if idxs(1) == 1
		idxs(1) = [];
	end
	% identify pairs of consecutive spike times
	[junk, bins2] = histc(spikesbytrial_n{i}(idxs),timebins); bins2(end) = [];
	[junk, bins1] = histc(spikesbytrial_n{i}(idxs-1),timebins); bins1(end) = [];
	if ~isempty(bins1) && ~isempty(bins2) && all(size(bins1)==size(bins2))
		if bins1(1) == 0
			bins1(1) = [];
			bins2(1) = [];
		end
		for j = 1:numel(bins2)
			rescaled(end+1,1) = BINWIDTH*sum(lambda_n_hat(bins1(j):bins2(j)));
		end
	end
	z_n = [z_n; 1 - exp(-rescaled)];
end
z_n = sort(z_n);

% K-S plot on rescaled intervals
figure;
set(gcf,'Color','w');
modelz_n = ((1:numel(z_n))'-0.5)/numel(z_n);
line(z_n,modelz_n,'LineWidth',3,'Color',color_n);
line(modelz_n,modelz_n,'LineStyle','--','Color',color_n);
conf95 = 1.36/sqrt(numel(modelz_n));
patch([modelz_n; flipud(modelz_n)],...
	[modelz_n + conf95; flipud(modelz_n-conf95)],...
	color_n,'FaceAlpha',0.2,'EdgeColor','none');
set(gca,'DataAspectRatio',[1 1 1],'XLim',[0 1],'YLim',[0 1],'Box','on',...
	'FontSize',30,'LineWidth',3,'TickDir','out','TickLength',[0.02 0.02]);
xlabel('observed quantiles');
ylabel('model predicted quantiles');
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save data for decoding phase
save('training.mat','testidx','trialtype','BINWIDTH','timebins','b_r','b_n','spline_LUT');

