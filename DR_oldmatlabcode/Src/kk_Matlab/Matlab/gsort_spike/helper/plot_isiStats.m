function plot_isiStats(spikes, clustnum, parameters)
%%% Similar to ISI quality: you are just plotting the statistics

if (isfield(spikes,'parameters'))
    tmin = spikes.parameters.tmin; tref = spikes.parameters.tref; 
elseif (length(parameters.tmin)~=0)
    tmin = parameters.min; tref = parameters.tref; 
else
    tmin=0.001; tref=0.002;
end

tmax=0.010; % 10ms
Fs=spikes.Fs; 
if Fs<100,
    Fs = Fs*1000; % convert to Hz if in kHz
end
useassigns = spikes.hierarchy.assigns;

members1 = find(useassigns == clustnum(1)); members2 = find(useassigns == clustnum(2));
unit1times = sort(spikes.swtimes(members1)); unit2times = sort(spikes.swtimes(members2));

%%%%%%%%%%%%%%%%%%%%% ISI QUALITY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful numbers
f = (tref - tmin) ./ (tmax - tmin);        % fraction of range that is below tref
bins = linspace(0, tmax, floor(tmax*Fs));  % using tmax*Fs bins means no info is lost to rounding
refractory_bin = sum((bins <= tref));      % index of refractory period

max_bin=sum((bins<=tmax));

% Convert from times to interspike intervals
isi1 = diff(sort(unit1times));
isi2 = diff(sort(unit2times));
isiT = diff(sort([unit1times; unit2times]));

% Calculate histograms for the intervals below tmax.
isi1hist = hist(isi1(isi1 < tmax), bins);
isi2hist = hist(isi2(isi2 < tmax), bins);
isiThist = hist(isiT(isiT < tmax), bins);

% Count total # intervals below tmax (set 0 counts to 1 as a courtesy to the next step)
isi1count = max(sum(isi1hist),1);
isi2count = max(sum(isi2hist),1);
isiTcount = max(sum(isiThist),1);

% Make cdfs from histograms (Fa, Fb and Fab in paper)
cdfs = zeros(3, length(bins));
cdfs(1,:) = cumsum(isi1hist) ./ isi1count;
cdfs(2,:) = cumsum(isi2hist) ./ isi2count;
cdfs(3,:) = cumsum(isiThist) ./ isiTcount;

% 
% if ((isi1count == 1) || (isi2count == 1)),  % if either original list had no violations
% 
% 	
% else
	% Compute the (scaled) difference between the initial cdfs and the
	% combined cdfs. (Sa in paper)
	diffs = zeros(2, length(bins));
	diffs(1,:) = sqrt((isi1count*isiTcount)./(isi1count+isiTcount)) .* (cdfs(3,:) - cdfs(1,:));
	diffs(2,:) = sqrt((isi2count*isiTcount)./(isi2count+isiTcount)) .* (cdfs(3,:) - cdfs(2,:));
	
	% Max value of this difference in the region shorter than the refractory period

	dstats(1) = max(diffs(1, 1:refractory_bin));
	dstats(2) = max(diffs(2, 1:refractory_bin));
	
	% based on Kolmogorov-Smirnoff statistics (see Fee, Mitra, Kleinfeld, 1996)
	lambda = 0:0.01:2;
	pdfDstat = (0.5 * (1 + erf(lambda./(sqrt(2*f*(1-f)))))) - ...
		(0.5 .* exp(-2.*lambda.^2) .* (1 - erf((1-2*f).*lambda./(sqrt(2*f*(1-f))))));
	%cutoff = lambda(max(find(pdfDstat < 0.95))); 
    cutoff = lambda(max(find(pdfDstat < 0.95))); %% SHANTANU - THIS RAISES CUTOFF
    %%cutff = 0.63 (tmin=1.2ms, tref=2ms), = 0.62 (tmin=1ms, tref=2ms) for 95% confidence interval

	% Does either d-statistic exceed the statistical cutoff?
	if (any(dstats > cutoff)),  allow = 0;
	else,                       allow = 1;
	end
% end

% The cdf at the tref bin is the fraction of spikes below tref; to get a score,
% we normalize by dividing by f from above.
scores =  (cdfs(:, refractory_bin)) ./ f; % AS IN PAPER, cdfs(:,tmax) normalized to 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FskHz=Fs/1000;
figure; redimscreen_hor;
subplot(1,2,1); hold on;
plot(cdfs(1,:),'Linewidth',2,'Linestyle','--'); 
plot(cdfs(2,:),'r','Linewidth',2,'Linestyle','--'); 
plot(cdfs(3,:),'k','Linewidth',2,'Linestyle','-');
axis([0 Fs*tmax 0 1]);
set(gca,'xtick',[0:FskHz:Fs*tmax], 'xticklabel',[{num2str([0:1:tmax*1000]')}]);
xlabel('ISI Time Interval (ms)'); ylabel ('Probability');
text(FskHz, 0.9,['isicount1: ' num2str(isi1count)]);
text(FskHz, 0.75,['isicount2: ' num2str(isi2count)]);
text(FskHz, 0.6,['isicountT: ' num2str(isiTcount)]);

subplot(1,2,2); hold on;
plot(diffs(1,:),'Linewidth',2,'Linestyle','--'); 
plot(diffs(2,:),'r','Linewidth',2,'Linestyle','--'); ;
axis([0 Fs*tmax -0.5 max([2,max(diffs)])]);
plot([0:FskHz:Fs*tmax], zeros(size([0:FskHz:Fs*tmax])), 'k','Linewidth',2 );
plot([0:FskHz:Fs*tmax], cutoff*ones(size([0:FskHz:Fs*tmax])), 'g','Linewidth',2,'Linestyle','--' );
plot(tref*Fs*ones(size([-0.5:0.1:max([2,max(diffs)])])), [-0.5:0.1:max([2,max(diffs)])], 'g','Linewidth',2 );
plot(tmin*Fs*ones(size([-0.5:0.1:max([2,max(diffs)])])), [-0.5:0.1:max([2,max(diffs)])], 'g','Linewidth',2,'Linestyle','--' );
set(gca,'xtick',[0:FskHz:Fs*tmax], 'xticklabel',[{num2str([0:1:tmax*1000]')}]);
xlabel('ISI Time Interval (ms)'); ylabel ('Probability');

text(FskHz*5, 1.5,['Scores:   ' num2str(roundn(scores',-2))]);


