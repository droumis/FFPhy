


%{

MEC theta phase at time of SWR across animals
- get corr coef of MU theta phase at swr onset as a function of swr #

- continue with what i was doing last friday
- get the theta range peak of MU in MEC
- convolve with a wavelet of that frequency
- get phase at swr start
- do a linear regression of theta phase as a function of swr # for each
animal

%}



animals = {'D10'};
filtfunction = 'dfa_riptrigspiking';
filenamesave = 'riptrigFRxdays';
env = 'wtrack';

loaddata = 1;
plotfigs = 0;
pausefigs = 0;
savefigs = 0;
plot_allntrodes = 0;

%% load data
if loaddata
    paths = make_paths(filtfunction, env);
    data = load_filter_output(paths.filtOutputDirectory, paths.filenamesave, ...
        animals);
    out = combine_riptrigspiking_perntrode(data);
end
%% For each ntrode, rip trig window of MU spiking, use the instaFR to run fft to
% get theta-range peak
% record that peak and the full fft with the data 
% create a wavelet at that peak and convolve the post-swr 
% do i want to run this on a larger window of data so I can get a more
% accurate estimate of peak frequency? I think maybe yes.. 
% i should also do the autocorr of all the spikes, >4cm/s and <4cm/s

% filtfunction = 'mua_calcxcorrmeasures';
% env = 'wtrack';
% paths = make_paths(filtfunction, env);
% muxcorr = load_filter_output(paths.filtOutputDirectory, paths.filenamesave, animals);
% i = muxcorr{1}.F.output{1}(229).ac1;
% 
% plot((i-mean(i)).^2/std(i).^2)
% xlim([480 520])
% right... it doesn't really make sense to do the autocorr with MU bc SO
% many spikes.. 
% so instead my using the instaFR from the riptrigspiking i think makes
% more sense.. just do that again using a larger spiking window -1 : 1 s
% alternatively i could get the peak theta instFR from all times <4cm/s,
% regardless of whether they are near a ripple or not.. I think just using
% a larger window right now is easier and will make me feel more confident
% in the peak.. 

% i should turn this whole process of finding peak theta, creating wavelet,
% convolving, create analytic signal for all rip trig traces

% so first start with a single epoch by hand:
% inputs would be index, excludetimes, datastructures
% look at the muxcorr script for guidance

% i could also just start with the theta phase that i think is already in
% the theta eeg, which is what i must have used to plot the phase
% precession

% GET THETA PHASE OF SWRs
% dfa_plotthetamod([2 2 11 12], [2 2 1], [], spikes, theta);
% iterator = 'singlecelleeganal';
% shantanu was getting theta modulation with 
% /home/droumis/Src/Matlab/DR_oldmatlabcode/DFSsj_plotthetamod_ver4.m

% Fp.animals = {'D10'};
% Fp = load_filter_params({'thetaphasemode', 'wtrack', '>4cm/s'});
% F = createfilter('animal',Fp.animals,'epochs',Fp.runepochfilter, 'excludetimefilter', ...
%     Fp.timefilter, 'cells', Fp.cellfilter,'eegtetrodes', Fp.eegfilter, 'iterator', ...
%     Fp.iterator);
% F = setfilterfunction(F,'dfa_plotthetamod',{'spikes','theta'});
%% 
animal = 'D10';
andef = animaldef(animal);

days = [1:12];
epochs = [2 4];

% day = 4;
% epoch = 4;

ntrode = 9;

excludetimes = [];

theta = loadeegstruct(andef{2}, andef{3}, 'theta', days, epochs, ntrode);
rips = loaddatastruct(andef{2}, andef{3}, 'ca1rippleskons', days);
%%
i = 0;
for d = 1:length(days)
    for e = 1:length(epochs)
        i = i +1;
        tind = [days(d) epochs(e) ntrode];
        riptimes = rips{tind(1)}{tind(2)}{1}.starttime;
        eegtime = geteegtimes(theta{tind(1)}{tind(2)}{tind(3)});
        thetaphase = theta{tind(1)}{tind(2)}{tind(3)}.data(:,2);
        ripphase = thetaphase(lookup(riptimes, eegtime));
        ripphase = double(ripphase) / 10000;
        indripphase{i} = [ripphase repmat(tind(1), length(ripphase),1) repmat(tind(2), length(ripphase),1) ...
            repmat(ntrode, length(ripphase),1)];
    end
end
r = cell2mat(indripphase');
r = [r (1:length(r(:,1)))'];
% r = [r; [r(:,1)+2*pi r(:,2:end)]];
% r = [r (1:length(r(:,1)))'];
% s = scatter(1:length(r),r(:,4),10,r(:,1), 'filled');
% s.LineWidth = 0.6;
% s.MarkerEdgeColor = 'k';
% s.MarkerEdgeAlpha = .2;
% s.MarkerFaceAlpha = .5;
%%

[counts, bin_centers] = hist3(r(:,[4 5]),[50 50],'CdataMode','auto');
xlabel('theta phase')
ylabel('swr#')
colorbar
view(2)
%%
figure
x_bin_centers = bin_centers{1};
y_bin_centers = bin_centers{2};
imagesc(x_bin_centers, y_bin_centers, counts)
hold on
contour(x_bin_centers, y_bin_centers, counts');
%%
% Generate some normally distributed data
x = randn(50,1);
y = randn(50,1);
% Estimate a continuous pdf from the discrete data
[pdfx xi]= ksdensity(r(:,4),'Bandwidth',0.3);
[pdfy yi]= ksdensity(r(:,5),'Bandwidth',0.3);
% Create 2-d grid of coordinates and function values, suitable for 3-d plotting
[xxi,yyi]     = meshgrid(xi,yi);
[pdfxx,pdfyy] = meshgrid(pdfx,pdfy);
% Calculate combined pdf, under assumption of independence
pdfxy = pdfxx.*pdfyy; 
% Plot the results
mesh(xxi,yyi,pdfxy)
set(gca,'XLim',[min(xi) max(xi)])
set(gca,'YLim',[min(yi) max(yi)])
view(2)
%%
imagesc(pdfxy)
contourf(xxi, yyi, pdfxy, 1)
%%
j = jet(12);

scatterhist(r(:,1), r(:,end),'Group',r(:,2),'Kernel','on','Location','NorthEast',...
    'Direction','out', 'Color', j(end:-1:1,:), 'MarkerSize', 7, 'Marker', '.'); %'Bandwidth', repmat(.6, 2, 12)
rotw = find(ismember(r(:,[2 3]), [9 4], 'rows'), 1);
line([-3.3 3.3], [rotw rotw], 'LineWidth', 2, 'color', 'k')
h = text(-3.5,rotw+50,'rotated wtrack');
set(h,'Rotation',90);
h = text(-3.5,1,'wtrack');
set(h,'Rotation',90);
set(gca,'Color', 'w')
xlabel('MEC theta phase')
ylabel('swr #')
title('rip eegthetaphase D10 nt9 wtrack')
axis tight
leg = legend();
title(leg,'day')

% get performance state at each swr time


%%
if ~isempty(riptimes)
    userips = ~isExcluded(riptimes, excludetimes);
else
    userips = [];
end

if ~isempty(userips)
    ripphase = thetaphase(lookup(riptimes, eegtime));
    ripphase = double(ripphase(userips)) / 10000;  % If no spikes, this will be empty
else
    ripphase = [];
end
numgoodrips = sum(userips);
scatter(riptimes, ripphase)
    
%%
% 
% 
% waveSet = 'riptrigMU_instantFR';
% paths = make_paths(waveSet, env);
% ntrodes = [out{1}.ntrode];
% alld = [];
% for n = 1:length(out{1})
%     nt_mu_keys = cell2mat({out{1}(n).mudata.index}');
%     nt_instaFR_alldayeps = permute(cell2mat({out{1}(n).mudata.instantFR}'), [3 2 1]);
%     eplengths = cellfun(@(x) length(x(:,1)), {out{1}(n).mudata.instantFR}, 'un', 1);
%     dayepvec = cell2mat(arrayfun(@(x,y,z) repmat([y z],x,1),eplengths',nt_mu_keys(:,1),nt_mu_keys(:,2),'un',0));
%     AS(n).ntrode = ntrodes(n);
%     AS(n).dayep_perrip = dayepvec;
%     
%     postwin = nt_instaFR_alldayeps(:,size(nt_instaFR_alldayeps,2)/2:end,:);
%     wp = getWaveParams(waveSet, postwin);
%     AS(n).postresult = computeAnalyticSignal(postwin, wp, animals{ian}, ...
%         paths.filenamesave, 'saveAnalyticSignal', 0);
%     
%     prewin = nt_instaFR_alldayeps(:,1:size(nt_instaFR_alldayeps,2)/2,:);
%     wp = getWaveParams(waveSet, prewin);
%     AS(n).preresult = computeAnalyticSignal(prewin, wp, animals{ian}, ...
%         paths.filenamesave, 'saveAnalyticSignal', 0);
%     
%     fullwin = nt_instaFR_alldayeps; %
%     wp = getWaveParams(waveSet, fullwin);
%     AS(n).result = computeAnalyticSignal(fullwin, wp, animals{ian}, ...
%         paths.filenamesave, 'saveAnalyticSignal', 0);
% end