function cmpRipples(rerun)
%  
%  Mark timesteps during which ripples occured.

showScatter= 0;
showDensity= 0;
showFrac= 0;
showMean= 1;

if nargin<1; rerun= 0; end

tetselid= 'most_cells'; % {'most_cells', 'EEGref'}

nstd= 3;  % threshold= n*std+baseline
min_suprathreshold_duration= 0.015; 
selectid= 'CA1PE';

allsel= selectid;
setRoot;

if rerun

    load(sprintf('%s/data/steve-ripple-filter', root));


    [tet, ep]= collectTet(selectid);

    nep= length(ep.rat);

    %keyboard

    for ie=1:nep
        rat= ep.rat{ie}; d= ep.num(ie,1); e= ep.num(ie,2);
        fprintf(1, '%s [%d %d]\n', rat,d,e);
        monitorProgress(ie, nep);

        if ie==1 | ~strcmp(rat, ep.rat{ie-1})
            load(sprintf('%s/%s/data2/ripples_%dstd_%dms.mat', root,rat,nstd,min_suprathreshold_duration*1000));
        end
        rip= ripples{d}{e};

        % load behavior data
        if ie==1 | ~strcmp(rat, ep.rat{ie-1}) | d~= ep.num(ie-1,1)
            bfile= sprintf('%s/%s/data2/behavdata%.2d.mat',root,rat,d);
            load(bfile);
        end
        time= behavdata{d}{e}.time;
        dt= mean(diff(time));
        nt= length(time);



        % find all recorded tetrodes
        ind= find(strcmp(rat, tet.rat) & d== tet.num(:,1)' & e== tet.num(:,2)');
        amp= zeros(rip.nrip, length(ind));

        % amplitudes of ripples on reference tetrode
        amp(:,1)= (rip.peakheight-rip.baseline)/rip.std; 
        j=1;

        % amplitudes of ripples on all other recorded tetrodes
        for itet= tet.num(ind,3)'
            if itet== ripples{d}{e}.tetrode; continue; end
            j= j+1;
        
            % load EEG
            eegfile= sprintf('%s/%s/data/EEG/%seeg%.2d-%d-%.2d.mat',root,rat,rat,d,e,itet);
            load(eegfile);
            eval(sprintf('eegstruct= %seeg{%d}{%d}{%d};', rat, d, e, itet));

            % band-pass filter EEG at high freq and detect ripples
            feeg= filtereeg(eegstruct, ripplefilter);
    %        rip= findripples(feeg, min_suprathreshold_duration, nstd, time);
            amp(:,j)= auxCalcAmp(feeg, time(rip.lo), time(rip.hi));
        end

        amplitudes.(rat){d}{e}= amp;
        save(sprintf('%s/work/ripple_coherence_%dstd_%dms.mat', ...
            root,nstd,min_suprathreshold_duration*1000), 'amplitudes');
    end

else
    % run analysis
    load(sprintf('%s/work/ripple_coherence_%dstd_%dms.mat', ...
        root,nstd,min_suprathreshold_duration*1000));

    if showDensity
        dens.x= [3:100];
        dens.y= [-1:100];
        dens= auxBin2d(dens);
    end

    rats={'kyl', 'ter', 'sta', 'fel'};
    ie= 0;
    for r=1:length(rats)
%    for r=2 %@@
        rat= rats{r};
        for d=1:length(amplitudes.(rat))
%        for d=8 %@@
            if isempty(amplitudes.(rat){d}); continue; end
            for e=1:length(amplitudes.(rat){d})
%            for e=4 %@@
                if isempty(amplitudes.(rat){d}{e}); continue; end
                ie= ie+1;
                a= amplitudes.(rat){d}{e};

                % make scatter plots of different amplitudes
                if showScatter
                    figure
                    for itet=2:min(7,size(a,2))
                        subplot(2,3,itet-1);
                        plot(a(:,1), a(:,itet),'.');
                    end
                    figname= sprintf('rip_coh_%s_%d_%d',rat,d,e);
                    set(gcf, 'Name', figname);
                    saveas(gcf, figname);

                    subplot(2,3,1); ylabel('other tetrodes');
                    subplot(2,3,4); ylabel('other tetrodes');
                    subplot(2,3,5); xlabel('reference tetrode');
                end

                % make summary density plot
                if showDensity
                    for itet=2:min(7,size(a,2))
                        dens= auxBin2d(dens,a(:,1), a(:,itet));
                    end
                end

                % calc fraction of p
                % calc mean of other tetrodes

                if showMean
                    ind= find(a(:,1) < 6); 
                    M{1}(ie)= mean(mean(a(ind,2:end)));
                    ind= find(a(:,1) >= 6 & a(:,1) < 9); 
                    M{2}(ie)= mean(mean(a(ind,2:end)));
                    ind= find(a(:,1) >= 9 & a(:,1) < 12); 
                    M{3}(ie)= mean(mean(a(ind,2:end)));
                end
            end
        end
    end
    if showMean
        M{1}= M{1}(isfinite(M{1}));
        M{2}= M{2}(isfinite(M{2}));
        M{3}= M{3}(isfinite(M{3}));
        figure

%        sel{1}.label= '3-6';
%        sel{2}.label= '6-9';
%        sel{3}.label= '9-12';
%        opt.plot= 'cdf'; opt.nbins= 15;
%        auxCmpDist2(M,sel,opt);

        opt.title= 'rip_coh_mean';
        opt.ylabel= {'mean peak amplitude', 'other tetrodes'};
        opt.xticklabels= {'3-6'  '9-12'  '9-12'}; 
        opt.plotcol= [0 0 0];
        for i=1:3
            Zm(i,1)= mean(M{i});
            Zdev(i,1)= std(M{i})/sqrt(length(M{i}));
        end
        auxBar(Zm, Zdev, opt);

%        Mall= []; g= [];
%        for i=1:3;
%            Mall= [Mall, M{i}];
%            g= [g, ones(1, length(M{i}))*i];
%        end
%        anovan(Mall, {g});

        keyboard
    end
    if showDensity
        auxBin2d(dens,'show');
    end
end

function amp = auxCalcAmp(ripplestruct, lot, hit)


% define the standard deviation for the Gaussian smoother which we
% apply before thresholding (this reduces sensitivity to spurious 
% flucutations in the ripple envelope)
SMOOTHING_WIDTH = 0.002; % 2 ms

% take the Hilbert-transform amplitude envelope of the ripple trace and
% smooth with a Gaussian with standard deviation of SMOOTHING_WIDTH
% which spans 4 SD (note that this uses the gausswin function
% of the Signal Processing toolbox)
gaussian_kernel = gausswin(ceil(8*SMOOTHING_WIDTH*ripplestruct.samprate),4);
% normalize
gaussian_kernel = gaussian_kernel/sum(gaussian_kernel);
% filter the ripple-filtered trace
smoothed_envelope = filtfilt(gaussian_kernel,[1],ripplestruct.data(:,3));

baseline= mean(smoothed_envelope);
dev= std(smoothed_envelope);

lo= floor((lot-ripplestruct.starttime)*ripplestruct.samprate)+1;
hi= floor((hit-ripplestruct.starttime)*ripplestruct.samprate)+1;

nrip= length(lo);
amp= zeros(nrip,1);
for ir=1:nrip
    amp(ir)= max((smoothed_envelope(lo(ir):hi(ir))-baseline)/dev);
end

