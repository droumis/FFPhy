%function plotRipples
%  
%  Mark timesteps during which ripples occured.

nstd= 3;  % threshold= n*std+baseline
min_suprathreshold_duration= 0.015; 
selectid= 'CA1PE';

setRoot;
load(sprintf('%s/data/steve-ripple-filter', root));
load(sprintf('%s/data/epochs-%s', root, selectid));
nep= length(epochs.rat);


oldrat= '';
startie= 2; %%@@
for ie=startie:nep
    tic
    rat= epochs.rat{ie}; d= epochs.num(ie,1); e= epochs.num(ie,2);
    fprintf(1, '%s [%d %d]\n', rat,d,e);
    if ie==startie | ~strcmp(rat, epochs.rat{ie-1})
        load(sprintf('%s/%s/data2/ripples_%dstd_%dms.mat', root,rat,nstd,min_suprathreshold_duration*1000));
    end

    % load behavior data
    if ie==startie | ~strcmp(rat, epochs.rat{ie-1}) | d~= epochs.num(ie-1,1)
        bfile= sprintf('%s/%s/data2/behavdata%.2d.mat',root,rat,d);
        load(bfile);
    end

    % loadd EEG
    ripples{d}{e}
    tet= ripples{d}{e}.tetrode;
    load(ripples{d}{e}.eegfile);
    eval(sprintf('eegstruct= %seeg{%d}{%d}{%d};', rat, d, e, tet));

    toc
    time= behavdata{d}{e}.time;
    dt= mean(diff(time));
    ntime= length(time);

    fs= eegstruct.samprate;
    mint= eegstruct.starttime;
    maxt= (length(eegstruct.data)-1)/fs+mint;
    t= linspace(mint, maxt, length(eegstruct.data))';
    nt= length(t);
%    eeg= interp1(t, eegstruct.data, time);

    % band-pass filter EEG at high freq and detect ripples
    feeg= filtereeg(eegstruct, ripplefilter);

    toc
    baseline= ripples{d}{e}.baseline;
    dev= ripples{d}{e}.std;


    lor= floor((time(ripples{d}{e}.lo)-mint)*eegstruct.samprate);
    hir= ceil((time(ripples{d}{e}.hi)-mint)*eegstruct.samprate);
    nrip= ripples{d}{e}.nrip;
    if 1 
        still= nan*zeros(1,ntime);
        still(behavdata{d}{e}.traj<0)= 1;
        rip= nan*zeros(1,ntime);
        rip(ripples{d}{e}.rip)= 1;

        snippet_len= 0.5; % [sec]
        for ir=1:nrip
%            if (ripples{d}{e}.peakheight(ir)-baseline)/dev > 5; continue; end %%@@

            pad= round((snippet_len*eegstruct.samprate-(hir(ir)-lor(ir)))/2);
            pad= max(0,pad);
            lo= lor(ir)-pad; hi= hir(ir)+pad;
            lo(lo<1)= 1; lo(lo>nt)= nt;
            hi(hi<1)= 1; hi(hi>nt)= nt;
            hline= plot(t(lo:hi),eegstruct.data(lo:hi));
            set(hline, 'color', [.7 .7 .7], 'linewidth', .5);
            hold on
%            plot(t(lo:hi),feeg.data(lo:hi,3), 'y');
            hline= plot(t(lo:hi),feeg.data(lo:hi,1));
            set(hline, 'color', 'k', 'linewidth', .5);
            axis tight
            maxy= get(gca, 'ylim');

            pad= round((snippet_len/dt-ripples{d}{e}.hi(ir)+ripples{d}{e}.lo(ir))/2);
            pad= max(0,pad);
            lot= ripples{d}{e}.lo(ir)-pad;
            hit= ripples{d}{e}.hi(ir)+pad;
            lot(lot<1)= 1; lot(lot>ntime)= ntime;
            hit(hit<1)= 1; hit(hit>ntime)= ntime;
            plot(time(lot:hit), maxy(1)*still(lot:hit), 'b-', 'LineWidth', 2);
            plot(time(lot:hit), maxy(2)*rip(lot:hit), 'r-', 'LineWidth', 2);
%            h=line([t(lor(ir)) t(lor(ir))], get(gca, 'YLim'));
%            set(h, 'Color', 'r');
%            h=line([t(hir(ir)) t(hir(ir))], get(gca, 'YLim'));
%            set(h, 'Color', 'r');
            xlabel('time (s)')
            ylabel('EEG')
            hold off
            set(gca, 'ylim', maxy*1.1);
            fprintf(1,'%.2f std above baseline, ', (ripples{d}{e}.peakheight(ir)-baseline)/dev);
            fprintf(1,'key ...');
            pause
            fprintf(1,'\b\b\b\n');
        end
    else
        plot(t,feeg.data(:,1), 'g');  % filtered EEG
        hold on
    %    plot(t,feeg.data(:,3), 'y');  % envelope
        plot(t,eegstruct.data, 'k'); % raw EEG

        % not moving (v<2.4cm/s)
        tmp= nan*zeros(1,nt);
        tmp(behavdata{d}{e}.traj<0)= 0;
        plot(time, tmp, 'b-', 'LineWidth', 2);

        % ripples
        tmp= nan*zeros(1,nt);
        tmp(ripples{d}{e}.rip)= 10;
        plot(time, tmp, 'r-', 'LineWidth', 2);

        xlabel('time (s)')
        ylabel('EEG')
        hold off
    end
end
