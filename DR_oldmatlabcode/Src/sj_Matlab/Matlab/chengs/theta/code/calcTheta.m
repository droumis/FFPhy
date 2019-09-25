function calcTheta
%function calcTheta
%  
%  Determine theta phase angle for all CA1 excitatory cells.
%  1. Determine EEG reference for each epoch as the tetrode with highest peak in
%     EEG power spectrum.
%  2. For every epoch determine phase angles at every timestep
%       a) filter EEG.
%       b) Hilbert transform.
%       c) calc phase
%  3. Determine offset between population activity of CA1 excitatory cell 
%     and EEG reference such that the max. of acivity occurs at 180deg. (This is
%     important for determining absolute phase) 


%auxEEGref;

selectid= 'CA1PE';
%selectid= 'pf7-fam';
mtm= 1;
xcwin= 0.080; % +/- [msec] window width for xcorr calc

setRoot;
epfile= sprintf('%s/data/epochs-%s', root, selectid);
load(epfile);
ep= epochs;
nep= length(ep.rat);

[cl, nc]= collectCellList(selectid);

load(sprintf('%s/data/thetafilter_6_10', root));
lf.b= thetafilter.tf.num;
lf.a= thetafilter.tf.den;

%freqz(lf.b,lf.a,512, 1500.7);
%keyboard

%@@
if(mtm)
    anafile= sprintf('%s/work/psana-pmtm-%s',root,selectid);
else
    anafile= sprintf('%s/work/psana-pwelch-%s',root,selectid);
end
load(anafile);

for ie=1:nep
    rat= ep.rat{ie}; d= ep.num(ie,1); e= ep.num(ie,2); 
    fprintf(1, 'calc. phase for %s [%d %d]\n', rat,d,e);

    if ie==1 | ~strcmp(rat, ep.rat{ie-1}) | d~= ep.num(ie-1,1)
        bfile= sprintf('%s/%s/data2/behavdata%.2d.mat',root,rat,d);
        load(bfile);
    end

    t= behavdata{d}{e}.time;
    dt= mean(diff(behavdata{d}{e}.time));
    nt= length(t);

    tet= ep.EEGref(ie,3);
    %@@
%    if ie==1 | ~strcmp(rat, ep.rat{ie-1}) 
%        load(sprintf('%s/%s/data/%sthetaelect', root, rat, rat));
%        eval(sprintf('tet= %sthetaelect;', rat));
%    end


    eegfile= sprintf('%s/%s/data/EEG/%seeg%.2d-%d-%.2d.mat',root,rat,rat,d,e,tet);
    load(eegfile);
    eval(sprintf('eegstruct= %seeg{%d}{%d}{%d};', rat, d, e, tet));

    eeg= eegstruct.data;
    fs= eegstruct.samprate;
    mint= eegstruct.starttime;
    maxt= (length(eeg)-1)/fs+mint;
    eegt= linspace(mint, maxt, length(eeg))';

%    b=fir1(200,.01);
%    eeg= filtfilt(b,1,eeg);

    Leeg= filtfilt(lf.b, lf.a, eeg);

    Heeg= hilbert(Leeg);
    phasetmp= angle(Heeg);
    phase= interp1(eegt,phasetmp,t);
    invalid= find(~isfinite(phase));
    if length(invalid)/nt>0.01; error('too many invalid timesteps in EEG'); end

%@@
%    eegstruct.data= Leeg;
%    phase= convertEEG2Phase(eegstruct, t);

%@@
%    thetafile= sprintf('%s/%s/data/EEG/%stheta%.2d-%d-%.2d.mat',root,rat,rat,d,e,tet);
%    load(thetafile);
%    eval(sprintf('theta= %stheta{%d}{%d}{%d};', rat,d,e,tet));
%    phase= convertEEG2Phase(theta, t);

%    clf
%    plot([Leeg, eeg]);
%    plot([Leeg, phasetmp/pi*200]);
%    plot([Leeg, phasetmp/pi*800, cos(phasetmp)*800]);
%    keyboard
%    bd= behavdata{d}{e};
%    plot(bd.time, [bd.phase, newphase]);

    clear eegstruct phasetmp
    phase(~isfinite(phase))= 0;

    mint= t(1);
    maxt= t(end)+dt-1e-10;

    if ie==1 | ~strcmp(rat, ep.rat{ie-1}) | d~= ep.num(ie-1,1)
        sdfile= sprintf('%s/%s/data2/spikedata%.2d.mat',root,rat,d);
        load(sdfile);
    end
    ind= find(strcmp(cl.rat, rat)' & cl.cellnum(:,1)==d & cl.cellnum(:,2)==e);

%    xcn= ceil(xcwin/dt);
%    count= zeros(nt, 1);
%    for i=ind'
%        st= cl.cellnum(i,3); sc= cl.cellnum(i,4); 
%        spind= find(behavdata{d}{e}.traj(spikedata{d}{e}{st}{sc}.index)>0);
%        if isempty(spind); continue; end
%        it= floor((spikedata{d}{e}{st}{sc}.time(spind)-mint)/dt);
%        it= it(it>0 & it<=nt);
%        count(it)= count(it)+1;
%    end
%    xc= xcorr(count, cos(phase)+1, xcn, 'unbiased');
%    [mxc, ixc]= max(xc);
%    dtxc= dt*(ixc-xcn-1);
%    fTheta= ana.theta(ep.EEGrefind(ie));
%    tTheta= 1/fTheta;
%    @@
%    tTheta= 0.125;
%    offset= 2*pi*dtxc/tTheta; 

    % calc. spike phases
    spikethetas= [];
    for i=ind'
        st= cl.cellnum(i,3); sc= cl.cellnum(i,4); 
        spind= spikedata{d}{e}{st}{sc}.index;
        spikephase= phase(spind);
        spikephase2= spikephase(behavdata{d}{e}.traj(spind)>0);
        spikethetas= [spikethetas; spikephase2];
    end
    offset= circstat(spikethetas);

%    figure; 
%    spikethetas= mod(spikethetas-offset+pi, 2*pi);
%    hist(spikethetas,20);
%    set(gca, 'XLim', [0,2*pi]);
%    set(gcf, 'Name', sprintf('%s %d %d', rat, d, e));
%    pause

    % so that peak of spiking occurs at pi
    behavdata{d}{e}.phase=  mod(phase-offset+pi, 2*pi);

    % recalc. spike phases
    for i=ind'
        st= cl.cellnum(i,3); sc= cl.cellnum(i,4); 
        spikedata{d}{e}{st}{sc}.phase= ...
            behavdata{d}{e}.phase(spikedata{d}{e}{st}{sc}.index);
    end


    if ie==nep | ~strcmp(rat, ep.rat{ie+1}) | d~= ep.num(ie+1,1)
        save(bfile, 'behavdata');
        save(sdfile, 'spikedata');
    end
 
end
