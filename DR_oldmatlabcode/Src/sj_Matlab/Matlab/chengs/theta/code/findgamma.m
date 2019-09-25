function findgamma
%  


nstd= 3;  % threshold= n*std+baseline
tetselid= 'most_cells'; % {'most_cells', 'EEGref'}
selectid= 'CA1PE';

allsel= selectid;
setRoot;

%  20-80Hz bandpass filter
load(sprintf('%s/data/gammafilter', root));

% count number of PE cells on CA1 tetrode
if strcmp(tetselid, 'most_cells')
    load(sprintf('%s/data/tetrodes-%s', root, allsel));
    alltet= tetrodes;
    ntet= length(alltet.rat);

    cl= collectCellList(allsel);
    ncl= length(cl.rat);
    ncells= zeros(ntet,1);

    for it=1:ntet
        rat= alltet.rat{it}; 
        d= alltet.num(it,1); e= alltet.num(it,2); t= alltet.num(it,3); 
        ind= find(strcmp(cl.rat', rat) & cl.cellnum(:,1)==d & ...
            cl.cellnum(:,2)==e & cl.cellnum(:,3)==t);
        ncells(it,:)= length(ind);
    end
end

epfile= sprintf('%s/data/epochs-%s', root, selectid);
load(epfile);
ep= epochs;
nep= length(ep.rat);


for ie=1:nep
    rat= ep.rat{ie}; d= ep.num(ie,1); e= ep.num(ie,2);
    fprintf(1, '%s [%d %d]\n', rat,d,e);
    monitorProgress(ie, nep);

    if ie==1 | ~strcmp(rat, ep.rat{ie-1})
        ripples= {};
    end
    % load behavior data
    if ie==1 | ~strcmp(rat, ep.rat{ie-1}) | d~= ep.num(ie-1,1)
        bfile= sprintf('%s/%s/data2/behavdata%.2d.mat',root,rat,d);
        load(bfile);
    end


    switch tetselid
    case 'most_cells'
        % find tetrode with most cells
        ind= find(strcmp(alltet.rat', rat) & alltet.num(:,1)==d & ...
                alltet.num(:,2)==e);
        [nmax,imax]= max(ncells(ind));
        tet= alltet.num(ind(imax),3);
    case 'EEGref'
        tet= ep.EEGref(ie,3);
    end

    % loadd EEG
    eegfile= sprintf('%s/%s/data/EEG/%seeg%.2d-%d-%.2d.mat',root,rat,rat,d,e,tet);
    load(eegfile);
    eval(sprintf('eegstruct= %seeg{%d}{%d}{%d};', rat, d, e, tet));

    time= behavdata{d}{e}.time;
    dt= mean(diff(time));
    nt= length(time);

    % band-pass filter EEG at high freq and detect ripples
    feeg= filtereeg(eegstruct, ripplefilter);
    rip= findripples(feeg, min_suprathreshold_duration, nstd, time);

    nrip= size(rip.data,1);
    fprintf(1, 'found %d ripples...\n', nrip);

    fs= eegstruct.samprate;
    mint= eegstruct.starttime;
    maxt= (length(eegstruct.data)-1)/fs+mint;
    t= linspace(mint, maxt, length(eegstruct.data))';

    lo=    floor((rip.data(:,1)-time(1))/dt);
    hi=     ceil((rip.data(:,2)-time(1))/dt);
    mid=   round((rip.data(:,3)-time(1))/dt);
    lo(lo<1)= 1; lo(lo>nt)= nt;
    hi(hi<1)= 1; hi(hi>nt)= nt;
    mid(mid<1)= 1; mid(mid>nt)= nt;
    tmp= false(nt,1);
    for ir=1:nrip
        tmp(lo(ir):hi(ir))= true;
    end

%    behavdata{d}{e}.ripple= tmp;
%    if ie==nep | ~strcmp(rat, ep.rat{ie+1}) | d~= ep.num(ie+1,1)
%        save(bfile, 'behavdata');
%    end
    ripples{d}{e}.selectid= selectid;
    ripples{d}{e}.min_suprathreshold_duration= min_suprathreshold_duration;
    ripples{d}{e}.nstd= nstd;
    ripples{d}{e}.eegfile= eegfile;
    ripples{d}{e}.tetselid= tetselid;
    ripples{d}{e}.tetrode= tet;
    ripples{d}{e}.nrip= nrip;
    ripples{d}{e}.lo= lo;
    ripples{d}{e}.hi= hi;
    ripples{d}{e}.mid= mid;
    ripples{d}{e}.rip= tmp;
    ripples{d}{e}.baseline= rip.baseline;
    ripples{d}{e}.threshold= rip.threshold;
    ripples{d}{e}.std= rip.std;
    ripples{d}{e}.peakheight= rip.data(:,4);

%    save(sprintf('%s/%s/data2/ripples_%dstd_%dms.mat', root,rat,nstd,min_suprathreshold_duration*1000), 'ripples');
    save(sprintf('%s/%s/data2/ripples_%s_%dstd_%dms.mat', root,rat,selectid,nstd,min_suprathreshold_duration*1000), 'ripples');
end
