function markRipples
%  
%  Mark timesteps during which ripples occured.

tetselid= 'most_cells'; % {'most_cells', 'EEGref'}

%threshold=400;
%baseline= 200;
min_suprathreshold_duration= 0.015; 
%if nargin<1; error('must specify selectid'); end
selectid= 'CA1PE';
%selectid= 'pf9-fam';
allsel= selectid;
setRoot;


% count number of PE cells on CA1 tetrode
%if strcmp(tetselid, 'most_cells')
%    load(sprintf('%s/data/tetrodes-%s', root, allsel));
%    alltet= tetrodes;
%    ntet= length(alltet.rat);

%    cl= collectCellList(allsel);
%    ncl= length(cl.rat);
%    ncells= zeros(ntet,1);

%    for it=1:ntet
%        rat= alltet.rat{it}; 
%        d= alltet.num(it,1); e= alltet.num(it,2); t= alltet.num(it,3); 
%        ind= find(strcmp(cl.rat', rat) & cl.cellnum(:,1)==d & ...
%            cl.cellnum(:,2)==e & cl.cellnum(:,3)==t);
%        ncells(it,:)= length(ind);
%    end
%end

epfile= sprintf('%s/data/epochs-%s', root, selectid);
load(epfile);
ep= epochs;
nep= length(ep.rat);

load(sprintf('%s/data/steve-ripple-filter', root));

for ie=1:nep
    rat= ep.rat{ie}; d= ep.num(ie,1); e= ep.num(ie,2);
    fprintf(1, '%s [%d %d]\n', rat,d,e);
    monitorProgress(ie, nep);

    % load behavior data
    if ie==1 | ~strcmp(rat, ep.rat{ie-1}) | d~= ep.num(ie-1,1)
        bfile= sprintf('%s/%s/data2/behavdata%.2d.mat',root,rat,d);
        load(bfile);
    end


    time= behavdata{d}{e}.time;
    dt= mean(diff(time));
    nt= length(time);

    % band-pass filter EEG at high freq and detect ripples
%    feeg= filtereeg(eegstruct, ripplefilter);
%    rip= findripples(feeg, threshold,baseline,min_suprathreshold_duration,time);
%    rip= findripples(feeg, min_suprathreshold_duration, time);

%    nrip= size(rip.data,1);
%    fprintf(1, 'found %d ripples...\n', nrip);

%    fs= eegstruct.samprate;
%    mint= eegstruct.starttime;
%    maxt= (length(eegstruct.data)-1)/fs+mint;
%    t= linspace(mint, maxt, length(eegstruct.data))';

%    lo= floor((rip.data(:,1)-time(1))/dt);
%    hi= ceil((rip.data(:,2)-time(1))/dt);
%    lo(lo<1)= 1; lo(lo>nt)= nt;
%    hi(hi<1)= 1; hi(hi>nt)= nt;
    ripple= false(nt,1);
    behavdata{d}{e}.ripple6= behavdata{d}{e}.ripple>0;


    if ie==nep | ~strcmp(rat, ep.rat{ie+1}) | d~= ep.num(ie+1,1)
        save(bfile, 'behavdata');
    end
end
