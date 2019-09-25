function jiggleUpto(selectid, timeid, stdtime, jigburst)
%function jiggleUpto(selectid, timeid, stdtime, jigburst)
%
% stdtime:  std of Gaussian noise added to spike times
% jigburst: noise added to bursts, not individual spikes

if nargin<3; stdtime= 0.060; end  % = 60ms
if nargin<4; jigburst= 1; end 

cutoff= 90; % [index into tm{iC}]
%cutoff= 1; % [index into tm{iC}]
ibi= 0.080; % [sec]
prefix= 'jiggledStep';
%prefix= 'jiggledStep2';
%prefix= 'jiggled1stPass';

global behavdata spikedata 

load(['/bach/theta/data/analist-' selectid])
nC= length(analist.rat);
tm= allTimes(timeid, selectid);
pf= allPlaceFields(selectid);

tmpn=0; %@@

oldrat= ''; oldd= -1;
for iC= 1:nC
    % set up data
    rat= analist.rat{iC};
    num= analist.cellnum(iC,:); d=num(1); e=num(2); t=num(3); c=num(4);
    if d~=oldd | ~strcmp(rat, oldrat)
        if ~isempty(oldrat)
            save(sprintf('/bach5/%s/data2/%s_%dms_spikes%.2d', oldrat, prefix, round(stdtime*1000), oldd), 'spikedata', 'stdtime');
        end
        oldd= d; oldrat= rat;
        % load new datasets 
        load(sprintf('/bach5/%s/data2/behavdata%.2d', rat, d));
        load(sprintf('/bach5/%s/data2/spikedata%.2d', rat, d));
        fprintf(1, '%s, day %2d\n', rat, d);
    end
    % set up variables
    if isempty(spikedata{d}{e}{t}{c}); continue; end

    time= behavdata{d}{e}.time;
    dt= (time(2)-time(1));

    sp= spikedata{d}{e}{t}{c}.time;
    nsp= length(sp);
    ind= spikedata{d}{e}{t}{c}.index;
    pos= spikedata{d}{e}{t}{c}.linpos;

    %@@ consistency check
    if(sum((sp-time(ind))<dt*(1+1e-5)) ~= nsp);
        error('spike indices assigned incorrectly');
    end

    %@@
%    if(length(tm{iC})<cutoff); continue; end
%    inpf= find(pf(iC,1)<= pos & pos <= pf(iC,2) & sp > tm{iC}(cutoff));

    maxt= min(length(tm{iC}),cutoff);
    if(maxt~= cutoff) tmpn=tmpn+1; end %@@
    inpf= find(pf(iC,1)<= pos & pos <= pf(iC,2) & sp < tm{iC}(maxt));

    insp= sp(inpf);
    nin= length(insp);


    % jiggle spiketimes
    if ~jigburst
        insp= insp + stdtime*randn(nin,1);
    else
        isi= diff(insp);
        bhi= [find(isi>= ibi)];
        blo= [1; bhi+1]; bhi= [bhi; nin];
        for i=1:length(blo)
            insp(blo(i):bhi(i))= insp(blo(i):bhi(i)) + stdtime*randn;
        end
    end

    % sort spiketimes, make sure times are valid
    sp(inpf)= insp;
    sp= sort(sp);
    sp(sp<time(1))= time(1);
    sp(sp>time(end))= time(end);

    % find new timeindices of spikes
    ind= findIndex(sp, time, 1);
    %@@ consistency check
    if(sum((sp-time(ind))<dt*(1+1e-5)) ~= nsp);
        error('spike indices assigned incorrectly');
    end

    % assign new spikes
    spikedata{d}{e}{t}{c}= [];
    spikedata{d}{e}{t}{c}.time= sp;
    spikedata{d}{e}{t}{c}.index= ind;
    spikedata{d}{e}{t}{c}.linpos= behavdata{d}{e}.linpos(ind);
    spikedata{d}{e}{t}{c}.phase = behavdata{d}{e}.phase(ind);

end

save(sprintf('/bach5/%s/data2/%s_%dms_spikes%.2d', rat, prefix, round(stdtime*1000), d), 'spikedata', 'stdtime');

fprintf(1, '%d/%d shorter than %d sec', tmpn, nC, cutoff);
