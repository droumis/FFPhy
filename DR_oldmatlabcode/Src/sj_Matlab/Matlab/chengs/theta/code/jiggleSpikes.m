function jiggleSpikes(stdtime, jigburst)
%function jiggleSpikes(stdtime, jigburst)
%
% stdtime:  std of Gaussian noise added to spike times
% jigburst: noise added to bursts, not individual spikes

if nargin<1; stdtime= 0.060; end  % = 60ms
if nargin<2; jigburst= 1; end 

ibi= 0.080; % [sec]

global behavdata spikedata fmaux

oldd= -1;
[d,e,t,c,ncells]= startCellList;
while ~isempty(d)

    if d~=oldd
        oldd= d;
        % load new datasets 
        loadVar('.', 'behavdata', d);
        loadVar('.', 'spikedata', d);
    end
    % set up variables
    time= behavdata{d}{e}.time;
    dt= (time(2)-time(1));
    if ~isempty(spikedata{d}{e}{t}{c})
        sp= spikedata{d}{e}{t}{c}.time;
        nsp= length(sp);
        ind= spikedata{d}{e}{t}{c}.index;

        %@@ consistency check
        if(sum((sp-time(ind))<dt*(1+1e-5)) ~= nsp);
            error('spike indices assigned incorrectly');
        end

        % jiggle spiketimes
        if ~jigburst
            sp= sp + stdtime*randn(nsp,1);
        else
            isi= diff(sp);
            bhi= [find(isi>= ibi)];
            blo= [1; bhi+1]; bhi= [bhi; nsp];
            nb= length(blo);
            for i=1:nb
                sp(blo(i):bhi(i))= sp(blo(i):bhi(i)) + stdtime*randn;
            end
        end

        % sort spiketimes, make sure times are valid
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

    [d,e,t,c]= getNextCell;
    if isempty(d) | (oldd~= d)
        save(sprintf('jiggled_%dms_spikes%.2d', round(stdtime*1000), oldd), 'spikedata', 'stdtime');
    end
end

