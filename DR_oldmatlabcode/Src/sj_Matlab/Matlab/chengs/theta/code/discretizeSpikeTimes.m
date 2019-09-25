function discretizeSpikeTimes(timestep, jigburst, stdtime)
%function discretizeSpikeTimes(timestep, jigburst)
%
% timestep  [msec]
% jigburst: noise added to bursts, not individual spikes

if nargin<1; time= 0.002; end  % = 2ms
if nargin<2; jigburst= 1; end 
if nargin<3; 
    jitter= 0; stdtime= []; 
else
    jitter= 1;
end 

ibi= 0.080; % [sec]
tiny= 1e-8;

global behavdata spikedata fmaux

fprintf(1, 'Working in %s\n', pwd);
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
    if ~isempty(spikedata{d}{e}{t}{c})
        sp= spikedata{d}{e}{t}{c}.time;
        sp= sp(find(time(1)<sp & sp<time(end)));
        nsp= length(sp);

        % jiggle spiketimes
        if ~jigburst
            if jitter; sp= sp + stdtime*randn(nsp,1); end
            sp= sp-mod(sp-time(1), timestep)+tiny;
        else
            isi= diff(sp);
            bhi= [find(isi>= ibi)];
            blo= [1; bhi+1]; bhi= [bhi; nsp];
            nb= length(blo);
            for i=1:nb
                if jitter; sp(blo(i):bhi(i))= sp(blo(i):bhi(i)) + stdtime*randn; end
                sp(blo(i):bhi(i))= sp(blo(i):bhi(i))-mod(sp(blo(i))-time(1), timestep)+tiny;
            end
        end

        % sort spiketimes, make sure times are valid
        sp= sort(sp);
        sp(sp<time(1))= time(1);
        sp(sp>time(end))= time(end);


        % find new timeindices of spikes
        ind= findIndex(sp, time, 1);

        % assign new spikes
        sd{d}{e}{t}{c}.time= sp;
        sd{d}{e}{t}{c}.index= ind;
        sd{d}{e}{t}{c}.linpos= behavdata{d}{e}.linpos(ind);
        sd{d}{e}{t}{c}.phase = behavdata{d}{e}.phase(ind);
    end

    [d,e,t,c]= getNextCell;
    if isempty(d) | (oldd~= d)
        spikedata= sd;
        if jitter
            fname= sprintf('disc_%dms_jit_%dms_spikes%.2d', ...
                round(timestep*1000), round(stdtime*1000), oldd);
            save(fname, 'spikedata', 'timestep', 'stdtime'); else
            fname= sprintf('discrete_%dms_spikes%.2d', round(timestep*1000), oldd);
            save(fname, 'spikedata', 'timestep');
        end
        load(fname); %% @@ test file integrity
        fprintf(1, 'saved file %s\n', fname); 
    end
end

