function [ccg, aux]= getCCG(sel, opt) 
%
%  Compute cross-correlograms for all pairs in sel.
% Code only works for timesteps of 2ms!!
%  Uses 

global behavdata spikedata 
setRoot

selectid= sel.selectid;
spsel= sel.spsel;
pairsel= sel.pairsel;

if isfield(opt, 'timeid'); 
    opt.t= sel.t; 
    workdir= sprintf('%s/work/CCG-%s/%s/%s/%s', ...
        root, opt.timeid, pairsel, spsel, selectid);
    fname= sprintf('%s/XC_%s_%s_%dms.mat', workdir, opt.scaling, sel.label, 1000*opt.dt);
else
    workdir= sprintf('%s/work/CCG/%s/%s/%s', ...
        root, pairsel, spsel, selectid);
    fname= sprintf('%s/XC_%s_%dms.mat', workdir, opt.scaling, 1000*opt.dt);
end

if  exist(fname, 'file')
    load(fname);
else
%    test whether directory exists
    if ~exist(workdir); error(['directory ' workdir ' does not exist']); end
    M= ceil(opt.T/opt.dt);
    aux.lags= [-M:M];
    scale= 0; pairs= 1;
    if strcmp(opt.scaling, 'coef'); scale=1 ; end
    if strcmp(sel.pairsel, 'none'); pairs=0 ; end
    aux.nrej= 0;

    fprintf(1, 'Calculating %s ...\n', fname);
    [cl, nc]= collectCellList(selectid);
    if(pairs) 
        [pl, np]= collectPairs(cl, pairsel, 0, strfind(spsel, 'placefield'));
    else
        pl= cl; np= nc;
    end
    aux.pairlist= pl; 
    ccg= nan*zeros(2*M+1,np);

    for ip=1:np

        % load spike times of 1st cell
        rat= pl.rat{ip};
        num= pl.cellnum(ip,:); d=num(1); e=num(2); tet=num(3); c=num(4);

        data2dir= fullfile(root,rat,'data2');
        loadVar(data2dir, 'behavdata', d);
        loadVar(data2dir, 'spikedata', d);

        sd{1}= spikedata{d}{e}{tet}{c};

        if(pairs)
            % load spike times of 2nd cell
            num2= pl.cellnum(ip,5:8); 
            d=num2(1); e=num2(2); tet=num2(3); c=num2(4);
            sd{2}= spikedata{d}{e}{tet}{c};
        end

%        [t, rej]= auxSelectSpikes(sd, behavdata{d}{e}, spsel, pl, ip, opt);
%        for i=1:2; t{i}= t{i}-behavdata{d}{e}.time(1); end
%        if(rej) continue; end
%        ccg(:,ip)= spikeCCG(t, M, opt.dt, scale)';
        if pairs; ncells= 2; else ncells= 1; end

%        valid= auxSelectSpikes2(sd, behavdata{d}{e}, spsel, pl, ip, opt);
        valid= auxSelectSpikes(rat, num, spsel, pl, ip, opt);
        for i=1:ncells
            tind{i}= sd{i}.index; 
            if strfind(spsel, 'min')
                nvalid= sum(valid(tind{i}));
                if nvalid<opt.minspikes;  
                    warning(sprintf('Reject pair with only  %d validspikes!',...
                            nvalid)); 
                    continue; 
                end
            end
        end
        ccg(:,ip)= spikeCCG(tind, M, scale, valid)';

%        plot(ccg(:,ip)); hold on, plot(ccg2(:,ip), 'r--'); hold off
%        pause
%        set(gca, 'xlim', [980 1022])
%        pause

        fprintf(1, '  %3d/%3d\n', ip, np);

    end
    aux.nrej= sum(isnan(ccg(1,:)));
    aux.npairs= np;
    save(fname, 'ccg', 'aux');
end

if isfield(sel, 'novelDay') & sel.novelDay
     ind= find(aux.pairlist.day==sel.novelDay);
    aux.npairs= length(ind);
else
    ind=[1:aux.npairs];
end

% clean-up
ind= ind(isfinite(ccg(1,ind)));
aux.ix= ind;
aux.nrej= aux.npairs-length(ind);
ccg= ccg(:,ind);
%keyboard

