function showRipples(pairs, rerun, ripStr, varargin)
%
% Analyze ripples and neural activity during ripples.
%
% pairs: 
%     -1:compare epoch variables, not dependent on cells
%     0: compare individual cell firing (defined in id)
%     1: compare cell pairs with overlaping PF's (defined in id)
% vargin:
% Shantanu: Sen Chengs from Lorens directory

%ripStr= 'ripple';
%ripStr= 'nongamma3';
%ripStr= 'ECripple';
id= 'pf9';
dtBurst= 0.010; % 10 ms
opt.nmin= 5;
minrips= 5;
alpha= 0.05;

plotCoincNorm= 0; plotCoincRate= 0; plotCoincFrac= 0; plotRipRate= 0;
plotRipLength= 0; plotRipFrac= 0; plotOcc= 0; plotSpikeRate= 0; 
plotBurstRate= 0; plotSpikes_Burst= 0; plotNRips= 0; plotNBursts= 0;
plotNSpikes= 0; plotN12= 0; plotDN12= 0; plotZ= 0; plotZd= 0; plotNPairs= 0; 
plotRipVel= 0; plotVel= 0; plotActivation= 0; plotBurst_active= 0;
plotSpikes_active= 0; plotRipSpikeRatio= 0; plotAllSpikes= 0;

v= varargin;
for i=1:length(v)
    switch v{i}
    case 'coincNorm'
        plotCoincNorm= 1;
    case 'coincFrac'
        plotCoincFrac= 1;
%        plotDN12= 1;
    case 'z'
        plotZ= 1;
    case 'zd'
        plotZd= 1;
    case 'nPairs'
        plotNPairs= 1;
    case 'activation'
        plotActivation= 1;
    case 'burst_active'     % bursts per ripple activation
        plotBurst_active= 1;
    case 'spikes_active'    % spikes per ripple activation
        plotSpikes_active= 1;
    case 'ripSpikeRatio'     % ratio of spike w/in ripples to total # of ripples
        plotRipSpikeRatio= 1;
    case 'allSpikes'     % number of all spikes
        plotAllSpikes= 1;
    case 'ripRate'
        plotRipRate= 1;
    case 'ripLength'
        plotRipLength= 1;
    case 'ripFrac'
        plotRipFrac= 1;
    case 'spikeRate'
        plotSpikeRate= 1;
    case 'burstRate'
        plotBurstRate= 1;
    case 'spikes_burst'
        plotSpikes_Burst= 1;
    case 'nRips'
        plotNRips= 1;
    case 'nBursts'
        plotNBursts= 1;
    case 'nSpikes'
        plotNSpikes= 1;
    case 'vel'
        plotVel= 1;
    case 'ripVel'
        plotRipVel= 1;
    case 'occ'
        plotOcc= 1;
    otherwise 
        error('unknown option');
    end
end


global armId
if pairs>=0
    if pairs
        armId= {'famArm', 'novelArm', 'cross'};
%        armId= {'famArm', 'novelArm'}; %%@@
    else
        armId= {'famArm', 'novelArm'};
    end
    sets{1}.novelDay= 1;
    sets{1}.label= 'day1';
    sets{1}.arm= 'pf';  % {'famArm', 'novelArm', 'pf', 'nonpf'}
    sets{2}.novelDay= 2;
    sets{2}.label= 'day2';
    sets{2}.arm= 'pf';
    sets{3}.novelDay= 3;
    sets{3}.label= 'day3';
    sets{3}.arm= 'pf';

    sets{4}.novelDay= 1;
    sets{4}.label= 'day1';
    sets{4}.arm= 'nonpf';
    sets{5}.novelDay= 2;
    sets{5}.label= 'day2';
    sets{5}.arm= 'nonpf';
    sets{6}.novelDay= 3;
    sets{6}.label= 'day3';
    sets{6}.arm= 'nonpf';
else
    id= 'CA1PE';
    armId= {'famArm', 'novelArm'};
    sets{1}.novelDay= 1;
    sets{1}.label= 'day1';
    sets{1}.arm= 'pf';
    sets{2}.novelDay= 2;
    sets{2}.label= 'day2';
    sets{2}.arm= 'pf';
    sets{3}.novelDay= 3;
    sets{3}.label= 'day3';
    sets{3}.arm= 'pf';
end

na= length(armId);
for ia=1:na; sel{ia}.label= armId{ia}; end
b= {[-75:5:0], [0:5:75]};
%b= {[-75:7.5:0], [0:7.5:75]}; %@@
x= [b{1}(1:end-1) b{2}(1:end-1)]+75;
dx= mean(diff(x));
x= x+dx/2;

TRAJ= {[0 1], [2 3]};

% parsing argument list
if nargin<2; rerun= 0; end
nb(1)= length(b{1})-1; nb(2)= length(b{2})-1;
setRoot
if pairs>0; pairStr= 'overlap_'; else pairStr= ''; end

for is=1:length(sets)
    novelDay= sets{is}.novelDay;
    if pairs>=0
        optstr{is}= sprintf('%s%s_%s_day%d', pairStr, sets{is}.arm, id, novelDay);
    else
        optstr{is}= sprintf('%s%s_day%d', pairStr, id, novelDay);
%        optstr{is}= sprintf('%s%s_7_5_day%d', pairStr, id, novelDay); %%@@
    end
    fname{is}= sprintf('%s/work/%s/data_%s.mat', root, ripStr, optstr{is});


    dirname= sprintf('%s/work/%s', root, ripStr);
    if ~exist(dirname, 'dir')
        error(['dir ' dirname ' does not exist']);
    end

    clear rip ripRate count1 count2 count12 coincRate coincFrac
    clear nPairs iepoch
    clear coincNorm nRips nBursts nSpikes ripocc occ allSpikes
    ripVel= [];

    if ~rerun & exist(fname{is}, 'file')== 2
        continue;
    end
    fname{is}

    for ia=1:na
        if strfind(armId{ia}, 'cross') 
            [cl, nc]= mixPairs([id '-novelArm'], [id '-famArm']);
        else
            if strfind(id, 'pf') 
                selectid= [id '-' armId{ia}];
            else
                selectid= id;
            end
            [cl, nc]= collectCellList(selectid);
            if(pairs>0) [cl, nc]= collectPairs(cl, 'overlap', 0, 0); end
        end

        ie= 0; ix= 0;
        startic= 1;
        oldrat= ''; oldd= -1; olde= -1;

        aux{ia}.IP= [];

        for ic=startic:nc
            rat= cl.rat{ic};

            if(cl.day(ic)~= novelDay); continue; end
            d= cl.cellnum(ic,1); e= cl.cellnum(ic,2);
            t= cl.cellnum(ic,3); c= cl.cellnum(ic,4);

            if ~strcmp(rat, oldrat) | oldd~= d
                fprintf(1, '%s day %d\n', rat, d);
            end

            global behavdata spikedata lindistpos vel info
            data2dir= fullfile(root,rat,'data2');
            datadir= fullfile(root,rat,'data');

            loadVar(data2dir, 'behavdata', d);
            loadVar(data2dir, 'spikedata', d);
            loadVar(datadir, 'lindistpos', d, rat, 1);
            loadVar(data2dir, 'info', 0);
            if pairs<0
                loadVar(datadir, 'vel', d, rat, 1);
            end

            if ~strcmp(rat, oldrat) | oldd~= d | olde~= e
                oldrat= rat; oldd= d; olde= e;

                ie= ie+1;
                if (strcmp(sets{is}.arm,'pf') & ...
                        (strcmp(armId{ia}, 'novelArm') | strcmp(armId{ia}, 'cross'))) |...
                    (strcmp(sets{is}.arm,'nonpf') & strcmp(armId{ia}, 'famArm'))
                    % novel trajs
                    if cl.newarm(ic) < 5;
                        itraj=1; 
                    else
                        itraj=2;
                    end
                else 
                    % fam trajs
                    if cl.newarm(ic) < 5;
                        itraj=2;
                    else
                        itraj=1;
                    end
                end
                trajs= TRAJ{itraj};

                bd= behavdata{d}{e};
                center= info{d}{e}.centerlinpos;
                x= bd.linpos-center;
                dt= mean(diff(bd.time));
                trajvalid= ismember(lindistpos{d}{e}.estinfo, trajs);

                if pairs<0
                    velt= interp1(vel{d}{e}.data(:,1), vel{d}{e}.data(:,2), lindistpos{d}{e}.data(:,1));
                end

                ripple= auxGetRipples(rat, d, e, ripStr);
                ind= find(ripple & trajvalid);
                [lo hi]= findcontiguous(ind);
                nrip= length(lo);
                xrip= x(lo);

                for j=1:2
                    occtmp= saveHistc(x(trajvalid), b{j});
                    ripocctmp{j}= saveHistc(x(ind), b{j});

                    tmp= saveHistc(xrip, b{j});
                    nRips{ia}{j}(ie,:)= tmp(1:end-1);

                    occ{ia}{j}(ie,:)= occtmp(1:end-1)*dt;
                    occtmp(occtmp==0)= nan;
                    tmp= ripocctmp{j}(1:end-1)./occtmp(1:end-1);
%                        tmp(~isfinite(tmp))= 0;
                    ripocc{ia}{j}(ie,:)= tmp;

                    tmp= nRips{ia}{j}(ie,:)./occtmp(1:end-1);
%                        tmp(~isfinite(tmp))= 0;
                    ripRate{ia}{j}(ie,:)= tmp/dt;

                    if pairs>0; nPairs{ia}{j}(ie,:)= zeros(1,nb(j)); end

                    armind= ind(b{j}(1)<= x(ind) & x(ind) <= b{j}(end));
                    if pairs< 0
                        velocity{ia}{j}(ie)= mean(velt(armind));
                        if length(ripVel)<ia ripVel{ia}= {}; end
                        if length(ripVel{ia})<j ripVel{ia}{j}= []; end
                        inArmRip= find(b{j}(1)<=xrip & xrip<= b{j}(end));
                        for ir=1:length(inArmRip)
                            ripVel{ia}{j}(end+1)=  mean(velt(lo(inArmRip(ir)):hi(inArmRip(ir))));
                        end
                    end

                    inarm= b{j}(1)<=xrip & xrip<= b{j}(end);
                    rip{ia}{j}.length{ie}= (hi(inarm)-lo(inarm)+1)*0.002;
                end
            end

            ix= ix+1;

            aux{ia}.IP(ix)= ic;
            sd= spikedata{d}{e}{t}{c};
            spind= find(trajvalid(sd.index) & ripple(sd.index));
            ind= sd.index(spind);
            allind= sd.index(find(trajvalid(sd.index)));
            if pairs>0
                sd2= spikedata{d}{e}{cl.cellnum(ic,7)}{cl.cellnum(ic,8)};
                ind2= sd2.index(trajvalid(sd2.index) & ripple(sd2.index));
            end
            for j=1:2
                iepoch{ia}{j}(ix)= ie;
                sptmp= saveHistc(x(ind), b{j});
                nSpikes{ia}{j}(ix,:)= sptmp(1:end-1);

                tmp= saveHistc(x(allind), b{j});
                allSpikes{ia}{j}(ix,:)= tmp(1:end-1);

                tmp= find(diff(sd.time(spind))<dtBurst);
                indBurst= spind;
                indBurst(tmp+1)= nan;
                indBurst= indBurst(isfinite(indBurst));
                tmp= saveHistc(x(sd.index(indBurst)), b{j});
                nBursts{ia}{j}(ix,:)= tmp(1:end-1);

                found= false(1, nrip); 
                for ir=1:nrip
                     if find(lo(ir)<= ind & ind <=hi(ir)) found(ir)=true; end
                end
                tmp= saveHistc(xrip(found), b{j});
                count1{ia}{j}(ix,:)= tmp(1:end-1);


                if pairs>0
                    nPairs{ia}{j}(ie,:)= nPairs{ia}{j}(ie,:) + 1;
                    sptmp= saveHistc(x(ind2), b{j});
                    spcount2{ia}{j}(ix,:)= sptmp(1:end-1);
                    found2= false(1, nrip);
                    for ir=1:nrip; if find(lo(ir)<=ind2 & ind2<=hi(ir)) found2(ir)=true; end;  end
                    tmp= saveHistc(xrip(found2), b{j});
                    count2{ia}{j}(ix,:)= tmp(1:end-1);
                    tmp= saveHistc(xrip(found & found2), b{j});
                    count12{ia}{j}(ix,:)= tmp(1:end-1);

                    tmp= occ{ia}{j}(ie,:);
                    tmp(tmp==0)= nan;
                    tmp=  count12{ia}{j}(ix,:)./tmp;
%                        tmp(~isfinite(tmp))= 0;
                    coincRate{ia}{j}(ix,:)= tmp;

                    tmp= nRips{ia}{j}(ie,:);
                    tmp(tmp<5)= nan;
                    tmp=  count12{ia}{j}(ix,:)./tmp;
                    tmp(~isfinite(tmp))= 0;
                    coincFrac{ia}{j}(ix,:)= tmp;
%                        coincFrac{ia}{j}(ix,:)= count12{ia}{j}(ix,:)/ sum(nRips{ia}{j}(ie,:));

                    tmp= ripRate{ia}{j}(ie,:);
                    tmp(tmp==0)= nan;
                    tmp=  count12{ia}{j}(ix,:)./tmp;
%                        tmp(~isfinite(tmp))= 0;
                    coincNorm{ia}{j}(ix,:)= tmp;
                end
            end

%            if(ie>=2); break; end %%@@
        end
    end
    switch pairs
    case 1
        save(fname{is}, 'rip', 'nPairs', ...
        'count1', 'count2', 'count12', 'coincRate', 'coincFrac', 'coincNorm', ...
            'nRips', 'nBursts', 'nSpikes', 'ripocc', 'occ', ...
            'iepoch', 'armId', 'b', 'aux');
    case 0
        save(fname{is}, 'rip', 'count1', 'nRips', 'nBursts', 'nSpikes',...
        'allSpikes', 'ripocc', 'occ', 'iepoch', 'armId', 'b', 'aux');
    case -1
        save(fname{is}, 'rip', 'count1', 'nRips', 'nBursts', 'nSpikes', 'ripocc', 'occ', 'velocity', 'ripVel', 'iepoch', 'armId', 'b', 'aux');
    end
end % for is; do analysis

for is=1:length(sets)
    load(fname{is});
%    keyboard

%    plotcol= {[[0 0 0]; [1 0 0]], [[0 0 0]; [0 1 0]], [[0 0 0]; [0 0 1]]};
    plotcol= {[[1 1 1]; [1 0 0]], [[1 1 1]; [0 1 0]], [[1 1 1]; [0 0 1]]};
    opt.plotcol= plotcol{sets{is}.novelDay};

    opt.outname= [ripStr '_' optstr{is}];
%    opt.xrange= [0 150];

    if plotCoincRate
        opt.label= 'coincidence rate'; opt.title= 'ripCoincRate';
        opt.yrange= [-0.01 0.2];
        auxImage(b, coincRate, opt)
    end
    if plotCoincFrac
        opt.label= 'coincidence frac'; opt.title= 'ripCoincFrac';
%        opt.yrange= [-0.005 0.36]; 
%        auxImage(b, coincFrac, opt)
        for ia=1:na
            f= [];
            IP{is}{ia}= [];
            for j=2
%                num= sum(count12{ia}{j});
%                den= sum(nPairs{ia}{j}.*nRips{ia}{j});
                num= sum(count12{ia}{j},2);
                den= sum(nRips{ia}{j}(iepoch{ia}{j},:),2);
                ind= den>=minrips;

%                cl= collectCellList([id '-' armId{ia}]);
%                pl= collectPairs(cl, 'overlap');
%                ind= ind & strcmp({pl.rat{aux{ia}.IP}}, 'fel')';
%                fprintf(1, '%d, %d, %d \n', length(den), length(aux{ia}.IP), sum(ind));

                den(den<minrips)= nan;
                frac{ia}{j}= num./den;
                N{ia}{j}= den;
%                ci{ia}{j}= norminv(1-alpha/2)*sqrt(frac{ia}{j}.*(1-frac{ia}{j})./den);
                ci{ia}{j}= sqrt(frac{ia}{j}.*(1-frac{ia}{j})./den);
                f= [f; frac{ia}{j}(ind)];
                IP{is}{ia}= [ IP{is}{ia}, aux{ia}.IP(ind)];
            end
            COINCFRAC{is}{ia}= f;
        end
%        auxImage(b, frac, opt, ci)
%        auxAnova1(COINCFRAC{is}, sets{is}, opt);
%        opt.plotcol(3,:)= 0.7*[1 1 1];
%        auxCmpDist2(COINCFRAC{is}, sel, opt);
        save('data_COINCFRAC', 'COINCFRAC', 'IP');
    end
    if plotCoincNorm
        opt.label= 'coincNorm'; opt.title= 'ripCoincNorm';
        opt.yrange= [-1 33];
        auxImage(b, coincNorm, opt)
    end

    if plotZ | plotZd
        Z= cell(1,na);
        for ia=1:na
            IP{is}{ia}= [];
            for j=2
%                N= nRips{ia}{j}(iepoch{ia}{j},:);
%                n12= count12{ia}{j};
%                n1= count1{ia}{j};
%                n2= count2{ia}{j};
                N= sum(nRips{ia}{j}(iepoch{ia}{j},:),2);
                n12= sum(count12{ia}{j},2);
                n1= sum(count1{ia}{j},2);
                n2= sum(count2{ia}{j},2);

                z{ia}{j}= nan*ones(size(N));
%                N(N<minrips)= nan; 
                N(N<2)= nan;  %%@@
                V= n1.*n2.*(N-n1).*(N-n2)./ (N.^2.*(N-1));
                V(V<=0)= nan;
                z{ia}{j}= (n12-n1.*n2./N)./sqrt(V);

                ind= isfinite(z{ia}{j});
                Z{ia}= [Z{ia}; z{ia}{j}(ind)];
                IP{is}{ia}= [ IP{is}{ia}, aux{ia}.IP(ind)];
            end
            Zd{is}{ia}= Z{ia};
        end
        opt.label= 'z'; opt.title= 'z';
        sel{1}.label= 'fam'; sel{2}.label= 'novel';
        if plotZ auxCmpDist2(Zd{is}, sel, opt); end
%        opt.yrange= [0 1];
%        auxImage(b, z, opt)
%        auxDist(b, z, sets, opt)
%        auxAnova1(Zd{is}, sets{is}, opt);
%        opt.plotcol(3,:)= 0.7*[1 1 1];
%        auxCmpDist2(Zd{is}, sel, opt);
        save('data_Zd', 'Zd', 'IP');
    end
        
    if plotDN12
        opt.label= 'n12-E[n12]'; opt.title= 'dn12';
%        opt.yrange= [0 1];
        auxImage(b, dn12, opt)
    end
        
    if plotNPairs
        opt.label= 'n pairs'; opt.title= 'nPairs';
        auxImage(b, nPairs, opt)
    end
        
    if 1 & plotRipRate %%@@
        for ia=1:na
            r= []; IC{is}{ia}= [];
            for j=2
                N= sum(nRips{ia}{j}, 2);
                O= sum(occ{ia}{j}, 2);
                ind= O>= 1;
                r= [r; N(ind)./ O(ind)];
                IC{is}{ia}= [ IC{is}{ia}, aux{ia}.IP(ind)];
            end
            RIPRATE{is}{ia}= r;
        end
        save('data_RIPRATE', 'RIPRATE', 'IC');
    end
    if plotRipRate
        ripRate= {};
        for ia=1:na
            for j=1:2
                N= nRips{ia}{j};
                O= occ{ia}{j};
                O(O<= 1)= nan;
                ripRate{ia}{j}= N./O;
            end
        end
        opt.label= 'ripple rate (Hz)'; opt.title= 'ripRate'; 
%        opt.yrange= [0 0.6];
        opt.yrange= [0 0.1];
        auxImage(b, ripRate, opt)
    end

    if plotRipLength
        for ia=1:na
            r= [];
            for j=2
                for ie=1:length(rip{ia}{j}.length)
                    r= [r; 1000*rip{ia}{j}.length{ie}];
                end
            end
            RIPLENGTH{is}{ia}= r;
        end
        opt.label= 'ripple length (ms)'; opt.title= 'ripLength';
        opt.plot='pdf'; opt.nbins= 20;
        auxCmpDist2(RIPLENGTH{is}, sel, opt);
    end

    if plotNRips
        for ia=1:na
            r= [];
            for j=2
                r= [r; sum(nRips{ia}{j}, 2)];
            end
            NRIPS{is}{ia}= r;
        end
%        opt.label= 'n rips';
%        auxCmpDist2(NRIPS{is}, sel, opt);
    end

    if plotRipFrac
        for ia=1:na
            r= [];
            for j=2
                ltmp= [];
                for ie=1:length(rip{ia}{j}.length)
                    ltmp(ie,:)= sum(rip{ia}{j}.length{ie});
                end
                O= sum(occ{ia}{j}, 2);
                r= [r; 100*ltmp./O];
            end
            RIPFRAC{is}{ia}= r;
        end
        opt.label= 'ripple frac (%)'; opt.title= 'ripFrac';
        opt.plot='pdf'; opt.nbins= 20;
        auxCmpDist2(RIPFRAC{is}, sel, opt);
    end

    if plotActivation
        for ia=1:na
            r= []; IC{is}{ia}= [];
            for j=2
                n1= sum(count1{ia}{j}, 2);
                N= sum(nRips{ia}{j}(iepoch{ia}{j},:), 2);
                ind= N>= minrips;

%                cl= collectCellList([id '-' armId{ia}]);
%                ind= ind & strcmp({cl.rat{aux{ia}.IP}}, 'kyl')';
%                fprintf(1, '%d, %d, %d \n', length(N), length(aux{ia}.IP), sum(ind));

                r= [r; n1(ind)./N(ind)];
                IC{is}{ia}= [ IC{is}{ia}, aux{ia}.IP(ind)];
            end
            ACTIVATION{is}{ia}= r;
        end
        opt.label= 'activation prob'; opt.title= 'activation';
        save('data_ACTIVATION', 'ACTIVATION', 'IC');
%        auxCmpDist2(ACTIVATION{is}, sel, opt);
    end

    if plotBurst_active
        for ia=1:na
            r= [];
            for j=2
                n1= sum(count1{ia}{j}, 2);
                b= sum(nBursts{ia}{j}, 2);
                ind= n1>= opt.nmin;
                r= [r; b(ind)./n1(ind)];
            end
            BURST_ACTIVE{is}{ia}= r;
        end
        opt.label= 'burst/activation'; opt.title= 'burst_active';
        auxCmpDist2(BURST_ACTIVE{is}, sel, opt);
    end

    if plotSpikes_active
        for ia=1:na
            r= [];
            for j=2
                n1= sum(count1{ia}{j}, 2);
                sp= sum(nSpikes{ia}{j}, 2);
                ind= n1>= opt.nmin;
                r= [r; sp(ind)./n1(ind)];
            end
            SPIKES_ACTIVE{is}{ia}= r;
        end
        opt.label= 'spikes/activation'; opt.title= 'spikes_active';
%        auxCmpDist2(SPIKES_ACTIVE{is}, sel, opt);
    end

    if plotRipSpikeRatio
        for ia=1:na
            r= [];
            for j=2
                sp= sum(nSpikes{ia}{j}, 2);
                n= sum(allSpikes{ia}{j}, 2);
                ind= n>= opt.nmin;
                r= [r; sp(ind)./n(ind)];
            end
            RIPSPIKERATIO{is}{ia}= r;
        end
        opt.label= 'HFE spikes/ all spikes'; opt.title= 'ripSpikeRatio';
        auxCmpDist2(RIPSPIKERATIO{is}, sel, opt);
    end

    if plotAllSpikes
        for ia=1:na
            r= [];
            for j=2
                n= sum(allSpikes{ia}{j}, 2);
                r= [r; n];
            end
            ALLSPIKES{is}{ia}= r;
        end
%        opt.label= 'all spikes'; opt.title= 'allSpikes';
%        auxImage(b, allSpikes, opt)
%        auxCmpDist2(ALLSPIKES{is}, sel, opt);
    end

    if plotSpikes_Burst
        opt.label= 'spikes/burst'; opt.title= 'ripSpikes_Burst';
        opt.yrange= [0.98 3];
        for ia=1:na
            r= [];
            for j=2
                n= sum(nSpikes{ia}{j},2);
                N= sum(nBursts{ia}{j},2);
                ind= N>= opt.nmin;
                r= [r; n(ind)./ N(ind)];
            end
            sp_burst{is}{ia}= r;
        end
%        auxImage(b, sp_burst, opt)
    end
    if plotSpikeRate
        opt.label= 'spikes/ripple'; opt.title= 'ripSpikeRate';
%        opt.yrange= [-0.5 35];
        for ia=1:na
            sp= [];
            for j=2
                N= sum(nRips{ia}{j}(iepoch{ia}{j},:), 2);
                n= sum(nSpikes{ia}{j}, 2);
                spikerate{ia}{j}= nan*ones(size(N));
                ind= N>= opt.nmin;
                spikerate{ia}{j}(ind)= n(ind)./ N(ind);
                sp= [sp; n(ind)./ N(ind)];
            end
            SPRATE{is}{ia}= sp;
        end
%        auxImage(b, spikerate, opt)
%        auxAnova1(SPRATE{is}, sets{is}, opt);
%        auxCmpDist2(SPRATE{is}, sel, opt);
    end
    if plotBurstRate
        opt.label= 'burst/ripple'; opt.title= 'ripBurstRate';
%        opt.yrange= [-0.25 25];
        burstrate= {};
        for ia=1:na
            for j=2
                N= nRips{ia}{j}(iepoch{ia}{j},:);
                n= nBursts{ia}{j};
                burstrate{ia}{j}= nan*ones(size(N));
                ind= N>= opt.nmin;
                burstrate{ia}{j}(ind)= n(ind)./ N(ind);
            end
        end
        auxImage(b, burstrate, opt)
    end

    if plotOcc
        for ia=1:na
            r= [];
            for j=2
                r= [r; sum(occ{ia}{j}, 2)];
            end
            OCC{is}{ia}= r;
        end
    end


%    if plotOcc
%        opt.label= 'occupancy (s)'; opt.title= 'occ';
%        opt.yrange= [0 130];
%        auxImage(b, occ, opt)
%    end
    if plotN12
        opt.label= 'coinc. count'; opt.title= 'n12';
%        opt.yrange= [0 130];
        auxImage(b, count12, opt)
    end
    if plotVel
        for ia=1:na
            v= [];
            for j=2
                v= [v, velocity{ia}{j}];
            end
            VEL{is}{ia}= v;
        end
        opt.label= 'velocity'; opt.title= 'vel';
        opt.plot= 'hist';
        figure
        auxCmpDist2(VEL{is}, sel, opt);
    end
    if plotRipVel
        RIPVEL{is}{ia}= [];
        for ia=1:na
            RIPVEL{is}{ia}= [RIPVEL{is}{ia} ripVel{ia}{2}(isfinite(ripVel{ia}{2}))];
        end
        opt.label= 'velocity during HFE'; opt.title= 'ripVel';
        opt.plot= 'hist'; opt.nbins= 25; opt.xrange=[0, 50];
        figure
        auxCmpDist2(RIPVEL{is}, sel, opt);
    end
end % for is; show analysis


range={[1:3], [4:6]};
for i=1:round(length(sets)/3)
    if pairs>=0
        opt.outname= sprintf('%s_%s%s_%s', ripStr, pairStr, sets{(i-1)*3+1}.arm, id);
    else
        opt.outname= sprintf('%s_%s%s', ripStr, pairStr, id);
    end
    if plotZd
        opt.ylabel= 'z-score'; opt.title= 'z'; 
        opt.plotsize=[1.5 2/3];
    %    opt.ylim= [-0.1 1.9];
        auxPlotBars(Zd(range{i}), opt);
%        auxAnova(Zd(range{i}), sets(range{i}), opt);
    end
    if plotCoincFrac
        opt.ylabel= 'coincidence frac'; opt.title= 'ripCoincFrac';
        opt.ylim= [0 0.3];
        auxPlotBars(COINCFRAC(range{i}), opt);
%        auxAnova(COINCFRAC(range{i}), sets(range{i}), opt);
    end

    if plotNRips
        for ir=range{i}
            fprintf(1, '%d\t%d\n', sum(NRIPS{ir}{1}), sum(NRIPS{ir}{2}));
        end
%        opt.ylabel= 'ripple count'; opt.title= 'nRips'; 
%        auxPlotBars(NRIPS(range{i}), opt);
    end

    if plotOcc
        opt.ylabel= 'occupancy (s)'; opt.title= 'occ'; 
        auxPlotBars(OCC(range{i}), opt);
    end

    if plotRipRate
        opt.ylabel= 'ripple rate (Hz)'; opt.title= 'ripRate'; 
    %    opt.yrange= [-0.05 0.9];
        auxPlotBars(RIPRATE(range{i}), opt);
%        auxAnova(RIPRATE(range{i}), sets(range{i}), opt);
    end

    if plotRipLength
        opt.ylabel= 'ripple length (ms)'; opt.title= 'ripLength'; 
        auxPlotBars(RIPLENGTH(range{i}), opt);
    end

    if plotRipFrac
        opt.ylabel= 'ripple frac (%)'; opt.title= 'ripFrac'; 
        auxPlotBars(RIPFRAC(range{i}), opt);
    end

    if plotSpikeRate
        opt.ylabel= 'spikes/ripple'; opt.title= 'ripSpikeRate'; 
        opt.plotsize=[1.5 2/3]; 
        opt.ylim= [0 1.3];
        auxPlotBars(SPRATE(range{i}), opt);
%        auxAnova(SPRATE(range{i}), sets(range{i}), opt);
    end
    if plotVel
        opt.ylabel= 'velocity (cm/s)'; opt.title= 'vel'; 
        opt.plotsize=[1.5 2/3]; 
    %    keyboard
        auxPlotBars(VEL(range{i}), opt);
    end
    if plotRipVel
        opt.ylabel= 'HFE velocity (cm/s)'; opt.title= 'ripVel'; 
        opt.plotsize=[1.5 2/3]; 
    %    keyboard
        auxPlotBars(RIPVEL(range{i}), opt);
    end
    if plotActivation
        opt.ylabel= 'activation prob'; opt.title= 'activation'; 
        opt.plotsize=[1.5 2/3]; opt.ylim= [0 0.6];
        auxPlotBars(ACTIVATION(range{i}), opt);
    end
    if plotBurst_active
        opt.ylabel= 'bursts/activation'; opt.title= 'burst_active'; 
        opt.plotsize=[1.5 2/3]; 
        auxPlotBars(BURST_ACTIVE(range{i}), opt);
    end
    if plotSpikes_active
        opt.ylabel= 'spikes/activation'; opt.title= 'spikes_active'; 
        opt.plotsize=[1.5 2/3]; 
        auxPlotBars(SPIKES_ACTIVE(range{i}), opt);
    end
    if plotRipSpikeRatio
        opt.ylabel= 'HFE spikes/ all spikes'; opt.title= 'ripSpikeRatio'; 
        opt.plotsize=[1.5 2/3]; opt.ylim=[0 .4];
        auxPlotBars(RIPSPIKERATIO(range{i}), opt);
    end
    if plotAllSpikes
%        for ir=range{i}
%            fprintf(1, '%d\t%d\n', sum(ALLSPIKES{ir}{1}), sum(ALLSPIKES{ir}{2}));
%        end
        opt.ylabel= 'all spikes'; opt.title= 'allSpikes'; 
        opt.plotsize=[1.5 2/3]; 
        auxPlotBars(ALLSPIKES(range{i}), opt);
    end

    if plotSpikes_Burst
        opt.ylabel= 'spikes/burst'; opt.title= 'spikes_burst'; 
        opt.plotsize=[1.5 2/3]; 
        auxPlotBars(sp_burst(range{i}), opt);
    end
end

if plotZd
    opt.ylabel= 'z-score'; opt.title= 'z'; 
%    auxAnova(Zd, sets, opt);
end

if plotCoincFrac
    opt.ylabel= 'coincidence frac'; opt.title= 'ripCoincFrac';
%    auxAnova(COINCFRAC, sets, opt);
end

if plotSpikeRate
    opt.ylabel= 'spikes/ripple'; opt.title= 'ripSpikeRate'; 
%    auxAnova(SPRATE, sets, opt);
end

if plotRipRate
    opt.ylabel= 'ripple rate (Hz)'; opt.title= 'ripRate'; 
%    auxAnova(RIPRATE, sets, opt);
end

function auxImage(b, z, opt, dev)

% calc mean
for ia=1:length(z)
    for j=1:length(z{ia})
        zm{ia}{j}= mean(z{ia}{j});
    end
end

if isfield(opt, 'yrange')
    minz= opt.yrange(1); maxz= opt.yrange(2);
else
    % find max value
    minz= 0; maxz= nan; 
    for ia=1:length(z)
        for j=1:length(z{ia})
            maxz= max(maxz, max(max(zm{ia}{j})));
        end
    end
end

for j=1:2
    ib{j}= 1:length(b{j});
end

thick= 2;

if(0)
figure
hold on
%novel home arm
imagesc(1:thick, -fliplr(ib{1}), zm{1}{1}'*ones(1,thick), [minz maxz])
%novel arm
imagesc(ib{2}+thick, 0:thick-1, ones(thick,1)*zm{1}{2}, [minz maxz])

%fam home arm
imagesc(-[1:thick], -fliplr(ib{1}), zm{2}{1}'*ones(1,thick), [minz maxz])
%fam arm
imagesc(-ib{2}-thick, 0:thick-1, ones(thick,1)*zm{2}{2}, [minz maxz])
axis off auto equal 
%colorbar
figname= ['geo_' opt.title '_' opt.outname];
set(gcf, 'Name', figname);
myprint(2*[1 2/3], figname);
end


x= [b{1}(1:end-1) b{2}(1:end-1)]+75;
dx= mean(diff(x));
x= x+dx/2;
%      fam arm           novel arm
y= {[z{1}{1} z{1}{2}], [z{2}{1} z{2}{2}]};
%novel arm
%plot([b{1}(1:end-1) b{2}(1:end-1)], [zm{1}{1} zm{1}{2}], 'r')
%fam arm
%plot([b{1}(1:end-1) b{2}(1:end-1)], [zm{2}{1} zm{2}{2}], 'k')

if(1)
figure
hold on
if nargin>= 4;
    d= {[dev{1}{1} dev{1}{2}], [dev{2}{1} dev{2}{2}]};
    h= auxPlotMeanDev(x,y,opt,d);
else
    h= auxPlotMeanDev(x,y,opt);
end
hold on
plot([75 75], get(gca, 'YLim'), 'k--');
%set(gca, 'xlim', [-75, 75]);
%legend({'novel', 'fam'}, 0)
xlabel('distance (cm)');
ylabel(opt.label);
figname= ['lin_' opt.title '_' opt.outname];
set(gcf, 'Name', figname);
myprint([1.5 1], figname);
end

if(0)
opt.ylabel= opt.label;
opt.plotsize=[1.5 2/3];
    opt.outname= ['lin_' opt.title '_' opt.outname];
    for i=1:2; 
        my(:,i)= nanmean(y{i})'; 
        n= sum(isfinite(y{i}));
        n(n<1)= nan;
        dy(:,i)= nanstd(y{i})'./sqrt(n)'; 
        my(~isfinite(dy(:,i)),i)= nan;
    end
    auxBar(my, dy, opt);
end

function h= saveHistc(x, b)
if isempty(x);
    h= zeros(1, length(b));
else
    h= histc(x, b);
end
if size(h,1)>1 & size(h,2)==1; h= h'; end

function auxDist(b, z, sets, opt)
x= [b{1}(1:end-1) b{2}(1:end-1)]+75;
dx= mean(diff(x));
x= x+dx/2;
%      fam arm           novel arm
y= {[z{1}{1} z{1}{2}], [z{2}{1} z{2}{2}]};
%figure
%subplot(2,1,1)
%plot(x,y{1}')
%title('fam')
%subplot(2,1,2)
%plot(x,y{2}')
%title('novel')
%myprint([1.5 1], figname);
%xlabel('position (cm)');
%ylabel('z score');
%orient tall
%print('-dpsc', figname);

hold on
for i=1:length(x)
    for j=1:2
        t= y{j}(:,i);
        t= t(isfinite(t));
        if length(t)>2 & sum(t) > 0
            [h, p]= ttest(t,0.5);
            if(p<0.01) 
                hs= plot(x(i)+[-1 1], mean(t)+.2*[1 1], '*'); 
                set(hs, 'Color', opt.plotcol(j,:));
            elseif(p<0.05) 
                hs= plot(x(i), mean(t)+.1, '*'); 
                set(hs, 'Color', opt.plotcol(j,:));
            end
        end
    end
end


%figname= ['lin-' opt.title '-' opt.outname];
%set(gcf, 'Name', figname);
%myprint('large', figname, [], 0);

function y= auxCollate(z)
%      fam arm           novel arm
y= {[z{1}{1} z{1}{2}], [z{2}{1} z{2}{2}]};

function auxPlotBars(Zd, opt)
for d=1:3
    for ia=1:length(Zd{d})
        Zd{d}{ia}= Zd{d}{ia}(isfinite(Zd{d}{ia}));
    end
end

for d=1:3
opt.plotcol={zeros(3,3), [[1 0 0]; [0 1 0]; [0 0 1]], [.7*ones(3,3)]};
    for ia=1:length(Zd{d})
        ja= ia;

% comparing novel vs mixed pairs
%opt.plotcol={[[1 0 0]; [0 1 0]; [0 0 1]], [.7*ones(3,3)]};
%    for ia=2:length(Zd{d})
%        ja= ia-1;

        Zm(d,ja)= nanmean(Zd{d}{ia});
        Zdev(d,ja)= nanstd(Zd{d}{ia})/sqrt(sum(isfinite(Zd{d}{ia})));
        if ja>1
%             [h,p(d,ia-1)]= ttest2(Zd{d}{ia-1}, Zd{d}{ia});
             [p(d,ja-1), h]= ranksum(Zd{d}{ia-1}, Zd{d}{ia});
%             [h,p(d,ia-1)]= kstest2(Zd{d}{ia-1}, Zd{d}{ia});
         end
    end
%    p(d,:)= p(d,:)* length(p(d,:));
%    opt.outname
end
opt.xticklabels= {'day1'  'day2'  'day3'}; 
auxBar(Zm, Zdev, opt, p);
%auxBar(Zm, Zdev, opt);
p

function auxAnova1(Z, set, opt)
return
figure
n= 0; g=[]; y=[];
for ia=1:length(Z)
    len= length(Z{ia});
    y(n+1:n+len)= Z{ia};
    g(n+1:n+len)= ia;
    n=n+len;
end
%[p, tab, stats]= anova1(y, g, 'off');
[p, tab, stats]= kruskalwallis(y, g, 'off');
%mc= multcompare(stats, 'alpha', 0.01);
mc= multcompare(stats, 'alpha', 0.01, 'ctype', 'bonferroni'); tab
%set(gcf, 'Name', [opt.label '-' set.arm])
%mc
%keyboard

function auxAnova(Z, sets, opt)

figure
ns= length(sets);
for is=1:ns
    switch sets{is}.arm
    case 'famArm';
        sets{is}.arm= 1;
    case 'novelArm';
        sets{is}.arm= 2;
    case 'cross';
        sets{is}.arm= 3;
    end
end
global armId
if ns>3
    names= strvcat('pf-loc', 'day', 'arm');
else
    names= strvcat('pf-loc', 'day');
end
n= 0; g={}; y=[];
for is=1:length(Z)
    for ia=1:length(Z{is})
%        if(sets{is}.arm==1) continue; end
        len= length(Z{is}{ia});
        y(n+1:n+len)= Z{is}{ia};
        g{1}(n+1:n+len)= ia;
        g{2}(n+1:n+len)= sets{is}.novelDay;
        if ns>3;  g{3}(n+1:n+len)= sets{is}.arm; end
        n=n+len;
    end
end

[p, tab, stats]= anovan(y, g, 'full', 3, names, 'off');
%mc= multcompare(stats, 'alpha', 0.01, 'ctype', 'bonferroni', 'dimension', [1:length(g)]);
mc= multcompare(stats, 'alpha', 0.05, 'dimension', [1:length(g)]);
tab
%mc
%keyboard
