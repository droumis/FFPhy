function show2dHist(plotFigs, bursts, num)

opt.plotFigs= plotFigs;

if nargin < 2
    opt.bursts= 0;
else
    opt.bursts= bursts;
end

if nargin < 3
    [d,e,t,c,ncells]= startCellList;
    while ~isempty(d)
        auxrun(d,e,t,c,opt);
        [d,e,t,c]= getNextCell;
        pause
    end
else
    d=num(1); e=num(2); t=num(3); c=num(4);
    auxrun(d,e,t,c,opt);
end

function auxrun(d,e,t,c,opt)

data=loadData([d e t c]);
small=1e-4;
minpos=min(data.linpos)-small; % around 6cm
maxpos=max(data.linpos)+small; % around 160cm
                                                                                
%cpx= [linspace(minpos,maxpos,51)];
%cpy= [linspace(-small,2*pi+small,31)];
cpx= [linspace(minpos,maxpos,31)];
cpy= [linspace(-small,2*pi+small,21)];

dt= data.time(2)-data.time(1);
for traj=0:3
    itraj= find(data.traj== traj)';
    b=[];
    b(:,1)=data.linpos(itraj);
    b(:,2)=data.phase(itraj);
    [occ,bx,by]= hist2(b,cpx,cpy);

    spikeitraj= find(data.traj(data.spikeindex)== traj);
    if opt.bursts
        spikeitraj= spikeitraj(findBursts(data.spiketimes(spikeitraj)));
        fprintf('%d bursts on traj %d\n', length(spikeitraj), traj);
    else
        fprintf('%d spikes on traj %d\n', length(spikeitraj), traj);
    end
    sp=[];
    if ~isempty(spikeitraj)
        sp(:,1)=data.linpos(data.spikeindex(spikeitraj));
        sp(:,2)=data.phase(data.spikeindex(spikeitraj));
        [nspikes,bx,by]= hist2(sp,cpx,cpy);
    else
        nspikes= zeros(length(cpy)-1, length(cpx)-1);
    end

    F= getFanoFactor(nspikes);
    [yhot, xhot]= find(F>2.5);
    F1d= sort(reshape(F, prod(size(F)), 1)); 

    [ibad, jbad]= find(occ==0);
    if ~isempty(ibad)
        occ(sub2ind(size(occ), ibad, jbad))= -0.5;
    end

    rate= nspikes./(occ*dt);
    if ~isempty(ibad)
        rate(sub2ind(size(rate), ibad, jbad))= -0.5;
    %    rate(sub2ind(size(rate), ibad, jbad))= 0;
    end

    [izero, jzero]= find(nspikes<=0);
%    nspikes(sub2ind(size(nspikes), izero, jzero))= -0.5;
%    rate(sub2ind(size(rate), izero, jzero))= -0.5;
%    keyboard
    cmap= colormap;
    cmap(1,:)=[0,0,0];
    colormap(cmap);
    
    if opt.plotFigs
%    figure(1)
        subplot(4,3,3*traj+1);
%        low= min(min(occ));
%        if(low < 10)
%            high= 30;
%        else
%            high= 3*low;
%        end
%        imagesc(bx,by,occ, [-0.5,high]); 
%        axis xy; colorbar
%        if traj==0; title('occupancy'); end

        plot(F1d, [1:length(F1d)]/length(F1d));
        axis tight
        if traj==0; title('Fano Factors'); end

        %    figure(2)
        subplot(4,3,3*traj+2);
        high= max(max(nspikes));
        if high<=0; high= 1; end
        imagesc(bx,by,nspikes, [0, high]);
        axis xy; colorbar
        hold on
        plot(bx(xhot),by(yhot), 'w*', 'MarkerSize', 2);
        if traj==0; title('spike count'); end

        %    figure(3)
        subplot(4,3,3*traj+3);
        high= max(max(rate));
        if high<=0; high= 1; end
        imagesc(bx,by,rate, [0, high]);
        axis xy; colorbar
        hold on
        plot(bx(xhot),by(yhot), 'w*', 'MarkerSize', 2);
        if traj==0; title('rate'); end
    else
%        occ
%        nspikes
        [bx(xhot)',by(yhot)',  F(sub2ind(size(nspikes), yhot, xhot))]
        [hF bF]= hist(F1d)
        length(xhot)/length(F1d)
%        sum(nspikes)
%        rate
%        mean(rate)
%        keyboard
        pause
    end
end
