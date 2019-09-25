function show2dHist(plotFigs, bursts, num)
% function show2dHist(plotFigs, bursts, num)

global fmaux
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
%        pause
    end
else
    d=num(1); e=num(2); t=num(3); c=num(4);
    fmaux.currentCell= getSelectId(num);
    auxrun(d,e,t,c,opt);
end

function auxrun(d,e,t,c,opt)

if opt.plotFigs
    figure(1); clf;
    figure(2); clf;
    figure(3); clf;
end
global select fmaux
data=loadData([d e t c]);
small=1e-4;
minpos=min(data.linpos)-small; % around 6cm
maxpos=max(data.linpos)+small; % around 160cm
                                                                                
%cpx= [linspace(minpos,maxpos,151)];
%cpy= [linspace(-small,2*pi+small,31)];
cpx= [linspace(minpos,maxpos,31)];
cpy= [linspace(-small,2*pi+small,21)];

for traj=1:4
    cumspikes{traj}= zeros(length(cpy)-1,length(cpx)-1);
end

cmap1max=100;
cmap1= jet(cmap1max);
cmap1=[[0,0,0]; cmap1];
cmap2max= 10;
cmap2= jet(cmap2max);
cmap2=[[0,0,0]; cmap2];

dt= data.time(2)-data.time(1);

n= fmaux.currentCell;
A= select.a{n};
%for ia=1:length(A)
%    traj= A{ia}.traj;
for ia=1:4
    traj= ia-1;
    itraj= find(data.traj== traj)';
    isptraj= find(data.traj(data.spikeindex)== traj);

    b=[];
    b(:,1)=data.linpos(itraj);
    b(:,2)=data.phase(itraj);
    [occall,bx,by]= hist2(b,cpx,cpy);
    if opt.bursts
        isptraj= isptraj(findBursts(data.spiketimes(isptraj)));
    end
    sp=[];
    sp(:,1)=data.linpos(data.spikeindex(isptraj));
    sp(:,2)=data.phase(data.spikeindex(isptraj));
    [nspikesall,bx,by]= hist2(sp,cpx,cpy);
    [ibad, jbad]= find(occall==0);
    if ~isempty(ibad)
        occ(sub2ind(size(occall), ibad, jbad))= nan;
    end
    rateall= nspikesall./(occall*dt);
%    rateall(sub2ind(size(rateall), ibad, jbad))= -0.5;


    passt= [data.time(1), select.x{n}{ia}.time, data.time(end)];
    npass= length(passt)-1;
    for ipass= 1:npass
        %  select valid timesteps
        itime= itraj(find(data.time(itraj) > passt(ipass) & data.time(itraj) < passt(ipass+1)));
        isptime= isptraj(find(data.spiketimes(isptraj) > passt(ipass) & data.spiketimes(isptraj) < passt(ipass+1)));

        b=[];
        b(:,1)=data.linpos(itime);
        b(:,2)=data.phase(itime);
        [occ,bx,by]= hist2(b,cpx,cpy);

        if opt.bursts
            isptime= isptime(findBursts(data.spiketimes(isptime)));
            fprintf('%d bursts on traj %d\n', length(isptime), traj);
        else
            fprintf('%d spikes on traj %d\n', length(isptime), traj);
        end
        sp=[];
        sp(:,1)=data.linpos(data.spikeindex(isptime));
        sp(:,2)=data.phase(data.spikeindex(isptime));
        [nspikes,bx,by]= hist2(sp,cpx,cpy);
        cumspikes{traj+1}= cumspikes{traj+1}+nspikes;

%        [ibad, jbad]= find(occ==0);
    %    occ(sub2ind(size(occ), ibad, jbad))= -0.5;

    %    rate= nspikes./(occ*dt);
    %    rate(sub2ind(size(rate), ibad, jbad))= -0.5;

    %    [izero, jzero]= find(nspikes<=0);
    %    nspikes(sub2ind(size(nspikes), izero, jzero))= -0.5;
    %    rate(sub2ind(size(rate), izero, jzero))= -0.5;
    %    keyboard
        
        if opt.plotFigs
%            subplot(4,3,3*traj+1);
%            figure(1)
%            colormap(cmap1);
%            imagesc(bx,by,occ, [0 cmap1max]);
%            plot(b(:,1),b(:,2),'.','MarkerSize',1);
%            axis xy; colorbar
%            if traj==0; title('occupancy'); end

%            figure(2)
            colormap(cmap2);
%            subplot(4,2,2*traj+1);
%            subplot(2,1,1);
%            imagesc(bx,by,nspikes, [0 cmap2max]); 
%            axis xy; colorbar
%            if traj==0; title(sprintf('spike count, pass %d/%d',ipass,npass)); end
%            hold on
%            plot(b(:,1),b(:,2),'.','MarkerSize',1,'MarkerEdgeColor',[.2 .2 .2]);

%            subplot(4,2,2*traj+2);
            subplot(4,1,traj+1);
%            subplot(2,1,2);
%            imagesc(bx,by,cumspikes{traj+1}, [0 cmap2max]); 
            imagesc(bx,by,nspikesall, [0 max(max(nspikesall))]); 
%            imagesc(bx,by,rateall, [0 max(max(rateall))]);
            axis xy; colorbar
            if traj==0; title('cum spike count'); end
            hold on
            plot(b(:,1),b(:,2),'.','MarkerSize',1,'MarkerEdgeColor',[.3 .3 .3]);
            plot(sp(:,1),sp(:,2),'+','MarkerSize',4,'MarkerEdgeColor','w');
            key= input('(s)kip to next traj', 's');
            if key=='s'; break; end
        else
%            keyboard
            cumspikes{traj+1}
            pause
        end
    end
    disp('next traj'); pause
end
