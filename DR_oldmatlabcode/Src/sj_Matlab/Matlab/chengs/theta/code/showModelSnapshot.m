function showModelSnapshot(num, traj, nt, scale, vert)
%function showModelSnapshot(num, [traj, nt, scale, vert])
%
% scale=0:    use same colorscale for all plots

if nargin < 2; traj= [0 1 2 3]; end
if nargin < 3; nt= 4; end
if nargin < 4; scale= 0; end
if nargin < 5; vert= 'time'; end

if nargin < 1 | isempty(num)
    [d,e,t,c,ncells]= startCellList;
    while ~isempty(d)
        auxrun([d,e,t,c], traj, nt, scale, vert);
        [d,e,t,c]= getNextCell;
    end
else
    auxrun(num, traj, nt, scale, vert);
end

function auxrun(num, traj, nt, scale, vert)

global adaptest fmaux select
data=loadData(num);
loadVar('.','adaptest',num(1));
m= adaptest{num(1)}{num(2)}{num(3)}{num(4)}.model;
opt.traj= traj;
mini= min(find(data.traj>=0))+1000;
maxi= max(find(data.traj>=0))-1000;
opt.t= linspace(data.time(mini),data.time(maxi), 4);
opt.nx= 200; opt.ny= 70; % small print quality
%opt.nx= 400; opt.ny= 140; % large print quality

out= eval2DModel(data, m, opt);
t= out.t; x= out.x; y= out.y; z= out.z;
%[xp, isi]= auxLoadModels(num);
%[t, x,y,z]= auxEvalModel(xp, opt.traj, opt.nt, opt.nx, opt.ny);

% load and set up selection criteria
load(fmaux.select);
selid= getSelectId(num);
nsel= length(select.a{selid});
seltraj=[];
for is=1:nsel;  seltraj(end+1)= select.a{selid}{is}.traj; end
offset= select.a{selid}{1}.offset
[tmp,cut]= min(mod(y-offset,2*pi));

ntraj= size(z,1);
% find max
if scale
    zmax= 0;
    for itraj=1:ntraj
        for it=1:nt
            tmp= max(max(z{itraj,it}));
            if tmp> zmax; zmax= tmp; end 
        end
    end
end

%save aa %@@

fh= figure;
set(fh, 'Name', sprintf('[%d %d %d %d], %s', num, pwd));
set(fh, 'Position', [82 161 1229 785]);
set(0, 'DefaultAxesFontSize', 8);
for itraj=1:ntraj
    for it=1:nt
        if strcmp(vert, 'traj')
            subplot(ntraj, nt, (itraj-1)*nt+it);
        else
            subplot(nt, ntraj, (it-1)*ntraj+itraj);
        end

        %@@
%        [i,j]= find(z{itraj,it} < 6); 
%        z{itraj,it}(sub2ind(size(z{itraj,it}), i, j))= 0;

        ztmp= z{itraj,it}';
        ytmp= y;
        ztmp= [ztmp(cut:end,:); ztmp(1:cut-1,:)];
%        keyboard
%        ztmp= [ztmp; ztmp];
%        ytmp= [y, 2*y];
        if scale
            imagesc(x,ytmp,ztmp, [0, zmax]);
        else
            imagesc(x,ytmp,ztmp);
        end
        axis xy
        colorbar

        % label plots
        if strcmp(vert, 'traj')
            if (it==1) 
                ylabel('phase (rad)')
            end
            if (itraj==1) auxAppendTitle(sprintf('%.0f sec', t(it))); end
            if(it==nt) auxAppendTitle(sprintf('traj %d', traj(itraj))); end
            if (itraj==ntraj) xlabel('position (cm)'); end
        else
            if (itraj==1) ylabel('phase (rad)'); end
            if (it==1) auxAppendTitle(sprintf('traj %d', traj(itraj))); end
            if(itraj==1) auxAppendTitle(sprintf('%.0f sec', t(it))); end
            if (it==nt) xlabel('position (cm)'); end
        end


        matchid= find(traj(itraj)== seltraj);
%        keyboard

        for im=1:length(matchid)
            minpos= select.a{selid}{matchid(im)}.linpos(1);
            maxpos= select.a{selid}{matchid(im)}.linpos(2);
            % show selection
            line([[minpos; minpos], [maxpos;maxpos]], [[y(1);2*y(end)], [y(1);2*y(end)]], 'color', [1,1,1]);
            % show mutual info on middle of HomeArm
            ind= find(minpos < x & x < maxpos);
            [mx,my,Vx,Vy,R,I]= stats2d(z{itraj,it}(ind,:),x(ind),y);
            auxAppendTitle(sprintf('I=%.3f', I));
        end
    end
end

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 3 4]); % [left bottom width height]
figname= sprintf('snapshots_%d_%d_%d_%d', num);
print('-dpng', figname)
saveas(gcf, figname)
%print('-dpng', sprintf('snapshots_%d_%d_%d_%d', num))
%print('-dtiff', sprintf('snapshots_%d_%d_%d_%d', num))
%print('-dill', sprintf('snapshots_%d_%d_%d_%d', num))
%keyboard

function auxAppendTitle(str)
tstr= get(get(gca, 'title'), 'string');
if ~isempty(tstr); tstr= [tstr ', ']; end
title([tstr str]);
