function calcDist(num)
%function calcDist(num)
%  

global fmaux

fname= [fmaux.data2dir '/info.mat'];
if exist(fname, 'file');
    load(fname);
else
    info= [];
end

oldd= -1;
olde= -1;
if nargin < 1
    [d,e,t,c]= startCellList;
    while ~isempty(d)
        if oldd~= d | olde~= e
            info= auxrun(info,d,e);
            oldd= d; olde= e;
        end
        [d,e,t,c]= getNextCell;
    end
else
    fmaux.currentCell= getSelectId(num);
    auxrun(num(1),num(2));
end

save(fname, 'info');

function info= auxrun(info,d,e)

global fmaux task behavdata
loadVar(fmaux.datadir, 'task', 0, fmaux.prefix, 1);
loadVar(fmaux.data2dir, 'behavdata', d);

z= task{d}{e}.centercoord*task{d}{e}.pixelsize;

n= length(behavdata{d}{e}.xpos);
[tmp, i]= min(sum(([behavdata{d}{e}.xpos, behavdata{d}{e}.ypos]-ones(n,1)*z).^2, 2));
info{d}{e}.centerlinpos= behavdata{d}{e}.linpos(i);

if sqrt(tmp)>0.5; error('Could not find close enough a position to center'); end
info{d}{e}.minlinpos= min(behavdata{d}{e}.linpos);
info{d}{e}.maxlinpos= max(behavdata{d}{e}.linpos);
fprintf(1, '[%d %d] center=%.2f, min=%.2f, max=%.2f\n', d, e, ...
    info{d}{e}.centerlinpos, info{d}{e}.minlinpos, info{d}{e}.maxlinpos);

%    keyboard
%i= find(sum(([behavdata{d}{e}.xpos, behavdata{d}{e}.ypos]-ones(n,1)*z).^2,2)<.25)
%minmax(behavdata{d}{e}.linpos(i))

