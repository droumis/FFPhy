function showModel(num, opt, bursts)

if nargin < 3; bursts= 0; end
if nargin < 2; opt= []; opt.playrate= 8000; end

if nargin < 1
    [d,e,t,c,ncells]= startCellList;
    while ~isempty(d)
        auxrun(d,e,t,c, bursts);
        [d,e,t,c]= getNextCell;
    end
else
    d=num(1); e=num(2); t=num(3); c=num(4);
    auxrun(d,e,t,c, opt, bursts);
end

function auxrun(d,e,t,c, opt, bursts)


if ~isfield(opt, 'cmap')
    global adaptest cmap
    load /home/chengs/AdaptFilter/vis/colormap;
    %load /home/chengs/AdaptFilter/vis/colormap2;
    opt.cmap= cmap;
end
%opt.cmap= jet(1024);
data=loadData([d e t c], bursts);
loadVar('.','adaptest',d);
visModel(data,adaptest{d}{e}{t}{c}.model,opt);