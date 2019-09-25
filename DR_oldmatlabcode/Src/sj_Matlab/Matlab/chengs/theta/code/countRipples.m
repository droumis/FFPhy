function countRipples
%  
%  Mark timesteps during which ripples occured.


%n= auxGetN('CA1PE', 'ripple');
%nEC= auxGetN('CA1PE', 'ECripple')
nnover= auxGetN('CA1PE', 'EC_nonoverlap_ripple')

keyboard



function n= auxGetN(selectid, ripStr)

setRoot;
epfile= sprintf('%s/data/epochs-%s', root, selectid);
load(epfile);
ep= epochs;
nep= length(ep.rat);
n= 0;

for ie=1:nep
    rat= ep.rat{ie}; d= ep.num(ie,1); e= ep.num(ie,2);
    rip= auxGetRipples(rat, d, e, ripStr);
    [lo hi]= findcontiguous(find(rip));
    n= n+ length(lo);
end
