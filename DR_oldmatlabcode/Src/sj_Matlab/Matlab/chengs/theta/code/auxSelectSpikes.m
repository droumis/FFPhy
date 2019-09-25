function valid= auxSelectSpikes(rat, num, spsel, pl, ip, opt);

global behavdata spikedata
setRoot
data2dir= fullfile(root,rat,'data2');
loadVar(data2dir, 'behavdata', num(1));
be= behavdata{num(1)}{num(2)};

ind= findstr(spsel, '_');
ncrit= length(ind)+1;
ind= [0, ind, length(spsel)+1];
crit= {};
for ic=1:ncrit
    crit{ic}= spsel(ind(ic)+1:ind(ic+1)-1);
    if strfind(crit{ic}, 'ripple')
        if strfind(crit{ic}, 'non')
            ripStr{ic}= crit{ic}(4:end);
            crit{ic}= '**nonRIP**';
        else
            ripStr{ic}= crit{ic};
            crit{ic}= '**RIP**';
        end
    end
end

N= length(be.time);

if nargin<6; opt= []; end
if isfield(opt, 'timeid')
    tm= getTimes(opt.timeid, pl.rat{ip}, pl.cellnum(ip,1), ...
        pl.cellnum(ip,2), pl.traj{ip}); 
    if isfinite(opt.t(end))
        tb= minmax(tm(opt.t));
    else
        tb= minmax(tm(opt.t(1):end));
    end
    valid= be.time>tb(1) & be.time<tb(end);
else
    valid= ones(N,1);
end


for ic=1:ncrit
    switch crit{ic}
    case 'all'
    case 'run'
        valid= valid & (be.traj>=0);
    case 'nonrun'
        valid= valid & (be.traj<0);
    case '**RIP**'
        rip= auxGetRipples(rat, num(1), num(2), ripStr{ic});
        valid= valid & rip;
    case '**nonRIP**'
        rip= auxGetRipples(rat, num(1), num(2), ripStr{ic});
        valid= valid & ~rip;
    case {'placefield'}
        x= be.linpos;
        npf= size(pl.pf{ip}, 2);
        within= zeros(N,1);
        for ipf=1:npf
            minx= min(pl.pf{ip}(:,ipf)); maxx= max(pl.pf{ip}(:,ipf));
            within= within | (be.traj==pl.traj{ip}(ipf) & minx<x & x<maxx);
        end
        valid= valid & within;
    case {'min'}
    otherwise
        error(['Unknown spike selection filter "' spsel '".']);
    end

end
