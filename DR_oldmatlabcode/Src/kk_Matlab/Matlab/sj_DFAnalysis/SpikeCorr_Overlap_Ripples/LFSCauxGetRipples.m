function [rip, frac]= auxGetRipples(rat, d, e, ripStr)
%function rip= auxGetRipples(rat, d, e, ripStr)
%
% Return array indicating whether a ripple occured during timestep.
% Shantanu: Sen Chengs from Lorens directory

if nargin < 4; ripStr= ''; end

global behavdata ripples
setRoot
datadir= fullfile(root,rat,'data');
data2dir= fullfile(root,rat,'data2');

invalid= [];

switch (ripStr)
case {'ripple40ms'}
    loadVar(data2dir,'ripples',0,'',0,'ripples_40ms');
    rip= ripples{d}{e}.data;
case {'ripple3'}
    loadVar(data2dir,'ripples',0,'',0,'ripples_3std');
    rip= ripples{d}{e}.data;
case {'ripple35'}
    loadVar(data2dir,'ripples',0,'',0,'ripples_3_5std');
    rip= ripples{d}{e}.data;
case {'ripple4'}
    loadVar(data2dir,'ripples',0,'',0,'ripples_4std');
    rip= ripples{d}{e}.data;
case {'ripple5'}
    loadVar(data2dir,'ripples',0,'',0,'ripples_5std');
    rip= ripples{d}{e}.data;
case {'ripple6'}
    loadVar(data2dir,'ripples',0,'',0,'ripples_6std');
    rip= ripples{d}{e}.data;
case {'', 'ripple'}
    loadVar(data2dir, 'behavdata', d);
    rip= behavdata{d}{e}.ripple;
case 'xripple'  % 3StdOver, <2.8cm/s
    loadVar(data2dir, 'behavdata', d);
    rip= behavdata{d}{e}.ripple;
    invalid= behavdata{d}{e}.traj>=0;
case 'yripple'  % 3StdOver, <4cm/s
    loadVar(data2dir, 'behavdata', d);
    valid= behavdata{d}{e}.traj<0;
    nlen= length(valid);
    [lo,hi]= findcontiguous(find(valid));
    pad= 250;
    lo= lo-pad; lo(lo<1)= 1;
    hi= hi+pad; hi(hi>nlen)= nlen;
    valid(makecontiguous(lo,hi))=1;
    invalid= ~valid;
    rip= behavdata{d}{e}.ripple;
case 'eripple'
    global lindistpos
    loadVar(datadir, 'lindistpos', d, rat, 1);
    valid=  lindistpos{d}{e}.estinfothetavel>=0;
    nlen= length(valid);
    [lo,hi]= findcontiguous(find(valid));
    pad= 1200;
    lo= lo-pad; lo(lo<1)= 1;
    hi= hi+pad; hi(hi>nlen)= nlen;
    valid(makecontiguous(lo,hi))=1;
    invalid= ~valid;
    loadVar(data2dir, 'behavdata', d);
    rip= behavdata{d}{e}.ripple;
case 'fripple'
    global lindistpos
    loadVar(datadir, 'lindistpos', d, rat, 1);
    valid=  lindistpos{d}{e}.estinfothetavel>=0;
    nlen= length(valid);
    [lo,hi]= findcontiguous(find(valid));
    pad= 250;
    lo= lo-pad; lo(lo<1)= 1;
    hi= hi+pad; hi(hi>nlen)= nlen;
    valid(makecontiguous(lo,hi))=1;
    invalid= ~valid;
    loadVar(data2dir, 'behavdata', d);
    rip= behavdata{d}{e}.ripple;
case 'gripple'
    global lindistpos
    loadVar(datadir, 'lindistpos', d, rat, 1);
    valid=  lindistpos{d}{e}.estinfothetavel>=0;
    invalid= ~valid;
    loadVar(data2dir, 'behavdata', d);
    rip= behavdata{d}{e}.ripple;
case 'iripple'
    global lindistpos
    loadVar(datadir, 'lindistpos', d, rat, 1);
    valid=  lindistpos{d}{e}.estinfothetavel>=0;
    nlen= length(valid);
    [lo,hi]= findcontiguous(find(valid));
    pad= 250;
    lo= lo-pad; lo(lo<1)= 1;
    hi= hi+pad; hi(hi>nlen)= nlen;
    valid(makecontiguous(lo,hi))=1;
    invalid= valid;
    loadVar(data2dir, 'behavdata', d);
    rip= behavdata{d}{e}.ripple;
case 'tripple'
    global lindistpos
    loadVar(datadir, 'lindistpos', d, rat, 1);
    loadVar(data2dir, 'behavdata', d);
    invalid= behavdata{d}{e}.traj>=0 | lindistpos{d}{e}.estinfothetavel>=0;
    rip= behavdata{d}{e}.ripple;
case '6std'
    rip= behavdata{d}{e}.ripple2;
case 't6ripple'
    global lindistpos
    loadVar(datadir, 'lindistpos', d, rat, 1);
    invalid= behavdata{d}{e}.traj>=0 | lindistpos{d}{e}.estinfothetavel>=0;
    rip= behavdata{d}{e}.ripple2;
case 'ECripple'
    global ripples
    loadVar(data2dir,'ripples',0,'',0,'ripples_ECPE_3std_15ms');
    if length(ripples)>=d & length(ripples{d})>=e & ~isempty(ripples{d}{e})
        rip= ripples{d}{e}.rip;
    else
        loadVar(data2dir, 'behavdata', d);
        rip= false(length(behavdata{d}{e}.time),1);
    end
case 'EC_nonoverlap_ripple'
    loadVar(data2dir, 'behavdata', d);
    global ripples
    loadVar(data2dir,'ripples',0,'',0,'ripples_ECPE_3std_15ms');
    if length(ripples)>=d & length(ripples{d})>=e & ~isempty(ripples{d}{e})
        rip= ripples{d}{e}.rip;
    else
        rip= false(length(behavdata{d}{e}.time),1);
    end
    invalid= behavdata{d}{e}.ripple;
case 'nongamma3'
    loadVar(data2dir, 'behavdata', d);
    rip= behavdata{d}{e}.ripple;

    global gamma
    loadVar(data2dir,'gamma',0,'',0,'gamma_3std');
    if length(gamma)>=d & length(gamma{d})>=e
        invalid= gamma{d}{e}.elev;
    else
        invalid= false(length(behavdata{d}{e}.time),1);
    end
case 'gamma3'
    loadVar(data2dir, 'behavdata', d);
    rip= behavdata{d}{e}.ripple;

    global gamma
    loadVar(data2dir,'gamma',0,'',0,'gamma_3std');
    if length(gamma)>=d & length(gamma{d})>=e
        invalid= ~gamma{d}{e}.elev;
    else
        invalid= true(length(behavdata{d}{e}.time),1);
    end
otherwise
    error('unknown ripStr');
end


if isempty(invalid); return; end

[lo,hi]= findcontiguous(find(rip));
nrip= length(lo);
nrej= 0;
for ir=1:nrip
%    @@
    if any(invalid(lo(ir):hi(ir)))
%    if all(invalid(b(ir,1):b(ir,2)))
        rip(lo(ir):hi(ir))= false;
        nrej= nrej+1;
    end
end
if nrip==0 
    frac= nan;
else
    frac= nrej/ nrip;
end

