function sp= allSpikeData(selectid, loadtheta, nsel)
%function sp= allSpikeData(selectid, loadtheta, nsel)
%
% Return structure that contains all information about spikes.

if nargin<2; loadtheta= 0; end

setRoot
load([root '/data/analist-' selectid])
nC= length(analist.rat);
oldrat= ''; oldd= -1;
sp= [];
for iC= 1:nC
    if nargin>=3 & ~ismember(iC, nsel); continue; end
    % set up data
    rat= analist.rat{iC};
    num= analist.cellnum(iC,:); d=num(1); e=num(2); tet=num(3); c=num(4);
    if ~strcmp(oldrat, rat) | oldd ~= d
        oldrat= rat; oldd= d;
        load(sprintf('/home/chengs/theta/%s/data2/spikedata%.2d.mat', rat, d));
        if loadtheta
            load(sprintf('/home/chengs/theta/%s/data2/behavdata%.2d.mat', rat, d));
        end
    end
    sp(iC).time= spikedata{d}{e}{tet}{c}.time;
    sp(iC).linpos= spikedata{d}{e}{tet}{c}.linpos;
    if loadtheta
        sp(iC).theta= behavdata{d}{e}.phase(spikedata{d}{e}{tet}{c}.index);
    end
end
spikeData= sp;
save(['spikeData-' selectid], 'spikeData')
