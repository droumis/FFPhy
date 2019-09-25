function getlinmeanvel(num)
% function getlinmeanvel(num)
%
% compute mean velocity [pixel/s] vs. linear distance [pixel]
%
% num= [day epoch] (optional) 
% otherwise run all

global fmaux linmeanvel

if nargin < 1
    [d,e,t,c]= startCellList;
    while ~isempty(d)
	auxrun(d,e);
	[d,e]= getNextEpoch;
    end
else
    auxrun(num(1),num(2));
end

save([fmaux.data2dir 'linmeanvel'], 'linmeanvel');
clear global linmeanvel

function auxrun(d,e)

global fmaux linmeanvel vel lindistpos linbehav

loadVar(fmaux.datadir, 'linbehav', d, fmaux.prefix, 1);
loadVar(fmaux.datadir, 'lindistpos', d, fmaux.prefix, 1);
loadVar(fmaux.datadir, 'vel', d, fmaux.prefix, 1);

% match binning in linbehav [@@hack]
t=1;
found= 0;
while isempty(linbehav{d}{e}{t})
    t= t+1;
end
c=1;
while isempty(linbehav{d}{e}{t}{c});
    c= c+1;
end
if isempty(linbehav{d}{e}{t}{c}.data{1,2})
    error('data set without traj 1->2');
end
centers= linbehav{d}{e}{t}{c}.data{1,2}(:,1);
nbins= size(centers,1);
if nbins < 2; error('list of lindistpos must have 2 or more entries'); end
delta= centers(2)-centers(1);
edges= centers- (delta/2);
edges(nbins+1)= centers(nbins)+delta;

ldp= lindistpos{d}{e}.data;
times= vel{d}{e}.data(:,1);
v= vel{d}{e}.data(:,2); 
%max(v)
vi= findIndex(times,ldp(:,1));
lindist= interp1(ldp(:,1), ldp(:,2), times);
type= lindistpos{d}{e}.estinfothetavel(vi);
traj=[1, 2; 2 1; 1 3; 3 1];

mtmp=[];
mtmp(:,1)= centers;
velos=cell(4,1);
for i=1:4  % over trajectories
    ti= find(type== i-1);
    vtmp= cell(nbins,1);
    bi= findIndex(lindist(ti), edges); % find bins
    for k=1:size(ti,1)    % run over velocity array
        vtmp{bi(k)}(end+1)= v(ti(k));
    end
    for k=1:nbins
        if ~isempty(vtmp{k})
            mtmp(k,2)= mean(vtmp{k});
        else
            mtmp(k,2)= nan;
        end
    end
    velos{i}= mtmp;
end    


linmeanvel{d}{e}= [];
linmeanvel{d}{e}.descript= 'traj (X) mean velocity vs. linearized position';
linmeanvel{d}{e}.data= velos;
