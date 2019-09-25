function debugStat(rat, din, ein, tin, cin, selectid)

d=str2num(din); e= str2num(ein); t= str2num(tin); c= str2num(cin);

%num=[1,4,25,7];
%rat= 'ter';
%selectid= 'placefields4-satellite';
%d=num(1); e= num(2); t= num(3); c= num(4);

num= [d e t c];

load(['/home/chengs/theta/' rat '/data2/select-' selectid]);
load(['/home/chengs/theta/' rat '/data2/behavdata' sprintf('%.2d', d)]);
load(['/home/chengs/theta/' rat '/data2/spikedata' sprintf('%.2d', d)]);

data= behavdata{d}{e};
data.spiketimes= spikedata{d}{e}{t}{c}.time;
data.spikeindex= spikedata{d}{e}{t}{c}.index;

load(['adaptest' sprintf('%.2d', d)]);
model= adaptest{d}{e}{t}{c}.model;

cellId=1;
ncells= size(select.cellnum,1);
while(cellId <= ncells & sum(select.cellnum(cellId,:)==num)~= 4)
        cellId=cellId+1;
end
if(cellId > ncells)
    %    num
        error('could not find requested cell');
    cellId= [];
end

ana.a= select.a{cellId};
ana.x= select.x{cellId};




%%%%%%%%%%%%%%%%%%%

nTimeSteps= 20;

sopts{1}.name= 'MutualInfo';
sopts{1}.nbinx= 20;
sopts{1}.nbiny= 20;

sopts{2}.name= 'LinCorr';
sopts{2}.nbinx= 100;
sopts{2}.nbiny= 100;

sopts{3}.name= 'RescaledSpikes';

sopts{4}.name=  'Integral2d';
sopts{4}.dx= .5;
sopts{4}.dy= .05;

sopts{5}.name= 'Moments2d';
sopts{5}.nbinx= 100;
sopts{5}.nbiny= 100;
sopts{5}.circy= 1;

small= 1e-2;
mint= data.time(1)+small;
maxt= data.time(end)-small;
for i=1:length(ana.a)
    ana.a{i}.phase= [0;2*pi];
    if isempty(ana.x) | (length(ana.x)<i) | ~isfield(ana.x{i},'time')
        ana.x{i}.time=linspace(mint, maxt, nTimeSteps);
    end
end

%%%%%%%%%%%%%%%%%%%
s= calcStatistics(data, model, sopts, ana);

s.Moments2d{1}
%keyboard
