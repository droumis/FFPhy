function runStats(num)
%function runStats(num)
%
% Run statistics calculation. Use local options via getLocalStatOpts.m

global fmaux stats select
stats= [];

load(fmaux.select);

% load stats file if it already exists
fname= ['stats-' fmaux.selectid '.mat'];
if exist(fname,'file')
    load(fname);
end

if nargin < 1
    [d,e,t,c]= startCellList;
    while ~isempty(d)
        auxrun(d,e,t,c);
        [d,e,t,c]= getNextCell;
    end
else
    fmaux.currentCell= getSelectId(num);
    auxrun(num(1),num(2),num(3),num(4));
end

save(fname, 'stats');


function auxrun(d,e,t,c)

global stats fmaux select adaptest

if isempty(fmaux.currentCell) | ~isfield(select,'a') | ~isfield(select,'x')
    ana= [];
else
    ana.a= select.a{fmaux.currentCell};
    ana.x= select.x{fmaux.currentCell};

    if(isempty(ana.a)) return; end
end
data= loadData([d e t c]);
loadVar('.','adaptest', d);
%model= adaptest{d}{e}{t}{c}.model.xp;
model= adaptest{d}{e}{t}{c}.model;
[sopts, ana]= getLocalStatOpts(data,ana);

s= calcStatistics(data, model, sopts, ana);

if ~isempty(s)
    names= fieldnames(s);
    for i=1:length(names)
        stats{d}{e}{t}{c}.(names{i})= s.(names{i});
    end
end
stats{d}{e}{t}{c}.a= ana.a;
stats{d}{e}{t}{c}.x= ana.x;
