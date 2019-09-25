function runAdaptFilter(num, bursts, iter)
%function runAdaptFilter(num, bursts, iter)
%
% run adaptive filter algorithm for particular cell
%   bursts: only analyze bursts, 0 for all spikes
%   num=[day epoch tetrode cellnumber]
% if num is not given, all cells will be run via startCellList and getNextCell
%
% getLocalFilterModel is called to set the model class and parameters


if nargin < 2 | isempty(bursts)
    opt.bursts= 0;
else
    opt.bursts= bursts;
end
if nargin < 3 | isempty(iter)
    opt.iter= 0;
else
    opt.iter= iter;
end

global adaptest
if nargin < 1 | isempty(num)
    oldday= -1;
    n=1;
    changed_day= 0;
    [d,e,t,c,ncells]= startCellList;
    while ~isempty(d)
        monitorProgress(n, ncells)
        if (oldday~= d)
            auxload(d);
            fname= sprintf('adaptest%.2d',d);
            oldday= d;
            changed_day= 0;
        end
        ran= auxrun(d,e,t,c, opt);
        if ran
            changed_day= 1;
        end
        [d,e,t,c]= getNextCell;
        if (isempty(d) | oldday~= d) & changed_day
            save(fname, 'adaptest');
            changed_day= 0;
            disp('saved');
        end
        n= n + 1;
    end
else
    auxload(num(1));
    ran= auxrun(num(1),num(2),num(3),num(4), opt);
    if ran
        fname= sprintf('adaptest%.2d', num(1));
        save(fname, 'adaptest');
    end
end


function auxload(d)
global adaptest
fname= sprintf('adaptest%.2d.mat',d);
if(exist(fname,'file'))
    loadVar(pwd,'adaptest',d);
else
    fprintf(1,'%s/%s does not exist yet...\n', pwd, fname);
    adaptest= {};
end

function run= auxrun(d,e,t,c, opt)

%fprintf(1,'cell [%d %d %d %d]\n',d,e,t,c);
global adaptest

done= 0;
run= 0;

if length(adaptest) >= d & ~isempty(adaptest{d}) & ...
        length(adaptest{d}) >= e & ~isempty(adaptest{d}{e}) & ...
        length(adaptest{d}{e}) >= t & ~isempty(adaptest{d}{e}{t}) & ...
        length(adaptest{d}{e}{t}) >= c & ~isempty(adaptest{d}{e}{t}{c}) & ...
        isfield(adaptest{d}{e}{t}{c}, 'model')
    struct_exists= 1;
else
    struct_exists= 0;
end


if ~done
    data= loadData([d e t c], opt.bursts);
    data.traj(data.ripple==1)= -1; %%@@
%    keyboard
    if opt.iter==1 | ~struct_exists
        [model, fopts]= getLocalFilterModel(data);
    else
        model= adaptest{d}{e}{t}{c}.model;
        fopts= adaptest{d}{e}{t}{c}.filter;
    end
    fopts.niter= opt.iter;
    [adaptest{d}{e}{t}{c}.model adaptest{d}{e}{t}{c}.filter]= ...
        adaptFilter(data, model, fopts);
    run= 1;
end