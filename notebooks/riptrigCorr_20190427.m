

andef = animaldef('JZ1');
animaldir = andef{2};
animalprefix = andef{3};
datatype = 'linpos';
linpos = loaddatastruct(animaldir, animalprefix, datatype);
%%
bs = load([animaldir, 'JZ1BehaveState.mat']);
%%
bs.BehaveState
%%

day = 7;
ep = 2;

ax1 = subaxis(2,1,1);
bs_statechangeseq_time = bs.BehaveState.statechanges{day}{ep}.statechangeseq(:,3);
bs_statechangeseq_correct_incorrect = bs.BehaveState.statechanges{7}{2}.statechangeseq(:,[7,8]);
plot(bs_statechangeseq_time, bs_statechangeseq_correct_incorrect, '.');

hold on
linpos_time = linpos{7}{2}.statematrix.time;
linpos_segmentIndex = linpos{7}{2}.statematrix.segmentIndex;
plot(linpos_time, linpos_segmentIndex)
ylabel('segment index / correct')
xlabel('time')
hold off

ax2 = subaxis(2,1,2);
bs_statespace_time = sort(bs.BehaveState.statespace.allepsMat(:,3));
bs_statespace_allbound = bs.BehaveState.statespace.allbound;
ssmode = bs_statespace_allbound(:,1);
sslower = bs_statespace_allbound(:,2);
sshigher = bs_statespace_allbound(:,3);
% day = bs_statespace_allbound(:,5);
% epoch = bs_statespace_allbound(:,6);
errfillAll = fill([bs_statespace_time; flipud(bs_statespace_time)], ...
    [sslower; flipud(sshigher)],[0 0 1],'linestyle','none');
set(errfillAll, 'FaceAlpha', .1)
hold on
plot(bs_statespace_time, ssmode, '.')
ylabel('statespace prob correct')
hold off
        
linkaxes([ax1,ax2],'x')

