directoryname = '/data14/jai/H2_/';
fileprefix = 'H2';
days = 1:10;
tetrodes = -1; % specify all tetrodes
mindur = 0.015; % 15 ms minimum duration
nstd = 2; % 2 std
for d = days
    extractripples_JY(directoryname, fileprefix, d, tetrodes, mindur, nstd);
end