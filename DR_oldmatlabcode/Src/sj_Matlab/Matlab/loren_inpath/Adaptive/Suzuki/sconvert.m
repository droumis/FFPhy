% write out all of the data to binary files
% spikes (used for viewing movies)
%cd c:/Loren_stuff
pushd /data4/loren/Adaptive/Suzuki
fid = fopen('spiketimes.dat', 'w');
%fwrite(fid, spiketimes{1}, 'float');
fwrite(fid, spiketimes, 'float');
fclose(fid);

%Time
fid = fopen('time.dat', 'w');
fwrite(fid, timestep, 'float');
fwrite(fid, fixationtimes, 'float');
fclose(fid);



fid = fopen('trialinfo.dat', 'w');
fwrite(fid, ntrials, 'int');
fwrite(fid, testID, 'int');
fwrite(fid, ntimesteps, 'int');
fwrite(fid, (cobj.response == 0), 'int');
fclose(fid);

% estimated field
%tmp = [CP.x(1) CP.x CP.x(end)];
fid = fopen('thetahate.dat', 'w');
fwrite(fid, length(CP.x), 'float');
fwrite(fid, length(CP.t), 'float');
fwrite(fid, CP.x, 'float');
fwrite(fid, CP.t, 'float');
fwrite(fid, thetainit, 'float');
fwrite(fid, thetahat, 'float');
clear t;
fclose(fid);

popd
