function setrootdir( newrootdir )

load('/home/ckemere/research/matlab/Lorenlab/DataAnalysis/DataFilter/rootdir.mat');

if nargin > 0
  rootdir = newrootdir;
else
  rootdir = defaultrootdir;
end

save('/home/ckemere/research/matlab/Lorenlab/DataAnalysis/DataFilter/rootdir.mat','rootdir','defaultrootdir');
