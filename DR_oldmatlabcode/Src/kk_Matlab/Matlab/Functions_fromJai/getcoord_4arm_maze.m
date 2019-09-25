function coords = getcoord_4arm_maze(directoryname,pos,index)
%this program is called by CREATETASKSTRUCT to produce the trajectory
%coodinates for a four arm maze.
%Click locations in the following order:
%
%   1                 2
%    \
%
%
%
%
%   1    2    3
%   |    |    |
%   |    |    |
%   |    |    |
%   |    |    |
%   4----5----6
%

fid = figure;
plot(pos(:,2),pos(:,3));
[x,y] = ginput(6);
lincoord{1} = [x([2 5 4 1]) y([2 5 4 1])];
lincoord{2} = [x([2 5 6 3]) y([2 5 6 3])];
    
numtimes = size(pos,1);
for i = 1:length(lincoord)
    coords{i} = repmat(lincoord{i},[1 1 numtimes]);
end
close(fid);