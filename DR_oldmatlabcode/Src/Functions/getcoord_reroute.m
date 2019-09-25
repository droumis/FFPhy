function coords = getcoord_4arm_maze(directoryname,pos,index)
%this program is called by CREATETASKSTRUCT to produce the trajectory
%coodinates for a reroute maze.
%Click locations in the following order:
%
%   1-------2
%   |\     /|
%   | \   / |
%   |   5   |
%   |  / \  |
%   | /   \ |
%   4-------3 

fid = figure;
plot(pos(:,2),pos(:,3));
[x,y] = ginput(5);
lincoord{1} = [x([1 2]) y([1 2])];
lincoord{2} = [x([2 3]) y([2 3])];
lincoord{3} = [x([3 4]) y([3 4])];
lincoord{4} = [x([4 1]) y([4 1])];    
lincoord{5} = [x([1 5]) y([1 5])]; 
lincoord{6} = [x([5 3]) y([5 3])]; 
lincoord{7} = [x([2 5]) y([2 5])]; 
lincoord{8} = [x([5 4]) y([5 4])]; 


numtimes = size(pos,1);
for i = 1:length(lincoord)
    coords{i} = repmat(lincoord{i},[1 1 numtimes]);
end
close(fid);