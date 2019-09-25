function coords = getcoord_reroute_arms(directoryname,pos,index)
%this program is called by CREATETASKSTRUCT to produce the trajectory
%coodinates for a reroute maze.
%Click locations in the following order:
% 6           7
%  \         /  
%   1-------2
%   |\     /|
%   | \   / |
%   |   5   |
%   |  / \  |
%   | /   \ |
%   4-------3 
%  /         \ 
% 9           8   

fid = figure;
plot(pos(:,2),pos(:,3));
xlim([0 160]); % defines figure dimensions for specifying the points below, make bigger than the position data to see entire arena
ylim([0 120]); % resolution * pixel dimension
[x,y] = ginput_jai(9);
lincoord{1} = [x([1 2]) y([1 2])];
lincoord{2} = [x([2 3]) y([2 3])];
lincoord{3} = [x([3 4]) y([3 4])];
lincoord{4} = [x([4 1]) y([4 1])];    
lincoord{5} = [x([1 5]) y([1 5])]; 
lincoord{6} = [x([5 3]) y([5 3])]; 
lincoord{7} = [x([2 5]) y([2 5])]; 
lincoord{8} = [x([5 4]) y([5 4])];
lincoord{9} = [x([1 6]) y([1 6])];
lincoord{10} = [x([2 7]) y([2 7])];
lincoord{11} = [x([3 8]) y([3 8])];
lincoord{12} = [x([4 9]) y([4 9])];

numtimes = size(pos,1);
for i = 1:length(lincoord)
    coords{i} = repmat(lincoord{i},[1 1 numtimes]);
end
close(fid);