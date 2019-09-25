function coords = getcoord_utrack(directoryname,pos,index)
%this program is called by CREATETASKSTRUCT to produce the trajectory
%coodinates for a wtrack.
%Click locations in order:
%   

fid = figure;
plot(pos(:,2),pos(:,3));
[x,y] = ginput(6);
lincoord{1} = [x(:) y(:)];

    
numtimes = size(pos,1);
for i = 1:length(lincoord)
    coords{i} = repmat(lincoord{i},[1 1 numtimes]);
end
close(fid);
