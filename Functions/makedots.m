function [list]= makedots( xdim, ydim, xspace, yspace )
% makes a x y matrix for a grid of points on a xy plane

list=[];
xlist=[];

i=1;


while i<xdim;
    ylist=[];
    for j=1:yspace:ydim;
        currylist=[];
        currylist(1,1)=i;
        currylist(1,2)=j;
        ylist=[ylist;currylist];
        j=j+yspace;
    end
    list=[list;ylist];
    i=i+xspace;
end

end

