function trajplot(pos, statevector,index)
%trajplot(pos, statevector,index)

%Plots the positions of each trajectory
%
%statevector - the output of GETBEHAVESTATE. This is a vector with the traj
%              number for each position time (1 based). -1 values signify
%              invalid times and are not used.
%index - [day epoch]



pos = pos{index(1)}{index(2)}.data;
%statevector = statevector{index(1)}{index(2)};
timestep = pos(2,1) - pos(1,1);
trajnum = max(statevector);

try        
   trajdata = cell(1,trajnum);
   goodlocationind = (find(statevector ~= -1 ));
   goodlocations = [pos(goodlocationind,[2 3]) statevector(goodlocationind)]; %the x, y locations and trajnum at all valid times
end

fig1 = figure;
for i = 1:length(trajdata)
   
   %get all the linear locations when the animal was on the ith
   %trajectory
   switch mod(i,6)
    case 1
      col = 'b';
    case 2
      col = 'c';
    case 3
      col = 'r';
    case 4
      col = 'm';
    case 5
      col = 'g';
    case 0
      col = 'y';
   end
   dot = '.';
   if ((i > 6) & (i <= 12))
        dot = '+';
   elseif (i>=13)
        dot = '*';
   end
   figure(fig1);
   tmploc = goodlocations(find(goodlocations(:,3) == i),[1 2]);
   plot(tmploc(:,1),tmploc(:,2),[col,dot])
   hold on
   figure
   plot(tmploc(:,1),tmploc(:,2),[col,dot])
   

end




