function coordinates= get_coord_yuri(directoryname,pos,index)
number_of_times=size(pos,1); k=2;
load([directoryname, 'diode_locations.mat'])

% 'diodes(:,:,k)' specify the diode positions for every 1<k<number_of_moves;   
% characteristics(k,:)=[t(flash_frames(k)) flash_frames(k) k] characterize the positions
%                              ^                 ^
%                              |                 |
%                        time of flashes    frame_numbers

 times_of_frames=characteristics(:,1);         % the times at which actual frames occured
 coordinate_array=zeros(11,2,number_of_times); % At each time, specify the track configuration 
 
diode_shift=zeros(11,2); %scale=[(300/241.55) (377/321.60)]; 
scale=[1 1];
diode_shift(1,:)=scale.*[0 22.5]; diode_shift(8,:)=scale.*[0 -22.5];
diode_shift(2,:)=scale.*[0 22.5]; diode_shift(9,:)=scale.*[0 -22.5];
diode_shift(3,:)=scale.*[0 22.5]; diode_shift(10,:)=scale.*[0 -22.5];
diode_shift(4,:)=scale.*[0 22.5]; diode_shift(11,:)=scale.*[0 -22.5];
diode_shift(5,:)=scale.*[-7.0 0]; 
diode_shift(6,:)=scale.*[-7.0 0]; 
diode_shift(7,:)=scale.*[-7.0 0];
  
  for j=1:number_of_times
      
        while(times_of_frames(k)<pos(j,1))
            k=k+1;
        end
        
        coordinate_array(:,:,j)=diodes(:,:,k-1)+diode_shift(:,:);
  end
   
coordinates{1,1}=coordinate_array;