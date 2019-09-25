function [ coorind ] = lenscorrtransf( ind )
% Uses indices from correction algorithm to return coordinates
% Compares indices against template matrix
% template matrix located in /home/jai/Src/ImageCorrection/

% reads transformation table
fid=fopen('/home/jai/Src/ImageCorrection/transformation_codes/output_coordinates_after20120514.dat','r');

% text file has 4 columns
% 1st and 2nd: the original pixel coordinates, original top left of image
% 3rd and 4th: transformed pixel coordinates


% for before20120514.dat
% t=textscan(fid,'%f %f', 'delimiter',' ');
% temp=[t{1,1} t{1,2}];
% coorind=temp(ind,1:2); 


% % for after20120514.dat
t=textscan(fid,'%f %f %f %f', 'delimiter',' ');


temp=[t{1,1} t{1,2} t{1,3} t{1,4}]; % temp contains the transformed coordinates for each pixel defined by the order they are loaded

temp=temp(1:76800,:); % the output file has transformation coordinates for all 3 color channels, 3*320*240 values, only need every 3rd

% flip along y axis to get the correct order of pixels
midpointy=120;
temp(:,2)=-(temp(:,2)-midpointy)+midpointy;

midpointtransformedy=120;
temp(:,4)=-(temp(:,4)-midpointtransformedy)+midpointtransformedy;
% 
% % since the new correction algorithm generates corrected pixels not in a linear order,
% % we use the pixel indices in the 1st and 2nd columns to sort the results
% % into the correct order expected by the Matlab code
% 
temp=sortrows(temp,[2 1]);
coorind=temp(ind,3:4); % ind has the pixel order information of the original points


% old version
% for i=1:size(ind,1);
%     currcoord=temp(ind(i,1),:);
%     coorind=[coorind;currcoord];
%     i=i+1;
% end

%end

%max(temp(:,3))