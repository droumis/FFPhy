%%
% function cmperpix = cmperpixer(mpegfile,minute,interval,noframes,objectincm)

% Calculates cmperpix based off of pixel length of a measurable object.

% minute : minute after open all files to start grabbing frames
% interval (in seconds) to skip between grabbed frames
% objectincm : size of object being measured in centimeters

cd '/data19/droumis/bob/bob13/';
% mpegfile = 'bob13';
% cmperpixel = cmperpix('bob13.mpeg',2,10,3,33);
clear cmperpix;
clear cmperpixel;

pixellength = [];
minute = 16;
interval = 30;
noframes = 3;
cmperpix = [];
objectincm = 76.8;

for i=1:noframes
%     M=mpgread(mpegfile,(minute*60*30+(i-1)*interval*30),'truecolor');
    
    M=mpgread('bob13.mpeg',(minute*60*30+(i-1)*interval*30),'truecolor');
    
    a=frame2im(M);
    a = rgb2gray(a);
    a = imadjust(a, [.1 .2], []);
    image(a);
%  imcontrast(a);

    coords=ginput(2)
    pixellength=[pixellength ; dist(coords(1,:),coords(2,:))]
    cmperpix = [cmperpix ; objectincm/pixellength(end)]
    cmperpixel(i,1) = cmperpix
    cmperpix = [];

end


