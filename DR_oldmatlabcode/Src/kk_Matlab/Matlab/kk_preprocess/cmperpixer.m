function cmperpix = cmperpixer(mpegfile,minute,interval,noframes,objectincm)

% Calculates cmperpix based off of pixel length of a measurable object.

% minute : minute after open all files to start grabbing frames
% interval (in seconds) to skip between grabbed frames
% objectincm : size of object being measured in centimeters

pixellength = [];
cmperpix = [];

for i=1:noframes
    M=mpgread(mpegfile,(minute*60*30+(i-1)*interval*30),'truecolor');
    a=frame2im(M);
    image(a)

    coords=ginput(2)
    pixellength=[pixellength ; dist(coords(1,:),coords(2,:))]
    cmperpix = [cmperpix ; objectincm/pixellength(end)]

end

end

