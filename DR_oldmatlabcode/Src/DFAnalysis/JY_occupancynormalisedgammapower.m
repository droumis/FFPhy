function [out] = JY_occupancynormalisedgammapower(index, excludetimes,data,linpos, meaneegspectrograms,gamma,tetinfo, varargin)

% plots the occupancy normalised z-scored envelope amplitude for the
% frequency of interest on the track

% index - [day epoch tetrode cell]

appendindex = 0;

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end
warning('OFF','MATLAB:divideByZero');

% get the data
% position
pos=data{index(1)}{index(2)}.Pos.correcteddata(:,2:3);
postimeindex=data{index(1)}{index(2)}.Pos.correcteddata(:,1);
velocity=data{index(1)}{index(2)}.Pos.correcteddata(:,5);

% eeg data
eegdata=meaneegspectrograms{index(1)}{index(2)}{index(4)}.data;
eegsamprate=meaneegspectrograms{index(1)}{index(2)}{index(4)}.samprate;
eegstart=meaneegspectrograms{index(1)}{index(2)}{index(4)}.starttime;


% filter eeg for noise

% define threshold for noise
% 5 for ACC
% 4 for HP

thresholdstd=5;
% define baseline for cuttoff time
% 1 for HP
% 2 for ACC
baseline=1;
event=JY_findhighamplitudeeeg(eegdata,eegsamprate,eegstart,thresholdstd,baseline,0.001);

eegtimeindex=(eegstart:1/eegsamprate:(size(eegdata,1)-1)*1/eegsamprate+eegstart)';
eegdataincludeindex=find(isExcluded(eegtimeindex,event.excludetimes)==0);


% plot where noise happens

badpos=find(isExcluded(postimeindex,event.excludetimes)==1);
plot(pos(:,1),pos(:,2),'Color',[0.7 0.7 0.7]);
hold on;
plot(pos(badpos,1),pos(badpos,2),'.r');





% band envelope amplitude
band=double(gamma{index(1)}{index(2)}{index(4)}.data(:,1));
bandsamprate=gamma{index(1)}{index(2)}{index(4)}.samprate;
bandstart=gamma{index(1)}{index(2)}{index(4)}.starttime;
bandtimeindex=(bandstart:1/bandsamprate:(size(band,1)-1)*1/bandsamprate+bandstart)';
bandenvelope=double(gamma{index(1)}{index(2)}{index(4)}.data(:,3));
bandenvelopemean=mean(bandenvelope);
bandenvelopestd=std(bandenvelope);
bandenvelopez=(bandenvelope-bandenvelopemean)./bandenvelopestd;

%
% figure;
% plot(bandenvelope,'r')
% hold on;
% plot(band,'b');
%
%
figure;
subplot(3,1,1)
plot(velocity);
subplot(3,1,2)
plot(eegdata);

subplot(3,1,3)
plot(bandenvelopez);

% generate time index for envelope data
samprate=double(gamma{index(1)}{index(2)}{index(4)}.samprate);
starttime = gamma{index(1)}{index(2)}{index(4)}.starttime;
timeindex=(starttime:1/samprate:size(bandenvelope,1)*1/samprate+starttime)';

% look up time indices in pos
posindex=find(isExcluded(postimeindex,event.excludetimes)==0);
%posindex=1:length(pos);
includebandenvelop=bandenvelope(isExcluded(bandtimeindex,event.excludetimes)==0);

figure;

plot(bandenvelope);
hold on;
plot(bandenvelope(includebandenvelop),'r');

% get all the value into a 2d histogram of track space
bins=75;

[data1,min1,max1] = rescale(pos(posindex,1),1,bins); % SJ 05/31/10
[data2,min2,max2] = rescale(pos(posindex,2),1,bins); % SJ 05/31/10
xind=round(data1);
yind=round(data2);


% plot wells
cmperpix=data{index(1)}{index(2)}.Pos.cmperpixel;
wellpos=(linpos{index(1)}{index(2)}.wellSegmentInfo.wellCoord);

[wellx,min1,max1] = rescale([pos(:,1); wellpos(:,1)],1,bins); % SJ 05/31/10
[welly,min2,max2] = rescale([pos(:,2); wellpos(:,2)],1,bins); % SJ 05/31/10

wellx=round(wellx(end-1:end));
welly=round(welly(end-1:end));



linindex=sub2ind([100 100],xind',yind');

% mean theta power for each location on track
meanthetapower=accumarray(linindex',bandenvelopez(includebandenvelop),[],@(x)max(x));
% spread of means for each locaion
meanthetastd=accumarray(linindex',bandenvelopez(includebandenvelop),[],@(x)std(x));


arrayindex=unique(linindex)';

[rex rey]=ind2sub([100 100],arrayindex);

meanthetapower2d=accumarray([rex rey],meanthetapower(arrayindex));
meanthetastd=accumarray([rex rey],meanthetastd(arrayindex));

% plot spread of means to see where biggest variation in z score is
% plotstd=1;
% g = gaussian2(plotstd,(6*plotstd));
% smoothed2d = filter2(g,(meanthetastd)); % is this the right filter?
%
% figure;
% imagesc(smoothed2d)
% axis xy;
% axis square;
% colorbar;


%meanthetapower2d(meanthetapower2d==0)=0;

% plot the mean z score on track
figure; set(gcf,'position',[0 0 600 600]); set(gcf,'PaperPositionMode','auto');
plotstd=1;
g = gaussian2(plotstd,(6*plotstd));
smoothed2d = filter2(g,(meanthetapower2d)); % is this the right filter?


cmap=colormap(jet(256));

cmap(1,:)=[0 0 0];


imagesc(meanthetapower2d,[0 2])

hold on;
for wind=1:size(wellx,1);
    plot(wellx(wind),welly(wind),'og','MarkerSize',15);
    
end


%imagesc(meanthetapower2d)
axis xy;
axis square;
colorbar;




close







out.epochmeanzscore = meanthetapower2d;

if appendindex
    out.index = index;
end
warning('ON','MATLAB:divideByZero');
end

