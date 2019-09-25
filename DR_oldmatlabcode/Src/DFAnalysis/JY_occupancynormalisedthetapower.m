function [out] = JY_occupancynormalisedthetapower(index, excludetimes,data,meaneegspectrograms,theta,tetinfo, varargin)

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
baseline=2;
event=JY_findhighamplitudeeeg(eegdata,eegsamprate,eegstart,thresholdstd,baseline,0.001);

% filter eeg data for high amplitude noise

eegdatamean=mean(eegdata);
eegdatastd=std(eegdata);
eegdataz=(eegdata-eegdatamean)./eegdatastd;

% exlcude times when eegdataz is beyond 3 stds

eegtimeindex=(eegstart:1/eegsamprate:(size(eegdata,1)-1)*1/eegsamprate+eegstart)';
eegdataincludeindex=find(isExcluded(eegtimeindex,event.excludetimes)==0);

% band envelope amplitude
band=double(theta{index(1)}{index(2)}{index(4)}.data(:,1));
bandenvelope=double(theta{index(1)}{index(2)}{index(4)}.data(:,3));
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
samprate=double(theta{index(1)}{index(2)}{index(4)}.samprate);
starttime = theta{index(1)}{index(2)}{index(4)}.starttime;
timeindex=(starttime:1/samprate:size(bandenvelope,1)*1/samprate+starttime)';


% look up time indices in pos
posindex=find(isExcluded(postimeindex,event.excludetimes)==0);

includebandenvelop=lookup(postimeindex(posindex),timeindex);




% get all the value into a 2d histogram of track space
bins=100;

[data1,min1,max1] = rescale(pos(posindex,1),1,bins); % SJ 05/31/10
[data2,min2,max2] = rescale(pos(posindex,2),1,bins); % SJ 05/31/10

xind=round(data1);
yind=round(data2);

linindex=sub2ind([100 100],xind',yind');

% mean theta power for each location on track
meanthetapower=accumarray(linindex',bandenvelopez(includebandenvelop),[],@(x)mean(x));
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


imagesc(smoothed2d,[-0.5 1])
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

