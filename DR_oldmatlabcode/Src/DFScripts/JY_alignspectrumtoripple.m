function [out] = JY_alignspectrumtoripple(index, excludetimes,eeg,meaneegspectrograms,ripples,tetinfo, varargin)

%
% get the times of ripples as specified by filtering
% excludetime -- times of ripples
% does spectral analysis time locked to ripples
% subtracts mean of the epoch


% index - [day epoch tetrode ]



appendindex = 0;
std = 1;
binsize = 1;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'std'
                std = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
                  case 'tetrodereference'
            tetref=varargin{option+1};
                
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end
warning('OFF','MATLAB:divideByZero');

% get times of ripples
excludeperiods=excludetimes;

rippletimes = JY_getripples(index, ripples,'excludeperiods', excludeperiods,'minstd',5,'tetrodereference',tetref);
  
%out.ripplestart=rippletimes;

% sets up parameters

days=index(1);
epoch=index(2);
tetrodes=index(4);


eegstruct=eeg;



%% use mtspectrumtrigc to get spectrogram around trigger times

params = {};
params.Fs = 1500;
params.fpass = [10 300];
params.trialave = 0;
movingwin = [125 10]/1000;
win=[1000 1000]/1000;
params.tapers = [3 5];

data = meaneegspectrograms{days}{epoch}{tetrodes}.data;

%data=meaneegspectrograms{days}{epoch}{tetrodes}.data;

starttime=meaneegspectrograms{days}{epoch}{tetrodes}.starttime;

% filter eeg data to remove high amplitude noise
eegsamprate=meaneegspectrograms{days}{epoch}{tetrodes}.samprate;



% filter eeg for noise

% define threshold for noise
% 5 for ACC
% 4 for HP

thresholdstd=5;
% define baseline for cuttoff time
% 1 for HP
% 2 for ACC
baseline=2;
event=JY_findhighamplitudeeeg(data,eegsamprate,starttime,thresholdstd,baseline,0.001);




endtime = (length(data)-1) * (1 / params.Fs);

% Define triggering events as the start of each ripple
triggers = rippletimes(:,1)-starttime;



% remove triggers within window time from star[S,t,f] = mtspecgramtrigc(data,triggers(startindex:endindex),win,movingwin,params);t and end
startindex=find(triggers>win(1),1,'first');

endindex=find((abs(triggers-endtime))>win(2),1,'last');

%% Calculate the event triggered spectrogram
[Smt,tmt,fmt] = mtspecgramtrigc(data,triggers(startindex:endindex),win,movingwin,params);

%% z score each spectrogram

% precalculated mean epoch spectrograms

epochmean=meaneegspectrograms{days}{epoch}{tetrodes}.epochmean;
epochstd=meaneegspectrograms{days}{epoch}{tetrodes}.epochstd;
for r=1:size(Smt,3)
    zscorespectrograms(:,:,r)=bsxfun(@minus,Smt(:,:,r),epochmean);
    zscorespectrograms(:,:,r)=bsxfun(@rdivide,Smt(:,:,r),epochstd);
end

%epocheegmean=mean(zscorespectrograms,3);
epocheegmean=mean(Smt,3);

% for ii=1:10
%                     figure;
%                     imagesc(tmt-win(1),fmt,zscorespectrograms(:,:,ii)');
%                     %imagesc(tmt-win(1),fmt,epocheegmean',[0, 20]);
%                     set(gca,'YDir','normal');
%                     string = sprintf('%d',tmt);
%                     %title(string);
%                     colormap('jet')
%                     colorbar
% end
% figure;
% imagesc(tmt-win(1),fmt,epocheegmean');
% set(gca,'YDir','normal');
%                     string = sprintf('%d',tmt);
%                     %title(string);
%                     colormap('jet')
%% Generate output

out.epochmeanzscore=epocheegmean;
out.epochmean=epochmean;
out.params=params;
out.times=tmt-win(1);
out.freq=fmt;




warning('ON','MATLAB:divideByZero');