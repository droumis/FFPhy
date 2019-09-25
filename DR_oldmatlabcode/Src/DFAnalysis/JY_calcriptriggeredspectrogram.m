function out = JY_calcriptriggeredspectrogram(index, excludeperiods, data,eeg,meaneegspectrograms, ripples,ripple,tetinfo, varargin)
% out = calcripspectrum(index, excludeperiods, eeg,ripples,cellinfo, options)
% Based on Maggie's calcriptriggeredspectrogram 
% Computes the spectrogram around the middle of each decoded event.
% eeg - eeg data
% meaneegspectrum - dereferenced eeg data
% ripples - extracted ripple data
% tetinfo - tetrode information
%  Options:
%       appendindex-- Determines whether index is included in output vector
%           Default: 1
%       fpass-- Determines the frequency range for computing spectrum.
%           Default: [2 350]
%       average_trials-- Determines if events are averaged or not.
%           Default: 0
%       spectrum_window-- Determines the sliding window used to compute
%           the event triggered spectrogram. Default: [0.1 0.01]
%       event_window--Determines the size of the window around each
%           triggering event. Default: [0.2 0.2]
%       cellfilter--Determines which tetrodes to use for ripple extraction.
%           Default is 'isequal($area, ''CA1'') & numcells>1)'
%  out is a structure with the following fields
%       S-- This is a MxNxT matrix of the spectrogram for each tetrode
%           M is time relative to triggering event, N is frequency, T is event
%       F-- Frequency vector
%       T-- time relative to triggering event
%       fit-- This is the fit based on the spectrum computed for the entire
%           epoch to normalize S. To reconstruct S without normalization,
%           add log10(frequency)*fit(2)+fit(1) to S.
%       index-- Only if appendindex is set to 1 (default)


%parse the options
params = {};
params.Fs = 1500;
params.fpass = [0 400];
params.trialave = 0;
cwin = [250 10]/1000;
win = [500 500]/1000;
appendindex = 1;
cellfilter = [];
params.tapers = [3 5];

for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'aveargetrials'
                params.trialave = varargin{option+1};
            case 'spectrum_window'
                cwin = varargin{option+1};
            case 'event_window'
                win = varargin{option+1};
            case 'cellfilter'
                case 'tetrodereference'
            tetref=varargin{option+1};
                
                
                cellfilter = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% if isempty(cellfilter)
%     riptimes = getripples(index, ripples, cellinfo, 'cellfilter', '(isequal($area, ''CA1''))','excludeperiods', excludeperiods,'minstd',3);
%         
%     %riptimes = getriptimes(index, ripples, cellinfo, 'cellfilter', '(isequal($area, ''CA1''))','excludeperiods', excludeperiods,'minstd',3);
% 
% 
% else
%     riptimes = getripples(index, ripples, cellinfo, 'cellfilter', cellfilter,'excludeperiods', excludeperiods,'minstd',3);
% end


% get position from data

pos=data{index(1)}{index(2)}.Pos.correcteddata(:,2:3);

% use riptimes to get ripples
%riptimes=JY_getripples(index, ripples,'excludeperiods', excludeperiods,...
%    'minstd',5,'tetrodereference',tetref);

% get ripples manually from one tetrode

getripindex=find(ripples{index(1)}{index(2)}{tetref}.maxthresh>=5);
riptimes=ripples{index(1)}{index(2)}{tetref}.midtime(getripindex)*10000;
ripplesize=ripples{index(1)}{index(2)}{tetref}.maxthresh(getripindex);
rippleposind=ripples{index(1)}{index(2)}{tetref}.posind(getripindex);


if ~isempty(riptimes)
    
% Define EEG
e = meaneegspectrograms{index(1)}{index(2)}{index(4)}.data';
starttime = meaneegspectrograms{index(1)}{index(2)}{index(4)}.starttime;
endtime = (length(e)-1) * (1 / params.Fs);
clear eeg cellinfo


% remove events during recording noise
eegsamprate=meaneegspectrograms{index(1)}{index(2)}{index(4)}.samprate;
% filter eeg for noise

% define threshold for noise
% 5 for ACC
% 4 for HP

thresholdstd=5;
% define baseline for cuttoff time
% 1 for HP
% 2 for ACC
baseline=2;
event=JY_findhighamplitudeeeg(e,eegsamprate,starttime,thresholdstd,baseline,0.001);




% Define triggering events as the start of each ripple
triggers = riptimes(find(~isExcluded(riptimes(:,1)/10000,event.excludetimes)),1)/10000-starttime;

%Remove triggering events that are too close to the beginning or end
while triggers(1)<win(1)
    triggers(1) = [];
end
while triggers(end)> endtime-win(2)
    triggers(end) = [];
end

% Calculate the event triggered spectrogram
[S,t,f] = mtspecgramtrigc(e,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);

% Compute a z-scored spectrogram using the mean and std for the entire session
P = mtspecgramc(e,[cwin(1) cwin(1)],params);
%clear %e %triggers %riptimes
meanP = mean(P);
stdP = std(P);
clear P
for i = 1:size(S,1)
    for j = 1:size(S,3)
    	S(i,:,j) = (S(i,:,j) - meanP)./stdP;
    end
end


% get raw eeg trace


% % 
% for ii=1:size(S,3);
% 
% figure;
% 
% set(gcf,'position',[100 100 900 600]);
% 
% subplot(2,3,1);
% 
% plot(pos(:,1),pos(:,2),'color',[0.7 0.7 0.7]);
% hold on;
% 
% plot(pos(rippleposind(ii),1),pos(rippleposind(ii),2),'+r','MarkerSize',20);
% 
% axis square;
% 
% 
% subplot(2,3,2:3);
% 
% eegstartind=ceil((triggers(ii)-win(1))*1500);
% eegendind=ceil((triggers(ii)+win(1))*1500);
% 
% plot(e(eegstartind:eegendind));
% xlim([0 abs(eegstartind-eegendind)]);
% 
% hold on;
% 
% line([0.5*abs(eegstartind-eegendind) 0.5*abs(eegstartind-eegendind)],...
%     [-500 500],'color','r');
% 
% colorbar;
% 
% subplot(2,3,5:6);
% 
%                     %imagesc(t-win(1),f,S(:,:,ii)');
%                      imagesc(t-win(1),f,S(:,:,ii)');
%                     
%                     %imagesc(tmt-win(1),fmt,epocheegmean',[0, 20]);
%                     set(gca,'YDir','normal');
%                     string = sprintf('%d',t);
%                     %title(string);
%                     colormap('jet')
%                     colorbar
%                     
%                h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
%                         1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%                     
%                     text(0.5, 1,sprintf('EEG and spectrum aligned to SWR, ripple %s %s std', ...
%                         num2str(ii),num2str(ripplesize(ii),2)),'HorizontalAlignment'...
%                         ,'center','VerticalAlignment', 'top','FontSize',8);       
%             
%                     close
%                     
% end


% get eeg around ripple

epochripple=ripple{index(1)}{index(2)}{tetref};

%rippleripple=cell(1,size(riptimes,1));

halfwindow1=win(1)*10000;
halfwindow2=win(2)*10000;
eegsamprate=1500;

%riptimes=riptimes;
i=1;
while i<length(riptimes)
    % go through each tetrode
    tstart=riptimes(i,1)-halfwindow1;
    tend=riptimes(i,1)+halfwindow2;
    %tend=rippletimes(i,2);+halfwindow;
    
      %ripple data
        ripplestart=epochripple.starttime*10000;
         
        tstartindex=ceil((tstart-ripplestart)/(10000/eegsamprate));
        if tstartindex <1 % check for mismatch between epoch start times defined by times.mat and position reconstruction
            tstartindex=1;
        end
        tendindex=ceil((tend-ripplestart)/(10000/eegsamprate));
        
        if tendindex>size(epochripple.data,1);
            tendindex=size(epochripple.data,1);
           i=length(riptimes)+1;
        end
        
        rippleripple{i}=double(epochripple.data(tstartindex:tendindex,3));
        i=i+1;
    
end

meanripple=mean(cell2mat(rippleripple),2);

out.epochmeanzscore = mean(S,3);
out.times = t - win(1);
out.freq = f;
out.meanrippleeeg=meanripple;
out.epochrawzscore=S;
    

else
  meanripple=ones(10,10);

out.epochmeanzscore = ones(10,1);
out.times = 10;
out.freq = 10;
out.meanrippleeeg=ones(10,1);
out.epochrawzscore=ones(10,10);  
end

if appendindex
	out.index = index;
end
end