%% aligns oscillatory events, averages, and plots

function out = averageEVENTS(eegstruct,filterstruct,eventstruct,day,epoch,tetrodes,reftetrode,windowsize)

%% eegstruct can be the raw eeg, OR the filtered eeg (consolidated struct)
%% filterstruct data is used to obtain the time index of the peak of the oscillatory event
%% eventstruct contains times of ripples, gammal, hightheta, etc.
%% tetrodes is a vector of tetrodes
%% windowsize is in seconds


Fs=eegstruct{day}{epoch}{1}.samprate;
no_ref_events=length(eventstruct{day}{epoch}{reftetrode}.midind);
halfwin=floor(windowsize*Fs/2);            %% window in samples
searchwin=0.025;                            %% search for peaks within 25 ms
searchwin=floor(searchwin*Fs);

%% plot events
if 0==1
for e=1:1
    startindex=eventstruct{day}{epoch}{reftetrode}.midind(e);                     %% begin search at midpoint of event
    index = findnearestpeak(filterstruct{day}{epoch}{reftetrode}.data,startindex,searchwin);        %% index at cycle peak
    figure
    hold on
    for t=tetrodes
        if t==10            %% DG
%             plot((1:(1+2*halfwin))/Fs, ... 
%                 eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin)), ...
%                 'm','LineWidth',2);
        elseif t==12        %% CA3
            plot((1:(1+2*halfwin))/Fs, ... 
                eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin)), ...
                'r','LineWidth',2);
        elseif t==9        %% CA2
            plot((1:(1+2*halfwin))/Fs, ... 
                eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin)), ...
                'g','LineWidth',2);
        elseif t==8 || t==13      %% near CA2
            plot((1:(1+2*halfwin))/Fs, ... 
                eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin)), ...
                'Color',[0 0.5 0],'LineWidth',2);            
        elseif t==11 || t==14                %% CA1
            plot((1:(1+2*halfwin))/Fs, ...  
                eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin)), ...
                'k','LineWidth',2);
        else
        end
    end
    ylim([-200 200]);
end
end

%% average events

average_eeg=cell(length(tetrodes),1);

for t=tetrodes
    dummy=[];
    for e=1:no_ref_events
        startindex=eventstruct{day}{epoch}{reftetrode}.midind(e);                                       %% begin search at midpoint of event
        index = findnearestpeak(filterstruct{day}{epoch}{reftetrode}.data,startindex,searchwin);        %% index at cycle peak
        dummy = [dummy ; eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin))];
    end
    average_eeg{t}=mean(dummy,1);
end

figure
hold on
title('event-triggered total averages')
no_ref_events
  for t=tetrodes
        if t==10            %% DG
             plot((-halfwin:halfwin)/Fs, ... 
                average_eeg{t}, ...
                 'm','LineWidth',2);
        elseif t==12        %% CA3
             plot((-halfwin:halfwin)/Fs, ... 
                average_eeg{t}, ...
                 'r','LineWidth',2);
        elseif t==9        %% CA2
             plot((-halfwin:halfwin)/Fs, ... 
                average_eeg{t}, ...
                 'g','LineWidth',2);
        elseif t==8 || t==13      %% near CA2
             plot((-halfwin:halfwin)/Fs, ... 
                average_eeg{t}, ...
                 'Color',[0 0.5 0],'LineWidth',2);          
        elseif t==11 || t==14                %% CA1
             plot((-halfwin:halfwin)/Fs, ... 
                average_eeg{t}, ...
                 'k','LineWidth',2);
        else
        end
    end




end



