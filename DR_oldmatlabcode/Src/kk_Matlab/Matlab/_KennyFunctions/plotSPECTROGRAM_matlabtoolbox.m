
%% Plots continuous spectrogram (z-scored to entire day).
   % ingredients
        % 1. eegstruct
        % 2. pos (to overlay velocity on spectrogram)

%% Set these parameters manually.

loaddata = 0;
plot_only = 1;
reward_flag = 1;
days = 4;
epochs = [2 4 6];
tetrodes = [17];

winduration = 15;  % in seconds

if loaddata
    eegstruct = loadeegstruct('/data12/mari/Egy/','egy','eeg',days);
    pos = loaddatastruct('/data12/mari/Egy/','egy','pos',days);
end
if reward_flag
    rewardinfo = loaddatastruct('/data12/mari/Egy/','egy','rewardinfo',days);
end

Fs = 1500;                                      

% matlab sig proc spectrogram f'n params
specgramwindow = 1*Fs;   % in samples
noverlap = 0.9*Fs;  % in samples

tetno=tetrodes(end);
epno=epochs(end);


%% First, calculate the continuous spectrogram for all epochs for the day.
%% At the same time, concatenate the spectrograms from each epoch and
%% calculate the mean and std.

if ~plot_only

spectrograms={};
meandayspectrum={};
stddayspectrum={};
timevec={};

for d=days
    for t=tetrodes                                
        dummy=[];
        for e=epochs
            [~,freqs,timevec{d}{e}{t},S] = spectrogram(eegstruct{days}{e}{t}.data,specgramwindow,noverlap,0:.1:14,Fs);
            spectrograms{d}{e}{t}=S;
            dummy=[dummy S];
        end
    meandayspectrum{d}{t}=mean(dummy,2);
    stddayspectrum{d}{t}=std(dummy,0,2);
    end
end

% Second, z-score the spectrograms.

for d=days
for e=epochs
    for t=tetrodes
        %spectrograms{d}{e}{t}=spectrograms{d}{e}{t} - repmat(meandayspectrum{d}{t},1,size(spectrograms{d}{e}{t},2));   
        %spectrograms{d}{e}{t}=spectrograms{d}{e}{t} ./ repmat(stddayspectrum{d}{t},1,size(spectrograms{d}{e}{t},2));
        spectrograms{d}{e}{t}=bsxfun(@minus,spectrograms{d}{e}{t},meandayspectrum{d}{t});   
        spectrograms{d}{e}{t}=bsxfun(@rdivide,spectrograms{d}{e}{t},stddayspectrum{d}{t});
    end
end
end

end

    
% Third, plot z-score spectra, all

for d = days
    
    for e=2
           
        epochduration = pos{d}{e}.data(end,1)-pos{d}{e}.data(1,1);
        
        for t=tetrodes
            
            if reward_flag
                eegstartstamp = eegstruct{d}{e}{t}.starttime*10000;
                outerrewardtimes = intersect(find(rewardinfo{d}{e}(:,3)==1),find(rewardinfo{d}{e}(:,4)==10));
                innerrewardtimes = intersect(find(rewardinfo{d}{e}(:,3)==1),find(rewardinfo{d}{e}(:,4)==11));
                innererrortimes = intersect(find(rewardinfo{d}{e}(:,3)==0),find(rewardinfo{d}{e}(:,4)==11));
                outererrortimes = intersect(find(rewardinfo{d}{e}(:,3)==0),find(rewardinfo{d}{e}(:,4)==10));      
                rewardinfo{d}{e}(:,5) = (rewardinfo{d}{e}(:,5) - eegstartstamp)/10000;      % time of entry (beam breaking)
            end
            
            starttime = 0;
            % plot whole epoch, one duration window at a time
            while starttime < epochduration

                endtime = starttime + winduration;

                
            %full screen
            figure('units','normalized','outerposition',[0 0 1 1])
            set(gca,'color','none')
            
            % pos
            subplot(6,1,1)
            postimes = pos{d}{e}.data(:,1)-pos{d}{e}.data(1,1);
            firstindex = lookup(starttime,postimes);
            lastindex = lookup(endtime,postimes);
            plot(postimes(firstindex:lastindex),pos{d}{e}.data(firstindex:lastindex,5),'k','LineWidth',3);
            set(gca,'color','none')
            axis tight
            ylim([0 30])
            
            % raw eeg
            subplot(6,1,2)
            eegtimevec = geteegtimes(eegstruct{d}{e}{t});
                eegtimevec = eegtimevec-eegtimevec(1);
            firstindex_eeg = lookup(starttime,eegtimevec);
            lastindex_eeg = lookup(endtime,eegtimevec);
            plot(eegtimevec(firstindex_eeg:lastindex_eeg),eegstruct{d}{e}{t}.data(firstindex_eeg:lastindex_eeg),'k','LineWidth',1);
            set(gca,'color','none')
            axis tight
            ylim([-450 450])
            
            % spectrogram
            subplot(6,1,3:6)
            firstindex_spec = lookup(starttime,timevec{d}{e}{t});
            lastindex_spec = lookup(endtime,timevec{d}{e}{t});
            imagesc(timevec{d}{e}{t}(firstindex_spec:lastindex_spec), ...
                    freqs, ...
                    spectrograms{d}{e}{t}(:,firstindex_spec:lastindex_spec),[0,2.5]);
            colormap hot
            set(gca,'YDir','normal');
            string = sprintf('%d',t);
            title(string);
            
            
            % reward plot
            if reward_flag
                
                eventrows = intersect(find(starttime < rewardinfo{d}{e}(:,5)),find(endtime > rewardinfo{d}{e}(:,5)));
                
                % add to pos plot
                subplot(6,1,1)
                hold on
                for r = eventrows
                    
                    timeofreward = rewardinfo{d}{e}(r,5);
                    
                    % reward
                    if rewardinfo{d}{e}(r,3)==1
                        % outbound
                        if rewardinfo{d}{e}(r,4)==10
                            plot([timeofreward timeofreward],[0 30],'y-','LineWidth',4)
                        % inbound
                        elseif rewardinfo{d}{e}(r,4)==11
                            plot([timeofreward timeofreward],[0 30],'y-','LineWidth',2)
                        end
                    % error
                    elseif rewardinfo{d}{e}(r,3)==0
                        if rewardinfo{d}{e}(r,4)==10
                            plot([timeofreward timeofreward],[0 30],'r-','LineWidth',4)
                        elseif rewardinfo{d}{e}(r,4)==11
                            plot([timeofreward timeofreward],[0 30],'r-','LineWidth',2)
                        end                        
                    end
                end
                
                % add to eeg plot
                subplot(6,1,2)
                hold on
                for r = eventrows
                    
                    timeofreward = rewardinfo{d}{e}(r,5);
                    
                    % reward
                    if rewardinfo{d}{e}(r,3)==1
                        % outbound
                        if rewardinfo{d}{e}(r,4)==10
                            plot([timeofreward timeofreward],[-450 450],'y-','LineWidth',4)
                        % inbound
                        elseif rewardinfo{d}{e}(r,4)==11
                            plot([timeofreward timeofreward],[-450 450],'y-','LineWidth',2)
                        end
                    % error
                    elseif rewardinfo{d}{e}(r,3)==0
                        if rewardinfo{d}{e}(r,4)==10
                            plot([timeofreward timeofreward],[-450 450],'r-','LineWidth',4)
                        elseif rewardinfo{d}{e}(r,4)==11
                            plot([timeofreward timeofreward],[-450 450],'r-','LineWidth',2)
                        end                        
                    end
                end
                
                % add to spectrogram plot
                subplot(6,1,3:6)
                hold on
                for r = eventrows
                    
                    timeofreward = rewardinfo{d}{e}(r,5);
                    
                    % reward
                    if rewardinfo{d}{e}(r,3)==1
                        % outbound
                        if rewardinfo{d}{e}(r,4)==10
                            plot([timeofreward timeofreward],[0 30],'y-','LineWidth',4)
                        % inbound
                        elseif rewardinfo{d}{e}(r,4)==11
                            plot([timeofreward timeofreward],[0 30],'y-','LineWidth',2)
                        end
                    % error
                    elseif rewardinfo{d}{e}(r,3)==0
                        if rewardinfo{d}{e}(r,4)==10
                            plot([timeofreward timeofreward],[0 30],'r-','LineWidth',4)
                        elseif rewardinfo{d}{e}(r,4)==11
                            plot([timeofreward timeofreward],[0 30],'r-','LineWidth',2)
                        end                        
                    end
                end
                
            end
            
            starttime = starttime + winduration;

            pause
            close
            
            end
        end
    end
end







        
        
        
        
        
        
        
        
        
        
        
