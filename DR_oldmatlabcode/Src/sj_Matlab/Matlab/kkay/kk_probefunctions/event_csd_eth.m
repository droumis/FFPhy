


%% THREE ingredients:
  %%  1. eeg struct     (all epochs for a day, use Mattias' loadeegstruct)
  %%  2. <event> struct  (after <event>dayprocess and <event>extract)
  %%  3. sitemap     (NOT rawpos)

%% Set these parameters manually.
animalprefix = 'eth';                              % to label graphs later
eventtype = 'ripples';                              % must be the same name as the events' data structure
days = 14;                                           % days to analyze
ref_tetrode = 2;                                   % reference tetrode for events
reward_bits = 6:8;
epochs = [1:2];                                         % epochs to analyze

sitemap=[16:32]';
channels=1:31;

velocity_state_1 = 2;                               % state 1 (immobile)
velocity_state_2 = 8;                               % state 2 (run)
windowsize_sec = .2;                               %        

eegstruct=loadeegstruct('/data13/anna/Ann/',animalprefix,'eeg',days,epochs,channels);

        % rows correspond to shanks
        % within row, left:right::dorsal:ventral

% calculated for you
nosites=max(max(sitemap));
noshanks=size(sitemap,1);
halfwindow_samp = (windowsize_sec*1500)/2;
epno=epochs(end);
filename='';
rewardflag=0;

switch eventtype
    case 'reward'
        events=DIO;
        rewardflag=1;
    case 'errortrials'
        events=errortimes;
        rewardflag=2;
    case 'rewardaudio'
        events=rewardtimes;
        rewardflag=3;
    case 'ripples'
        events=ripples;
    case 'gammal'
        events=gammal;
    case 'gammah'
        events=gammah;
    case 'gammaff'
        events=gammaff;
    case 'slowripples'
        events=slowripples;
    case 'hightheta'
        events=hightheta;
    otherwise
        disp('this event type is not recognized..'); 
end

%% First, collects all event timestamps from chosen ref. tetrode.

eventtimes=cell(days(end),epno);

if rewardflag == 1         %%% because reward is different than EEG events
    for d=days
        for e=epochs
            for w=reward_bits
                eventtimes{d,e}=[eventtimes{d,e}; events{d}{e}{w}.pulsetimes(:,1)/10000];
            end
            eventtimes{d,e}=sort(eventtimes{d,e});
        end
    end
elseif rewardflag==2          %error trials, well was triggered but no reward delivered
    for d=days
        for e=epochs
            eventtimes{d,e}=[eventtimes{d,e}; events{d,e}(:,3)/10000];
        end
    end
elseif rewardflag==3
    for d=days
        for e=epochs
            eventtimes{d,e}=[eventtimes{d,e}; events{d,e}(:,3)/10000];
        end
    end
else                    %% eeg events like ripples, gamma, etc.
    for d=days
        for e=epochs
            eventtimes{d,e}=[eventtimes{d,e}; events{d}{e}{ref_tetrode}.midtime];
        end
    end
end

%% Fourth, smooth spatially, then calculate all individual CSDs.

smootheeg=nestcell(days(end),epno,nosites);
csd=nestcell(days(end),epno,nosites);           % holds all events' spectograms

for d=days
for e=epochs
    for h=1:size(sitemap,2)
       % triangular smoothing (Freeman-Nicholson-1974) for internal sites
       for z=2:(length(sitemap(:,h))-1)   % internal sites
           mid=eegstruct{d}{e}{sitemap(z,h)}.data;
           top=eegstruct{d}{e}{sitemap(z-1,h)}.data;
           bot=eegstruct{d}{e}{sitemap(z+1,h)}.data;
           smootheeg{d}{e}{sitemap(z,h)}.data=(top+2*mid+bot)/4;   % triangle
       end
       % take second spatial derivative
       for z=3:(length(sitemap(:,h))-2)
          mid=smootheeg{d}{e}{sitemap(z,h)}.data;
          top=smootheeg{d}{e}{sitemap(z-1,h)}.data;
          bot=smootheeg{d}{e}{sitemap(z+1,h)}.data;
          csd{d}{e}{sitemap(z,h)}.data=(top+bot-2*mid);     % given 100 um site spacing and uV / um^2
       end
    end
end
end


%% Second, copy the csd / smootheeg data into windows aligned to the event starttime.

csdwin=nestcell(days(end),epno,nosites);
smootheegwin=nestcell(days(end),epno,nosites);

for d=days
    for e=epochs
        for h=1:size(sitemap,2)
            for z=1:size(sitemap,1)
                s=sitemap(z,h);
                epoch_starttime=eegstruct{d}{e}{s}.starttime;          % time, in seconds, since Open All Files
                Fs=eegstruct{d}{e}{s}.samprate;
                epoch_times=epoch_starttime:1/Fs:(epoch_starttime+(1/Fs)*length(eegstruct{d}{e}{s}.data));
                for r=2:(length(eventtimes{d,e})-1)                                   % iterate through events
                    centerindex=lookup(eventtimes{d,e}(r),epoch_times);
                    if ~isempty(csd{d}{e}{s})
                        csdwin{d}{e}{s}=cat(1,csdwin{d}{e}{s}, ...
                            csd{d}{e}{s}.data((centerindex-halfwindow_samp):(centerindex+halfwindow_samp))');
                    end
                    if ~isempty(smootheeg{d}{e}{s})
                        smootheegwin{d}{e}{s}=cat(1,smootheegwin{d}{e}{s}, ...
                            smootheeg{d}{e}{s}.data((centerindex-halfwindow_samp):(centerindex+halfwindow_samp))');
                    end
                end
            end
        end
    end
end




%% Calculate average event smootheegwin / CSD.

averagecsd=nestcell(days(end),epno,nosites);
averageeeg=nestcell(days(end),epno,nosites);

for d=days
    for e=epochs
        for h=1:size(sitemap,2)
            for z=1:size(sitemap,1)
                c=sitemap(z,h);
                if c~=0 && ~isempty(csd{d}{e}{c})
                    averagecsd{d}{e}{c}=mean(csdwin{d}{e}{c},1);
                end
                if c~=0 && ~isempty(csd{d}{e}{c})
                    averageeeg{d}{e}{c}=mean(smootheegwin{d}{e}{c},1);
                end
            end
        end
    end
end

%% Construct CSD 2D (transverse hippocampal matrix), fit polynomial to averaged CSD, then construct fine 2D CSD.

csdmatrix=nestcell(days(end),epno,size(sitemap,2));
csdpoly=nestcell(days(end),epno,size(sitemap,2));
csdfine=nestcell(days(end),epno,size(sitemap,2));

polyN=10;
zpixels=200;


for d=days
    for e=epochs
        for h=1:size(sitemap,2)    % construct a 2D matrix for each shank
            
            % initialize matrix
            dummy=NaN;
            while isnan(dummy)
                for z=1:size(sitemap,1)
                    if ~isempty(averagecsd{d}{e}{sitemap(z,h)})
                        dummy=z;
                    end
                end
            end
            tracelength=length(averagecsd{d}{e}{sitemap(dummy,h)});
            csdmatrix{d}{e}{h}=nan(length(sitemap(:,h)),tracelength);
            
            % copy into matrix
            for z=1:length(sitemap(:,h))
                c=sitemap(z,h);
                if ~isnan(averagecsd{d}{e}{c})
                    csdmatrix{d}{e}{h}(z,:)=averagecsd{d}{e}{c};
                end
            end
            
            % fit polynomial
            novertsites=min(sum(~isnan(csdmatrix{d}{e}{h}),1));
            for i=1:tracelength
                csdpoly{d}{e}{h}(:,i)=polyfit(1:novertsites, ...
                                                        csdmatrix{d}{e}{h}((1+2):(end-2),i)',...
                                                        5);
            end
            
            % construct fine 2D CSD by evaluating polynomial
            Z=linspace(1,novertsites,zpixels);
            for i=1:tracelength
                csdfine{d}{e}{h}(:,i) = polyval(csdpoly{d}{e}{h}(:,i),Z);
            end           
            
            
        end
    end
end




%% Plot csdfine.

d=14;
e=1;

figure
hold on

imagesc(csdfine{d}{e}{1})
                    set(gca,'YDir','reverse');

axis tight




%% Plot average CSD / EEG site traces.

d=14;
e=1;

yshift=350;
xshift=halfwindow_samp*2.4;
figure
hold on

for s=1:size(sitemap,2)
    for z=1:size(sitemap,1)
        c=sitemap(z,s);
        if c~=0 && ~isempty(averagecsd{d}{e}{c})
            plot((0:(halfwindow_samp*2))/Fs + xshift*(s-1), ...
                averageeeg{d}{e}{c}-z*yshift,'k','LineWidth',2)
        end
    end
end

axis tight





