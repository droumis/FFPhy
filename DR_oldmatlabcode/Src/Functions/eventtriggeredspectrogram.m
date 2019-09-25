function eventtriggeredspectrogram(directoryname,animalname,Block,chan)

% based on Maggie's riptriggeredspectrogram
% out = riptriggeredspectrogram(directoryname,fileprefix,days,varargin)
%  Computes and saves the spectrogram around each ripple.
%
%directoryname - example '/data99/user/animaldatafolder/', a folder 
%                containing processed matlab data for the animal
%
%animalname -    animal specific eg 'CML9'
%
%block -         data block eg 21
% chan -        channel to be analysed eg 8
%

%  Options:
%       fpass-- Determines the frequency range for computing spectrum.
%           Default: [2 350]
%       average_trials-- Determines if events are averaged or not.
%           Default: 0
%       spectrum_window-- Determines the sliding window used to compute
%           the event triggered spectrogram. Default: [0.1 0.01]
%       event_window--Determines the size of the window around each
%           triggering event. Default: [0.2 0.4]
%       nonreference-- Specifies whether to use EEG or EEGnonreference to
%           complete filtering. Default: 1
%       cellfilter--Determines which tetrodes are used to detect triggering
%           events. Default: 'isequal($area,''CA1'')'
%       tetfilter--Determines which tetrodes ripple triggered spectrum is
%           computed for. Default: 'isequal($area,''CA1'')|isequal($area,''CA3'')'

%parse the options
params = {};
params.Fs = 24414.1;
params.fpass = [150 350];
params.trialave = 1;
params.tapers = [];
win = [2 12];
cwin = [0.1 0.01];
daytetlist = [];
cellfilter = 'isequal($area,''CA1'')';
tetfilter = 'isequal($area,''CA1'')|isequal($area,''CA3'')';
nonreference = 1;
savedirectoryname = directoryname;

dMaxV = 5e-3; %from Serge's original file - for scaling purposes
iResolutionBits = 16;
dVoltsToIntScale = pow2(iResolutionBits-1)/dMaxV;

cd(directoryname);

for i=chan;

 % load datafile
    fn = sprintf('/data14/jai/%s/Block-%d/Jai_Block-%d_xWav_ch%d.sev',...
            animalname,Block,Block,i);
        fmt = 'float32';
        fid = fopen(fn, 'rb');

        e = fread(fid, 1e9, fmt);
        e = e * dVoltsToIntScale;
    fclose(fid);

% load stimulation time file
timefile=sprintf('/data14/jai/%s/LaserStim_Params_Jai_Blk%d.mat',...
            animalname,Block);
load(timefile);


        
% define trigger times
% set stim amp
stimamp=0.5;
% set stim duration
stimdur=10000;
% find stim times matching stimamp
stamp=find((round(Data.StimAmp*1000)/1000)==stimamp);
% find stim times matching stim duration

rDur=round(Data.PulseDur*1000)/1000;
stdur=find(rDur(stamp,1)==stimdur);


% get trigger times
timelist=load(sprintf('/data14/jai/%s/timelist.mat',animalname));
timelist=timelist.t;

for u=1:size(timelist,1);
    subplot(size(timelist,1),1,u);
    triggers = timelist(u,3);  
%triggers = Data.Pulse(stamp(stdur,1),1);

endtime = (length(e)-1) * (1 / params.Fs);

    %Remove triggering events that are too close to the beginning or end
    %while triggers(1)<win(1)
    %    triggers(1) = [];
    %end
    %while triggers(end)> endtime-win(2)
    %    triggers(end) = [];
    %end
    
    
    % Calculate the event triggered spectrogram
    [S,t,f] = mtspecgramtrigc(e,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);
    
    baseline = repmat(mean(S(1:6,:)),length(t),1);
   imagesc(t-win(1),f,S'./baseline'); 
   imagesc(t-win(1),f,S'); 
   %set(gca,'xtick',[]);
   axis xy;
   ylabel(sprintf('%sV',num2str(timelist(u,1),1)),'Rotation',0);
   
   
   
   
   
   hold on;
   %title(sprintf('chan %s',num2str(i)));

    
    % Compute a z-scored spectrogram using the mean and std for the entire session
%     P = mtspecgramc(e,[cwin(1) cwin(1)],params);
%     meanP = mean(P);
%     stdP = std(P);
% 
%     for i = 1:size(S,1)
%         for j = 1:size(S,3)
%             S(i,:,j) = (S(i,:,j) - meanP)./stdP;
%         end
%     end
%    
u=u+1;
end
set(gca,'xtickMode','auto');

filename=sprintf('%s_%s_%s.pdf',animalname,num2str(Block),num2str(i));
saveas(gcf,filename,'pdf');

close;

end
end
